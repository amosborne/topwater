"""Convex exclusion and manifold publication checks, independent of solvers."""

import numpy as np
import pytest
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation

from topwater._surface_export import export_surface, validate_surface_topology


def closed_box():
    points = np.array([[x, y, z] for x in (-1., 1.) for y in (-1., 1.) for z in (-1., 1.)])
    faces = ConvexHull(points).simplices.copy()
    for i, triangle in enumerate(points[faces]):
        if np.dot(np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0]), triangle.mean(axis=0)) < 0:
            faces[i] = faces[i, ::-1]
    return points, faces


def patch(z=0., subdivisions=1):
    x, y = np.meshgrid(np.linspace(-3., 3., subdivisions+1), np.linspace(-3., 3., subdivisions+1))
    vertices = np.column_stack((x.ravel(), y.ravel(), np.full(x.size, z)))
    index = np.arange(x.size).reshape(x.shape)
    a, b, c, d = index[:-1, :-1], index[:-1, 1:], index[1:, :-1], index[1:, 1:]
    faces = np.vstack((np.stack((a, b, d), axis=-1).reshape(-1, 3),
                       np.stack((a, d, c), axis=-1).reshape(-1, 3)))
    return vertices, faces


def area(vertices, faces):
    triangles = vertices[faces]
    return np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0],
                                   triangles[:, 2]-triangles[:, 0]), axis=1).sum()/2


def boundary_edges(faces):
    edges = np.sort(np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]])), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    return unique[counts == 1]


@pytest.mark.parametrize('subdivisions', [1, 4, 7])
def test_flat_water_hole_area_and_shared_actual_contact(subdivisions):
    hull, hull_faces = closed_box()
    source, faces = patch(subdivisions=subdivisions)
    original = source.copy()
    water, triangles = export_surface(source, faces, hull, hull_faces)
    np.testing.assert_array_equal(source, original)
    np.testing.assert_allclose(area(water, triangles), 32., atol=2e-9, rtol=0)
    np.testing.assert_array_equal(water[:, 2], np.zeros(len(water)))
    validate_surface_topology(water, triangles)
    edges = boundary_edges(triangles)
    boundary = water[edges]
    outer = np.any(np.all(np.isclose(abs(boundary[:, :, :2]), 3., atol=1e-9), axis=1), axis=1)
    contact = boundary[~outer]
    assert len(contact) >= 4
    # Each contact edge lies on a real box side, not just the AABB exterior.
    for edge in contact:
        assert any(np.allclose(edge[:, axis], side, atol=1e-9, rtol=0)
                   for axis in (0, 1) for side in (-1., 1.))
        assert np.max(abs(edge[:, :2])) <= 1.+1e-9
    centers = water[triangles].mean(axis=1)
    assert not np.any(np.all(abs(centers[:, :2]) < 1.-1e-9, axis=1))


def test_pose_and_scale_preserve_exclusion_geometry():
    hull, hull_faces = closed_box()
    source, faces = patch(subdivisions=4)
    rotation = Rotation.from_rotvec([.2, -.3, .7]).as_matrix()
    for scale in (1e-3, 1e3):
        shift = np.array([5., -2., 1.])
        posed_hull = (hull @ rotation.T+shift)*scale
        posed_surface = (source @ rotation.T+shift)*scale
        water, triangles = export_surface(posed_surface, faces, posed_hull, hull_faces)
        local = (water/scale-shift) @ rotation
        np.testing.assert_allclose(local[:, 2], 0., atol=1e-9)
        np.testing.assert_allclose(area(water, triangles), 32.*scale**2, rtol=1e-9)


@pytest.mark.parametrize('subdivisions,pose', [(1, False), (4, False), (7, True)])
def test_registry_reuse_is_bitexact_to_per_polygon_materialization(monkeypatch, subdivisions, pose):
    """Compare complete contact meshes against the former allocation path."""
    import topwater.fluid_boundary as boundary

    hull, hull_faces = closed_box()
    source, faces = patch(subdivisions=subdivisions)
    if pose:
        rotation = Rotation.from_rotvec([.2, -.3, .7]).as_matrix()
        hull = (hull @ rotation.T+[5., -2., 1.])*1e-3
        source = (source @ rotation.T+[5., -2., 1.])*1e-3
    optimized_vertices, optimized_faces = export_surface(source, faces, hull, hull_faces)
    ears = boundary._ears
    calls = []

    def previous_materialization(indices, vertices, normal):
        # Before the optimization, each call converted a list of the complete
        # registry's row arrays to a fresh array. Recreate that same sequence.
        rows = [row for row in vertices]
        materialized = np.asarray(rows)
        np.testing.assert_array_equal(materialized, vertices)
        calls.append(len(indices))
        return ears(indices, materialized, normal)

    monkeypatch.setattr(boundary, '_ears', previous_materialization)
    previous_vertices, previous_faces = export_surface(source, faces, hull, hull_faces)
    assert len(calls) > 1  # Actual fragmented contact, not a no-clipping shortcut.
    assert optimized_vertices.tobytes() == previous_vertices.tobytes()
    assert optimized_faces.tobytes() == previous_faces.tobytes()


@pytest.mark.parametrize('transformed', [None, (.2, 1e-3), (-.3, 1e3)])
def test_contact_edge_vertex_retained_when_projected_ear_leaves_degenerate_sliver(transformed):
    """Reduced actual frame-277 polygon: repair its diagonal, never its rim."""
    from topwater._surface_export import _export_ears
    from topwater.fluid_boundary import _ears

    # Indices retain the original order relation B2406 < A5940 < C5941 < D6248.
    # D is a shared cut vertex on the constrained C-to-A contact segment.
    points = np.array([
        [.42790402925610144, .17039693147594254, -.03260405850092889],
        [.42797449851073016, .16844460462755206, -.03263331227243573],
        [.42600476914629437, .16835542952517174, -.03262962592722955],
        [.42699834260851194, .16840041134694653, -.032631485398299825],
    ])
    world = np.array([
        [2.760099464395641, 1.3193845747326938, -.36807457150284756],
        [2.760618817341613, 1.304996077145282, -.36829016953198696],
        [2.746102064559463, 1.304338863550891, -.36826300145347124],
        [2.7534246239842948, 1.3046703760917433, -.36827670561116926],
    ])
    normal = np.array([.00254541746973167, -.014890467616889126, .9998858909015844])
    polygon = [1, 0, 2, 3]
    if transformed is None:
        previous = np.asarray(_ears(polygon, points, normal))
        with pytest.raises(ValueError, match='degenerate triangles'):
            validate_surface_topology(world, previous)
    else:
        angle, scale = transformed
        rotation = Rotation.from_rotvec([angle, -.1, .7]).as_matrix()
        points = points @ rotation.T
        world = (world @ rotation.T+[5., -2., 1.])*scale
        normal = rotation @ normal
    faces = np.asarray(_export_ears(polygon, points, normal, world))
    validate_surface_topology(world, faces)
    assert set(faces.ravel()) == set(polygon)
    expected_boundary = {tuple(sorted((a, b))) for a, b in zip(polygon, np.roll(polygon, -1))}
    assert {tuple(edge) for edge in boundary_edges(faces)} == expected_boundary
    # Every constrained edge, including C-D and D-A, remains. The formerly
    # unconstrained A-C diagonal is replaced when that triangle is degenerate.
    triangle = world[faces]
    cross = np.cross(triangle[:, 1]-triangle[:, 0], triangle[:, 2]-triangle[:, 0])
    assert np.all(cross @ normal > 0)
    q = world[polygon]-world[polygon[0]]
    polygon_vector = np.cross(q, np.roll(q, -1, axis=0)).sum(axis=0)
    assert np.linalg.norm(cross.sum(axis=0)-polygon_vector) <= 256*np.finfo(float).eps*np.linalg.norm(cross, axis=1).sum()


def test_projected_exception_does_not_accept_reversed_concave_quad_fan(monkeypatch):
    from topwater._surface_export import _export_ears
    import topwater.fluid_boundary as boundary

    # Both fan triangles have nonzero area, but the first reverses orientation
    # and their unsigned area exceeds the concave parent's area.
    points = np.array([[0., 0., 0.], [.5, .5, 0.], [2., 0., 0.], [0., 2., 0.]])
    def projected_failure(*args):
        raise ValueError('degenerate final wet-facet triangle')
    monkeypatch.setattr(boundary, '_ears', projected_failure)
    with pytest.raises(ValueError, match='oriented parent coverage'):
        _export_ears([0, 1, 2, 3], points, np.array([0., 0., 1.]), points)


def test_above_hull_and_nonarea_tangency_keep_surface():
    hull, hull_faces = closed_box()
    source, faces = patch(z=2.)
    water, triangles = export_surface(source, faces, hull, hull_faces)
    np.testing.assert_array_equal(water, source)
    np.testing.assert_array_equal(triangles, faces)
    source = np.array([[1., -.5, 0.], [3., -.5, 0.], [1., .5, 0.]])
    water, triangles = export_surface(source, np.array([[0, 1, 2]]), hull, hull_faces)
    np.testing.assert_allclose(area(water, triangles), 1., atol=1e-10)


def test_coplanar_positive_area_contact_rejected():
    hull, hull_faces = closed_box()
    source, faces = patch(z=1.)
    with pytest.raises(ValueError, match='coplanar'):
        export_surface(source, faces, hull, hull_faces)


def test_surface_entirely_inside_hull_rejected():
    hull, hull_faces = closed_box()
    source = np.array([[-.2, -.2, 0.], [.2, -.2, 0.], [0., .2, 0.]])
    with pytest.raises(ValueError, match='no water triangles'):
        export_surface(source, np.array([[0, 1, 2]]), hull, hull_faces)


def test_water_winding_preserved_instead_of_silently_repaired():
    hull, hull_faces = closed_box()
    source, faces = patch()
    water, triangles = export_surface(source, faces[:, ::-1], hull, hull_faces)
    t = water[triangles]
    assert np.all(np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])[:, 2] < 0)


def test_opening_must_have_unambiguous_dry_rim_coverage():
    hull, hull_faces = closed_box()
    rim = np.array([1, 5, 7, 3])
    cap = np.flatnonzero(np.all(hull[hull_faces, 2] == 1., axis=1))
    opening = dict(opening_rim=rim, artificial_faces=cap)
    source, faces = patch(z=0.)
    export_surface(source, faces, hull, hull_faces, **opening)
    for z in (1., 2.):
        source, faces = patch(z=z)
        with pytest.raises(ValueError, match='opening rim'):
            export_surface(source, faces, hull, hull_faces, **opening)
    source = np.array([[-.1, -.1, 0.], [.1, -.1, 0.], [0., .1, 0.]])
    with pytest.raises(ValueError, match='coverage'):
        export_surface(source, np.array([[0, 1, 2]]), hull, hull_faces, **opening)


def test_rim_edge_crest_and_interior_cap_crossing_rejected():
    hull, hull_faces = closed_box()
    cap = np.flatnonzero(np.all(hull[hull_faces, 2] == 1., axis=1))
    opening = dict(opening_rim=np.array([1, 5, 7, 3]), artificial_faces=cap)
    source, faces = patch(subdivisions=6)
    # A crest on a rim-edge midpoint; the rim's four vertices remain dry.
    source[np.all(source[:, :2] == [0., 1.], axis=1), 2] = 1.2
    with pytest.raises(ValueError, match='opening rim'):
        export_surface(source, faces, hull, hull_faces, **opening)
    source, faces = patch(subdivisions=6)
    # Rim is dry, but the represented field crosses the internal lid inside it.
    source[np.all(source[:, :2] == [0., 0.], axis=1), 2] = 1.2
    with pytest.raises(ValueError, match='opening closure'):
        export_surface(source, faces, hull, hull_faces, **opening)


def test_topology_allows_disconnected_components_and_vertical_surface():
    first = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    vertical = np.array([[2., 0., 0.], [2., 1., 0.], [2., 0., 1.]])
    validate_surface_topology(np.vstack((first, vertical)), np.array([[0, 1, 2], [3, 4, 5]]))


def test_vertex_pinches_and_shared_edge_winding_rejected():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
                       [-1., 0., 0.], [0., -1., 0.]])
    with pytest.raises(ValueError, match='vertex|pinched'):
        validate_surface_topology(points, np.array([[0, 1, 2], [0, 3, 4]]))
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]])
    with pytest.raises(ValueError, match='winding'):
        validate_surface_topology(points, np.array([[0, 1, 2], [1, 2, 3]]))


def test_nonmanifold_edge_rejected():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., -1., 0.], [.5, 0., 1.]])
    with pytest.raises(ValueError, match='non-manifold edge'):
        validate_surface_topology(points, np.array([[0, 1, 2], [1, 0, 3], [0, 1, 4]]))
