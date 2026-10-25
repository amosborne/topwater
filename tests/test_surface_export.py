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
