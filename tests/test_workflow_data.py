"""Portable final-data contracts and physical-only boat ingestion."""

from dataclasses import replace
from zipfile import ZipFile

import numpy as np
import pytest

from topwater.boat_geometry import physical_hull
from topwater.workflow_data import Boat, Result, Sea


def open_box():
    vertices = np.array([[-1., -1., 0.], [1., -1., 0.], [1., 1., 0.], [-1., 1., 0.],
                         [-1., -1., 1.], [1., -1., 1.], [1., 1., 1.], [-1., 1., 1.]])
    faces = np.array([[0, 2, 1], [0, 3, 2], [0, 1, 5], [0, 5, 4],
                      [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6],
                      [3, 0, 4], [3, 4, 7]])
    return vertices, faces


def boat():
    vertices, faces = open_box()
    return Boat(vertices, faces, 1000., np.array([0., 0., .5]), np.diag([400., 400., 600.]))


def result():
    vertices, faces = open_box()
    # Different cadences and different surface topology, all face indices local.
    water0 = np.array([[-3., -3., 0.], [3., -3., 0.], [0., 3., 0.]])
    water1 = np.array([[-3., -3., .1], [3., -3., .1], [3., 3., .1], [-3., 3., .1]])
    return Result(time_s=np.array([0., .5, 1.]), position_m=np.array([[2., 3., .5]]*3),
                  orientation_wxyz=np.array([[1., 0., 0., 0.]]*3),
                  hull_vertices_m=vertices, hull_faces=faces, cg_m=np.array([0., 0., .5]),
                  surface_time_s=np.array([0., 1.]), surface_vertices_m=np.vstack((water0, water1)),
                  surface_faces=np.array([[0, 1, 2], [0, 1, 2], [0, 2, 3]]),
                  surface_vertex_offsets=np.array([0, 3, 7]), surface_face_offsets=np.array([0, 1, 3]))


def test_inferred_rim_and_internal_cap_preserve_physical_faces():
    source = boat()
    body = source.to_body()
    np.testing.assert_array_equal(body.mesh.faces[body.mesh.physical_faces], source.faces)
    np.testing.assert_array_equal(body.mesh.vertices[:len(source.vertices_m)], source.vertices_m)
    assert set(body.mesh.opening_rim) == {4, 5, 6, 7}
    assert len(body.mesh.artificial_faces) == 4
    assert body.mass == source.mass_kg
    # Internal cap vertices/panels never become part of saved physical geometry.
    assert len(body.mesh.vertices) == len(source.vertices_m)+1
    assert not source.vertices_m.flags.writeable


def test_canonical_hull_with_collinear_opening_vertices():
    from tests.reference_boat import mesh_arrays
    vertices, faces, rim, cap = mesh_arrays()
    physical = np.delete(faces, cap, axis=0)
    mesh = physical_hull(vertices, physical)
    assert set(mesh.opening_rim) == set(rim)
    np.testing.assert_array_equal(mesh.faces[mesh.physical_faces], physical)


def test_nonplanar_opening_caps_missing_envelope_facets_without_changing_hull():
    vertices, faces = open_box()
    # Remove the complete front side as well as the already absent top.
    physical = np.delete(faces, [2, 3], axis=0)
    mesh = physical_hull(vertices, physical)
    np.testing.assert_array_equal(mesh.vertices, vertices)
    np.testing.assert_array_equal(mesh.faces[mesh.physical_faces], physical)
    assert set(mesh.opening_rim) == {0, 1, 4, 5, 6, 7}
    caps = mesh.vertices[mesh.faces[mesh.artificial_faces]]
    assert all(np.allclose(t[:, 2], 1.) or np.allclose(t[:, 1], -1.) for t in caps)


def test_nonplanar_cap_preserves_collinear_rim_vertices():
    vertices, faces = open_box()
    physical = np.delete(faces, [2, 3], axis=0)
    vertices = np.vstack((vertices, (vertices[1]+vertices[5])/2))
    physical = np.array([f for f in physical if not np.array_equal(f, [1, 6, 5])]
                        + [np.array([1, 6, 8]), np.array([8, 6, 5])])
    mesh = physical_hull(vertices, physical)
    assert 8 in mesh.opening_rim
    assert 8 in mesh.faces[mesh.artificial_faces]
    np.testing.assert_array_equal(mesh.faces[mesh.physical_faces], physical)


def test_partial_nonplanar_cap_has_distinct_backend_limitation():
    from topwater.boat_geometry import UnsupportedClosureError
    vertices, faces = open_box()
    with pytest.raises(UnsupportedClosureError, match='valid hull'):
        physical_hull(vertices, np.delete(faces, 3, axis=0))


@pytest.mark.parametrize('change', ['closed', 'missing_face', 'winding', 'nonconvex', 'duplicate',
                                    'degenerate', 'unused', 'nonfinite', 'floating_faces'])
def test_bad_physical_hulls_fail(change):
    vertices, faces = open_box()
    if change == 'closed':
        faces = np.vstack((faces, [[4, 5, 6], [4, 6, 7]]))
    elif change == 'missing_face':
        faces = faces[1:]
    elif change == 'winding':
        faces[0] = faces[0, ::-1]
    elif change == 'nonconvex':
        vertices[0] = [0., 0., .5]
    elif change == 'duplicate':
        faces = np.vstack((faces, faces[0]))
    elif change == 'degenerate':
        faces[0] = [0, 0, 1]
    elif change == 'unused':
        vertices = np.vstack((vertices, [0., 0., 0.]))
    elif change == 'nonfinite':
        vertices[0, 0] = np.nan
    else:
        faces = faces.astype(float)
    with pytest.raises(ValueError):
        physical_hull(vertices, faces)


def test_extra_opening_rejected():
    from scipy.spatial import ConvexHull
    # A convex octahedron with two opposite physical triangles removed.
    vertices = np.array([[1., 0., 0.], [-1., 0., 0.], [0., 1., 0.],
                         [0., -1., 0.], [0., 0., 1.], [0., 0., -1.]])
    hull = ConvexHull(vertices)
    faces = hull.simplices.copy()
    for i, triangle in enumerate(vertices[faces]):
        if np.dot(np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0]), triangle.mean(axis=0)) < 0:
            faces[i] = faces[i, ::-1]
    opposite = next(i for i in range(1, len(faces)) if not set(faces[0]) & set(faces[i]))
    with pytest.raises(ValueError, match='one opening|exactly one'):
        physical_hull(vertices, np.delete(faces, [0, opposite], axis=0))


def test_geometric_overlap_and_t_junction_are_not_valid_adjacency():
    from topwater.boat_geometry import _check_intersections
    # Same support plane and distinct indices: one triangle lies within another.
    points = np.array([[0., 0., 0.], [2., 0., 0.], [0., 2., 0.],
                       [.2, .2, 0.], [.6, .2, 0.], [.2, .6, 0.]])
    with pytest.raises(ValueError, match='intersecting'):
        _check_intersections(points, np.array([[0, 1, 2], [3, 4, 5]]),
                             np.array([[0., 0., 1.]]*2), 1e-10)
    # A legitimate shared edge is allowed, but an unindexed edge contact is not.
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
                       [1., 1., 0.], [.5, .5, 0.], [2., 0., 0.], [2., 1., 0.]])
    _check_intersections(points, np.array([[0, 1, 2], [1, 3, 2]]),
                         np.array([[0., 0., 1.]]*2), 1e-10)
    with pytest.raises(ValueError, match='T-junction'):
        _check_intersections(points, np.array([[0, 1, 2], [4, 5, 6]]),
                             np.array([[0., 0., 1.]]*2), 1e-10)


def test_scale_and_transform_do_not_change_physical_identity():
    from scipy.spatial.transform import Rotation
    vertices, faces = open_box()
    rotation = Rotation.from_rotvec([.2, -.1, .7]).as_matrix()
    for scale in (1e-6, 1., 1e6):
        moved = (vertices @ rotation.T + [8., -6., 2.])*scale
        mesh = physical_hull(moved, faces)
        np.testing.assert_array_equal(mesh.faces[mesh.physical_faces], faces)


@pytest.mark.parametrize('factory', [boat, result, lambda: Sea([.1], [2.], [.3], [.4]),
                                      lambda: Sea([], [], [], [])])
def test_numeric_archive_roundtrip(factory, tmp_path):
    original = factory()
    path = tmp_path/'nested'/'data.npz'
    original.save(path)
    loaded = type(original).load(path)
    with np.load(path, allow_pickle=False) as archive:
        assert set(archive.files) == set(original.members)
        assert all(archive[name].dtype.kind in 'iuf' for name in archive.files)
    for name in original.members:
        np.testing.assert_array_equal(getattr(loaded, name), getattr(original, name))


def test_atomic_save_replacement_policy(tmp_path):
    path = tmp_path/'sea.npz'
    first, second = Sea([.1], [1.], [0.], [0.]), Sea([.2], [2.], [1.], [.3])
    first.save(path)
    with pytest.raises(FileExistsError):
        second.save(path)
    np.testing.assert_array_equal(Sea.load(path).amplitudes_m, [.1])
    second.save(path, overwrite=True)
    np.testing.assert_array_equal(Sea.load(path).amplitudes_m, [.2])
    assert list(tmp_path.iterdir()) == [path]


def test_failed_atomic_write_preserves_existing_output(tmp_path, monkeypatch):
    import topwater.workflow_data as module
    path = tmp_path/'sea.npz'
    sea = Sea([.1], [1.], [0.], [0.])
    sea.save(path)
    original = path.read_bytes()
    def fail(*args, **kwargs):
        raise OSError('injected archive write failure')
    monkeypatch.setattr(module.np, 'savez', fail)
    with pytest.raises(OSError, match='injected'):
        sea.save(path, overwrite=True)
    assert path.read_bytes() == original
    assert list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize('change', ['missing', 'extra', 'object', 'duplicate'])
def test_bad_archive_members_rejected(change, tmp_path):
    sea = Sea([.1], [1.], [0.], [0.])
    fields = {name: getattr(sea, name) for name in sea.members}
    path = tmp_path/'bad.npz'
    if change == 'missing':
        fields.pop('phases_rad')
    elif change == 'extra':
        fields['solver'] = 'wecsim'
    elif change == 'object':
        fields['phases_rad'] = np.array([object()], dtype=object)
    np.savez(path, **fields)
    if change == 'duplicate':
        with ZipFile(path, 'a') as archive:
            with pytest.warns(UserWarning, match='Duplicate'):
                archive.writestr('phases_rad.npy', b'invalid duplicate')
    with pytest.raises(ValueError):
        Sea.load(path)


def test_packed_surface_local_indices_and_independent_times():
    data = result()
    vertices, faces = data.surface_frame(1)
    assert vertices.shape == (4, 3)
    np.testing.assert_array_equal(faces, [[0, 1, 2], [0, 2, 3]])
    assert not vertices.flags.writeable and not faces.flags.writeable
    with pytest.raises(IndexError):
        data.surface_frame(-1)


def test_repeated_connectivity_still_checks_later_frame_geometry():
    data = result()
    vertices, faces = data.surface_frame(1)
    later = vertices.copy()
    later[1] = later[0]
    with pytest.raises(ValueError, match='surface frame 1 has degenerate triangles'):
        replace(data, surface_vertices_m=np.vstack((vertices, later)),
                surface_faces=np.vstack((faces, faces)),
                surface_vertex_offsets=np.array([0, 4, 8]),
                surface_face_offsets=np.array([0, 2, 4]))


@pytest.mark.parametrize('field,value', [
    ('time_s', [0., .5, .5]),
    ('time_s', [-.1, .5, 1.]),
    ('surface_time_s', [0., 2.]),
    ('surface_vertex_offsets', [0, 3, 6]),
    ('surface_vertex_offsets', [0, 0, 7]),
    ('surface_vertex_offsets', [0., 3., 7.]),
    ('surface_face_offsets', [0, 1, 2]),
    ('surface_faces', [[0, 1, 2], [0, 1, 2], [0, 2, 4]]),
    ('surface_faces', [[0, 1, 2], [0, 1, 2], [0, 1, 2]]),
    ('surface_faces', [[0, 0, 2], [0, 1, 2], [0, 2, 3]]),
    ('orientation_wxyz', [[1., 0., 0., 0.], [2., 0., 0., 0.], [1., 0., 0., 0.]]),
    ('position_m', [[0., 0., 0.], [0., np.nan, 0.], [0., 0., 0.]]),
])
def test_bad_result_rejected(field, value):
    with pytest.raises(ValueError):
        replace(result(), **{field: value})
