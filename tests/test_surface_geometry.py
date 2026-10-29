"""One shared contour kernel, unchanged water arrays, separate opening gates."""

import numpy as np
import pytest
from scipy.spatial import ConvexHull

from topwater.geometry.waterline import validate_opening_clearance
from topwater.geometry.waterline import waterline_segments


def box():
    vertices = np.array([[x, y, z] for x in (-1., 1.) for y in (-1., 1.) for z in (-1., 1.)])
    faces = ConvexHull(vertices).simplices.copy()
    for i, triangle in enumerate(vertices[faces]):
        if np.dot(np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0]), triangle.mean(axis=0)) < 0:
            faces[i] = faces[i, ::-1]
    cap = np.flatnonzero(np.all(vertices[faces, 2] == 1., axis=1))
    return vertices, faces, np.array([1, 5, 7, 3]), cap


def field(z=0.):
    x, y = np.meshgrid(np.arange(-3., 4.), np.arange(-3., 4.))
    vertices = np.c_[x.ravel(), y.ravel(), np.full(x.size, z)]
    index = np.arange(x.size).reshape(x.shape)
    a, b, c, d = index[:-1, :-1], index[:-1, 1:], index[1:, :-1], index[1:, 1:]
    faces = np.vstack((np.stack((a, b, d), axis=-1).reshape(-1, 3),
                       np.stack((a, d, c), axis=-1).reshape(-1, 3)))
    return vertices, faces


def test_render_and_producer_use_identical_kernel():
    from topwater.render.mesh import waterline_segments as renderer_contact
    assert renderer_contact is waterline_segments


def test_far_field_faces_leave_physical_contact_bit_identical():
    from scipy.spatial.transform import Rotation
    hull, faces, _, caps = box()
    hull = hull @ Rotation.from_euler('xyz', [.07, -.13, .2]).as_matrix().T
    water, connectivity = field()
    water[:, 2] = .08*water[:, 0]-.03*water[:, 1]
    expected = waterline_segments(hull, faces, water, connectivity, artificial_faces=caps)
    # Irrelevant triangles surround the hull on all three axes; adding them
    # must preserve both coordinates and contact ordering, without remeshing.
    distant = np.concatenate([water+100*sign*np.eye(3)[j]
                              for j in range(3) for sign in (-1., 1.)])
    distant_faces = np.concatenate([connectivity+(i+1)*len(water) for i in range(6)])
    actual = waterline_segments(hull, faces, np.vstack((water, distant)),
        np.vstack((connectivity, distant_faces)), artificial_faces=caps)
    np.testing.assert_array_equal(actual, expected)


def test_independent_waterline_preserves_full_rectangular_field_and_physical_hull():
    hull, faces, rim, caps = box()
    water, connectivity = field()
    before = [value.tobytes() for value in (hull, faces, water, connectivity)]
    validate_opening_clearance(water, connectivity, hull, faces, opening_rim=rim, artificial_faces=caps)
    segments = waterline_segments(hull, faces, water, connectivity, artificial_faces=caps)
    assert segments.shape[1:] == (2, 3)
    assert len(segments) > 0
    np.testing.assert_allclose(segments[:, :, 2], 0., atol=1e-12)
    assert np.all(np.max(abs(segments[:, :, :2]), axis=2) <= 1.+1e-12)
    assert np.all(np.isclose(np.max(abs(segments[:, :, :2]), axis=2), 1., atol=1e-12))
    np.testing.assert_allclose(segments[:, :, :2].min(axis=(0, 1)), [-1., -1.], atol=1e-12)
    np.testing.assert_allclose(segments[:, :, :2].max(axis=(0, 1)), [1., 1.], atol=1e-12)
    assert before == [value.tobytes() for value in (hull, faces, water, connectivity)]
    # Interior samples stay in the rectangular representation; the posed hull
    # excludes their physical use rather than changing connectivity or values.
    interior = np.all(abs(water[:, :2]) < 1., axis=1)
    assert np.any(interior) and len(water) == 49 and len(connectivity) == 72


def test_artificial_caps_never_define_the_waterline():
    hull, faces, _, caps = box()
    water, connectivity = field(z=1.)
    physical = np.delete(faces, caps, axis=0)
    expected = waterline_segments(hull, physical, water, connectivity)
    actual = waterline_segments(hull, faces, water, connectivity, artificial_faces=caps)
    assert expected.tobytes() == actual.tobytes()
    with pytest.raises(ValueError, match='Coplanar'):
        waterline_segments(hull, faces, water, connectivity)


@pytest.mark.parametrize('z', [1., 2.])
def test_rim_touch_or_immersion_remains_fatal_without_remeshing(z):
    hull, faces, rim, caps = box()
    water, connectivity = field(z=z)
    with pytest.raises(ValueError, match='opening rim'):
        validate_opening_clearance(water, connectivity, hull, faces, opening_rim=rim, artificial_faces=caps)


def test_mid_edge_crest_and_internal_cap_crossing_remain_fatal():
    hull, faces, rim, caps = box()
    for point, message in (([0., 1.], 'opening rim'), ([0., 0.], 'opening closure')):
        water, connectivity = field()
        water[np.all(water[:, :2] == point, axis=1), 2] = 1.2
        with pytest.raises(ValueError, match=message):
            validate_opening_clearance(water, connectivity, hull, faces, opening_rim=rim, artificial_faces=caps)


def test_physical_coplanar_contact_remains_fatal_separately_from_opening_check():
    hull, faces, rim, caps = box()
    water, connectivity = field(z=-1.)
    validate_opening_clearance(water, connectivity, hull, faces, opening_rim=rim, artificial_faces=caps)
    with pytest.raises(ValueError, match='Coplanar'):
        waterline_segments(hull, faces, water, connectivity, artificial_faces=caps)


def test_missing_cap_identity_and_uncovered_rim_fail_closed():
    hull, faces, rim, caps = box()
    water, connectivity = field()
    with pytest.raises(ValueError, match='together'):
        validate_opening_clearance(water, connectivity, hull, faces, opening_rim=rim)
    water = np.array([[-.1, -.1, 0.], [.1, -.1, 0.], [0., .1, 0.]])
    with pytest.raises(ValueError, match='coverage'):
        validate_opening_clearance(water, np.array([[0, 1, 2]]), hull, faces,
                                   opening_rim=rim, artificial_faces=caps)


def test_point_tangency_has_no_artificial_contour_segment():
    hull, faces, _, _ = box()
    water = np.array([[1., 1., 1.], [2., 1., 2.], [1., 2., 2.]])
    assert waterline_segments(hull, faces, water, np.array([[0, 1, 2]])).shape == (0, 2, 3)
