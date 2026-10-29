import numpy as np
from numpy.testing import assert_allclose
import pytest

from topwater import InitialState
from topwater.geometry.transforms import body_to_world
from topwater.geometry.transforms import propagation_angle_from_bearing
from topwater.geometry.transforms import rotation_matrix
from topwater.geometry.transforms import shift_moment
from topwater.geometry.transforms import world_to_body


@pytest.mark.parametrize(("axis", "expected"), [
    (0, [[1, 0, 0], [0, 0, -1], [0, 1, 0]]),
    (1, [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]),
    (2, [[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
])
def test_right_hand_quarter_turns(axis, expected):
    q = np.zeros(4)
    q[0] = q[axis + 1] = np.sqrt(0.5)
    assert_allclose(rotation_matrix(q), expected, atol=1e-15)
    # Rows of basis vectors must rotate just as individual vectors do.
    assert_allclose(body_to_world(np.eye(3), q), np.array(expected).T, atol=1e-15)
    assert_allclose(body_to_world([1, 0, 0], q), np.array(expected)[:, 0], atol=1e-15)


def test_arbitrary_rotation_against_rodrigues():
    axis = np.array([1, 2, 3]) / np.sqrt(14)
    angle = 0.7
    q = np.r_[np.cos(angle / 2), axis * np.sin(angle / 2)]
    v = np.array([2, -3, 4])
    expected = (v * np.cos(angle) + np.cross(axis, v) * np.sin(angle)
                + axis * np.dot(axis, v) * (1 - np.cos(angle)))
    assert_allclose(body_to_world(v, q), expected, atol=1e-14)
    assert_allclose(body_to_world(v, -q), expected, atol=1e-14)
    assert_allclose(world_to_body(expected, q), v, atol=1e-14)
    matrix = rotation_matrix(q)
    assert_allclose(matrix.T @ matrix, np.eye(3), atol=1e-15)
    assert_allclose(np.linalg.det(matrix), 1, atol=1e-15)


def test_point_transform_uses_cg_not_mesh_origin():
    q = [np.sqrt(0.5), 0, 0, np.sqrt(0.5)]
    cg = np.array([2, 3, 4])
    position = np.array([10, 20, 30])
    point = np.array([3, 3, 4])  # One metre toward the bow from CG.
    assert_allclose(position + body_to_world(point - cg, q), [10, 21, 30])


def test_moment_shift_and_rotation():
    force = [0, 0, -10]
    point = [2, 0, 0]
    assert_allclose(shift_moment(force, [0, 0, 0], point, [0, 0, 0]), [0, 20, 0])
    q = [np.sqrt(0.5), 0, 0, np.sqrt(0.5)]
    expected = [-20, 0, 0]
    assert_allclose(
        shift_moment(body_to_world(force, q), [0, 0, 0], body_to_world(point, q), [0, 0, 0]),
        expected, atol=1e-14,
    )


def test_moment_shift_preserves_mechanical_power():
    force = np.array([3, -4, 5])
    moment = np.array([6, 7, -8])
    from_point = np.array([2, 3, 4])
    to_point = np.array([-1, 5, 2])
    omega = np.array([0.2, -0.3, 0.4])
    velocity_from = np.array([1, 2, 3])
    velocity_to = velocity_from + np.cross(omega, to_point - from_point)
    shifted = shift_moment(force, moment, from_point, to_point)
    assert_allclose(
        force @ velocity_from + moment @ omega,
        force @ velocity_to + shifted @ omega,
    )
    assert_allclose(shift_moment(force, shifted, to_point, from_point), moment)


@pytest.mark.parametrize(("bearing", "expected"), [
    (0, 3 * np.pi / 2), (90, np.pi), (180, np.pi / 2), (270, 0),
    (360, 3 * np.pi / 2), (-90, 0), (45, 5 * np.pi / 4),
])
def test_coming_from_compass_bearings(bearing, expected):
    assert_allclose(propagation_angle_from_bearing(bearing), expected)


@pytest.mark.parametrize("q", [[0, 0, 0, 0], [2, 0, 0, 0], [1, 0, 0], [np.nan, 0, 0, 0]])
def test_invalid_quaternion_rejected_by_state_and_frames(q):
    with pytest.raises(ValueError):
        InitialState(orientation=q)
    with pytest.raises(ValueError):
        rotation_matrix(q)


@pytest.mark.parametrize("vectors", [[1, 2], [[[1, 2, 3]]], [1, np.inf, 3]])
def test_invalid_vector_shapes_and_values(vectors):
    with pytest.raises(ValueError):
        body_to_world(vectors, [1, 0, 0, 0])


def test_nonfinite_bearing_rejected():
    with pytest.raises(ValueError, match="finite"):
        propagation_angle_from_bearing(np.nan)
