"""SI, right-handed frames: world east/north/up; body bow/port/up.

Orientation is an active body-to-world unit quaternion in (w, x, y, z) order.
Positive rotations follow the right-hand rule; no Euler-angle convention is used.
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.spatial.transform import Rotation

from .._validation import array, quaternion


def rotation_matrix(orientation: ArrayLike) -> NDArray:
    """Return R such that a column vector satisfies v_world = R @ v_body."""
    return Rotation.from_quat(quaternion(orientation), scalar_first=True).as_matrix()


def body_to_world(vectors: ArrayLike, orientation: ArrayLike) -> NDArray:
    """Rotate a vector (3,) or rows of vectors (n, 3); apply no translation."""
    values = np.asarray(vectors)
    shape = (3,) if values.ndim == 1 else (None, 3)
    return array(values, shape, "vectors") @ rotation_matrix(orientation).T


def world_to_body(vectors: ArrayLike, orientation: ArrayLike) -> NDArray:
    """Inverse of body_to_world for free vectors, not absolute positions."""
    values = np.asarray(vectors)
    shape = (3,) if values.ndim == 1 else (None, 3)
    return array(values, shape, "vectors") @ rotation_matrix(orientation)


def shift_moment(
    force: ArrayLike,
    moment: ArrayLike,
    from_point: ArrayLike,
    to_point: ArrayLike,
) -> NDArray:
    """Shift a load's moment between points, with all arguments in one frame.

    M_to = M_from + (from_point - to_point) cross F. Force is unchanged.
    Rotate force, moment and point coordinates first when changing frames.
    """
    force = array(force, (3,), "force")
    moment = array(moment, (3,), "moment")
    origin = array(from_point, (3,), "from_point")
    target = array(to_point, (3,), "to_point")
    return moment + np.cross(origin - target, force)


def propagation_angle_from_bearing(coming_from_degrees: float) -> float:
    """Convert a coming-from compass bearing to a propagation angle in radians.

    Input is clockwise from north. Output is counterclockwise from world +x
    (east), in [0, 2*pi). Thus a wave from north travels toward south (3*pi/2).
    """
    bearing = array(coming_from_degrees, (), "coming_from_degrees").item()
    return float(np.deg2rad((270.0 - bearing % 360.0) % 360.0))
