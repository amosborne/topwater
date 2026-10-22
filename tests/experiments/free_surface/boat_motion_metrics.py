"""Requested motion accuracy and quaternion-derived attitude diagnostics."""
import numpy as np

from topwater.frames import rotation_matrix


ANGLE_ERROR_TARGET_DEGREES = .1
POSITION_ERROR_TARGET_METERS = .01


def attitude_degrees(orientation):
    """Return diagnostic roll/pitch/yaw; positive roll means port rises.

    Roll is world vertical in the body port/up plane. Yaw is bow heading
    counterclockwise from world east; positive ZYX pitch means bow down in
    the bow/port/up convention. These are diagnostic angles only;
    quaternion state evolution remains unchanged. Normalize a copy as for
    physical geometry reconstruction. Compare periodic angles modulo360.
    """
    q = np.asarray(orientation, dtype=float)
    if q.shape != (4,) or not np.isfinite(q).all() or np.linalg.norm(q) == 0:
        raise ValueError('finite nonzero wxyz quaternion required')
    matrix = rotation_matrix(q/np.linalg.norm(q))
    return np.degrees([np.arctan2(matrix[2, 1], matrix[2, 2]),
        np.arctan2(-matrix[2, 0], np.hypot(matrix[2, 1], matrix[2, 2])),
        np.arctan2(matrix[1, 0], matrix[0, 0])])
