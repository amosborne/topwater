"""Requested motion accuracy and quaternion-derived attitude diagnostics."""
import numpy as np
from scipy.interpolate import CubicHermiteSpline

from topwater.frames import rotation_matrix


ANGLE_ERROR_TARGET_DEGREES = .1
POSITION_ERROR_TARGET_METERS = .01
ACCURACY_DURATION_SECONDS = 60.


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


def body_poses(times, states, sample_times, *, dense_q=None):
    """Evaluate owned RK body polynomials; legacy RK23 from saved kinematics.

    Supplied dense_q stores SciPy accepted-step RK polynomial coefficients.
    Without it, only the legacy RK23 cubic can be recovered from kinematics.
    Position derivatives are saved velocities; raw quaternion derivatives are
    (0, world_omega) * q / 2. Do not canonicalize endpoint quaternion signs or
    interpolate Euler angles. The unique Hermite cubic matches RK23 dense
    output; this adds no trajectory accuracy and permits no extrapolation.
    """
    times, states, sample_times = map(np.asarray, (times, states, sample_times))
    if (times.ndim != 1 or len(times) < 2 or states.ndim != 2
            or states.shape[0] != len(times) or states.shape[1] < 13
            or sample_times.ndim != 1 or not np.isfinite(states).all()
            or not np.isfinite(times).all() or not np.isfinite(sample_times).all()
            or np.any(np.diff(times) <= 0) or np.any(sample_times < times[0])
            or np.any(sample_times > times[-1])):
        raise ValueError('finite increasing accepted states and contained sample times required')
    if dense_q is not None:
        coefficients = np.asarray(dense_q)
        if (coefficients.ndim != 3 or coefficients.shape[:2] != (len(times)-1, states.shape[1])
                or coefficients.shape[2] not in (3, 4) or not np.isfinite(coefficients).all()):
            raise ValueError('matching finite accepted RK23/RK45 dense coefficients required')
        steps = np.clip(np.searchsorted(times, sample_times, side='right')-1, 0, len(times)-2)
        h = times[steps+1]-times[steps]
        x = (sample_times-times[steps])/h
        powers = x[:, None]**np.arange(1, coefficients.shape[2]+1)
        return states[steps, :7]+h[:, None]*np.einsum('nck,nk->nc', coefficients[steps, :7], powers)
    q, omega = states[:, 3:7], states[:, 10:13]
    qdot = .5*np.column_stack((-np.sum(omega*q[:, 1:], axis=1),
        q[:, :1]*omega+np.cross(omega, q[:, 1:])))
    rates = np.column_stack((states[:, 7:10], qdot))
    return CubicHermiteSpline(times, states[:, :7], rates, extrapolate=False)(sample_times)


def motion_errors(reference, candidate, heading_degrees):
    """Per-axis candidate-minus-reference errors in one COMMON heading frame."""
    reference, candidate = np.asarray(reference), np.asarray(candidate)
    if (reference.shape != candidate.shape or reference.ndim != 2 or reference.shape[1] != 7
            or not np.isfinite(reference).all() or not np.isfinite(candidate).all()
            or not np.isfinite(heading_degrees)):
        raise ValueError('matching finite body poses and common heading required')
    delta = candidate[:, :3]-reference[:, :3]
    angle = np.deg2rad(heading_degrees)
    translation = np.column_stack((np.cos(angle)*delta[:, 0]+np.sin(angle)*delta[:, 1],
        -np.sin(angle)*delta[:, 0]+np.cos(angle)*delta[:, 1], delta[:, 2]))
    a = np.array([attitude_degrees(q) for q in reference[:, 3:7]])
    b = np.array([attitude_degrees(q) for q in candidate[:, 3:7]])
    return translation, (b-a+180)%360-180
