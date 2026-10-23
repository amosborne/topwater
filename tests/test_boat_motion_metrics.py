"""Roll reporting follows the same active body-to-world quaternion convention."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation
from scipy.integrate import RK23

from .experiments.free_surface.boat_motion_metrics import attitude_degrees, body_poses, motion_errors


@pytest.mark.parametrize('roll', [-30., -20., -.1, 0., .1, 20., 30.])
def test_roll_is_independent_of_heading_pitch_quaternion_sign_and_norm(roll):
    rotation = (Rotation.from_rotvec([0., 0., 1.7])
        *Rotation.from_rotvec([0., -.2, 0.])
        *Rotation.from_rotvec([np.deg2rad(roll), 0., 0.]))
    q = rotation.as_quat(scalar_first=True)
    assert_allclose(attitude_degrees(q), [roll, np.rad2deg(-.2), np.rad2deg(1.7)], atol=2e-14)
    assert_allclose(attitude_degrees(-2*q), attitude_degrees(q), atol=2e-14)


def test_heading_wrap_does_not_change_small_attitude_error():
    a = attitude_degrees(Rotation.from_rotvec([0., 0., np.deg2rad(179.98)]).as_quat(scalar_first=True))
    b = attitude_degrees(Rotation.from_rotvec([0., 0., np.deg2rad(-179.98)]).as_quat(scalar_first=True))
    assert_allclose((b-a+180)%360-180, [0., 0., .04], atol=2e-14)


def test_saved_rk23_dense_coefficients_recover_accepted_body_output():
    def rhs(t, y):
        q, omega = y[3:7], y[10:13]
        qdot = .5*np.r_[-omega@q[1:], q[0]*omega+np.cross(omega, q[1:])]
        return np.r_[y[7:10], qdot, [np.cos(t), -.3, .1*t], [.4, .2*t, -.3]]
    initial = np.r_[[.1, -.2, .3], [1., 0., 0., 0.], [.2, .3, -.1], [.1, -.2, .4]]
    solver = RK23(rhs, 0., initial, .2, first_step=.02, max_step=.04, rtol=1e-4, atol=1e-7)
    times, states, coefficients, samples, expected = [0.], [initial.copy()], [], [], []
    while solver.status == 'running':
        solver.step()
        interior = np.linspace(solver.t_old, solver.t, 5)
        samples.extend(interior)
        dense = solver.dense_output()
        expected.extend(dense(interior)[:7].T)
        coefficients.append(dense.Q.copy())
        times.append(solver.t); states.append(solver.y.copy())
    assert_allclose(body_poses(times, states, samples, dense_q=coefficients), expected, rtol=0, atol=3e-16)
    assert_allclose(body_poses(times, states, times, dense_q=coefficients), np.array(states)[:, :7], rtol=0, atol=3e-16)
    with pytest.raises(ValueError, match='contained'):
        body_poses(times, states, [-1e-9], dense_q=coefficients)


def test_motion_error_uses_one_declared_heading_frame():
    q = Rotation.from_rotvec([0., 0., np.pi/2]).as_quat(scalar_first=True)
    reference = np.r_[np.zeros(3), q][None]
    candidate = reference.copy()
    candidate[0, :3] = [-.003, .005, .002]
    translation, angles = motion_errors(reference, candidate, 90.)
    assert_allclose(translation, [[.005, .003, .002]], atol=1e-17)
    assert_allclose(angles, 0., atol=1e-13)
