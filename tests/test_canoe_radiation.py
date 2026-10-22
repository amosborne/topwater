"""Independent forcing calculus and radiation sign/extraction checks."""
import numpy as np
from tests.experiments.free_surface.canoe_radiation import heave, harmonic_fit, sponge


def test_quintic_heave_derivatives_and_smooth_join():
    amplitude, omega = .003, 2.7
    duration = 4*np.pi/omega
    for time in (.31*duration, .73*duration, duration, 1.4*duration):
        step = 1e-6
        value = heave(time, amplitude, omega)
        plus, minus = heave(time+step, amplitude, omega), heave(time-step, amplitude, omega)
        np.testing.assert_allclose((plus[:2]-minus[:2])/(2*step), value[1:], rtol=2e-7, atol=2e-11)
    np.testing.assert_array_equal(heave(0., amplitude, omega), 0.)
    expected = amplitude*np.array([np.cos(omega*duration), -omega*np.sin(omega*duration),
                                    -omega**2*np.cos(omega*duration)])
    np.testing.assert_allclose(heave(duration, amplitude, omega), expected, atol=1e-17)


def test_harmonic_force_sign_restoring_and_adaptive_time_sampling():
    omega, amplitude, restoring, added, damping = 3., .002, 950., 110., 12.
    # Independent rigid-fluid law: F = -A*zddot - B*zdot - C*z + DC.
    times = np.linspace(0, 1, 901)**1.4*4*np.pi/omega
    force = 70+amplitude*(added*omega**2-restoring)*np.cos(omega*times)+amplitude*omega*damping*np.sin(omega*times)
    result = harmonic_fit(times, force, omega, amplitude, restoring)
    np.testing.assert_allclose([result['added_mass'], result['radiation_damping']], [added, damping], rtol=2e-12)
    assert result['harmonic_residual_rms'] < 1e-12
    np.testing.assert_allclose(result['mean_fluid_work'], -.5*damping*(amplitude*omega)**2, rtol=1e-3)


def test_sponge_zero_in_pocket_and_smooth_bounded_outer_profile():
    points = np.array([[0., 0., 0.], [2*np.pi/3, np.pi/3, 0.], [6., 3., 0.], [-6., -3., 0.]])
    sigma, widths = sponge(points, 4.)
    np.testing.assert_array_equal(sigma[:2], 0.)
    np.testing.assert_array_equal(sigma[2:], 4.)
    assert np.all(widths > 0)


def test_all_prescribed_modes_have_correct_world_rigid_velocity_and_acceleration():
    from scipy.spatial.transform import Rotation
    from tests.experiments.free_surface.canoe_radiation import prescribed_state
    base = np.r_[[.4, -.2, .1], Rotation.from_rotvec([.4, -.3, .2]).as_quat(scalar_first=True),
                 np.ones(6), np.zeros(4)]
    fluid = np.array([.1, .2, .3, .4])
    body_point = np.array([.8, -.25, .17])
    time, step = .83, 1e-4
    for mode in range(6):
        def sample(t):
            state, acceleration = prescribed_state(base, fluid, heave(t, .03, 2.7), mode)
            offset = Rotation.from_quat(state[3:7], scalar_first=True).apply(body_point)
            return state, acceleration, state[:3]+offset, offset
        state, acceleration, point, offset = sample(time)
        plus, minus = sample(time+step)[2], sample(time-step)[2]
        expected_velocity = state[7:10]+np.cross(state[10:13], offset)
        expected_acceleration = (acceleration[:3]+np.cross(acceleration[3:], offset)+
                                 np.cross(state[10:13], np.cross(state[10:13], offset)))
        np.testing.assert_allclose((plus-minus)/(2*step), expected_velocity, atol=1e-9, rtol=1e-7)
        np.testing.assert_allclose((plus-2*point+minus)/step**2, expected_acceleration, atol=3e-7, rtol=2e-5)
        np.testing.assert_array_equal(state[13:], fluid)
        assert np.count_nonzero(state[7:13]) == 1


def test_full_radiation_column_extraction_includes_signed_couplings():
    from tests.experiments.free_surface.canoe_radiation import column_errors
    times = np.linspace(0., 4*np.pi/3, 501)
    added = np.array([2., -.7, 130., .4, -3., .09])
    damping = np.array([.1, -.03, 18., -.04, .7, -.01])
    restoring = np.array([0., 0., 1000., 2., -4., 0.])
    amplitude, omega = .002, 3.
    actual = []
    for row in range(6):
        forces = (row*10+amplitude*(omega**2*added[row]-restoring[row])*np.cos(omega*times)+
                  amplitude*omega*damping[row]*np.sin(omega*times))
        result = harmonic_fit(times, forces, omega, amplitude, restoring[row])
        actual.append([result['added_mass'], result['radiation_damping']])
    np.testing.assert_allclose(actual, np.column_stack((added, damping)), atol=1e-11)
    np.testing.assert_allclose(column_errors(added, added), 0.)
    assert column_errors(np.ones(6), np.zeros(6)) == [None]*6
    assert column_errors(np.ones(6), np.r_[np.ones(3), np.zeros(3)])[3:] == [None]*3
