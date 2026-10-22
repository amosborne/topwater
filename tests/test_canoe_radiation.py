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
