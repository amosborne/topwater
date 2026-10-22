"""Source-review controls; no time integration or factorization here."""
import numpy as np

from .experiments.free_surface.pml_strip import stretching
from .experiments.free_surface.sponge_strip import make_strip
from .experiments.free_surface.time_pml_strip import WIDTH, drive, left_profile, stationary_field


def test_pulse_is_zero_after_drive_but_excites_dc_memory():
    omega = 2.;period = 2*np.pi/omega
    times = np.linspace(0., 8*period, 16001)
    pulse = drive(times, omega)
    assert pulse[0] == 0.
    np.testing.assert_array_equal(pulse[times >= 6*period], 0.)
    # Integral of the chosen finite pulse is nonzero; an accidentally balanced
    # carrier alone would not exercise the zero-frequency memory poles.
    assert np.trapezoid(pulse, times) > .5*period
    np.testing.assert_allclose(drive(np.array([2., 3., 4.])*period, omega), 1.2)


def test_same_dirichlet_stationary_field_satisfies_transformed_pde_and_wall():
    strip = make_strip(2., np.pi/4, WIDTH, spacing=2., layers=4, quiet_wavelengths=.1)
    left = np.unique(strip.left)
    np.testing.assert_allclose(stationary_field(strip, strip.points[left]),
                               left_profile(strip)[left], atol=3e-15)
    # Independent finite differences of the analytic field test all three
    # transformed volume coefficients, not only the imposed boundary trace.
    point = np.array([[WIDTH*.43, -1.1]])
    delta = 2e-4
    ex = np.array([[delta, 0.]]);ez = np.array([[0., delta]])
    f = stationary_field(strip, point)
    fxp = stationary_field(strip, point+ex);fxm = stationary_field(strip, point-ex)
    fzp = stationary_field(strip, point+ez);fzm = stationary_field(strip, point-ez)
    _, s, _ = stretching(point[:, 0], WIDTH, 8., strip.omega)
    _, sp, _ = stretching(point[:, 0]+delta/2, WIDTH, 8., strip.omega)
    _, sm, _ = stretching(point[:, 0]-delta/2, WIDTH, 8., strip.omega)
    residual = ((fxp-f)/sp-(f-fxm)/sm)/delta**2
    residual += s*(fzp-2*f+fzm)/delta**2-strip.ky**2*s*f
    assert np.max(abs(residual)) < 2e-7
    wall = np.array([[WIDTH, -1.1]])
    wall_gradient = (stationary_field(strip, wall+ex)-stationary_field(strip, wall-ex))/(2*delta)
    assert np.max(abs(wall_gradient)) < 2e-8
