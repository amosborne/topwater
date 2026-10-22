"""Source-review controls; no time integration or factorization here."""
import numpy as np
import json
import pytest

from .experiments.free_surface.pml_strip import stretching
from .experiments.free_surface.sponge_strip import make_strip
from .experiments.free_surface.time_pml_strip import WIDTH, drive, left_profile, stationary_field


def test_pulse_is_zero_after_drive_but_excites_dc_memory():
    omega = 2.;period = 2*np.pi/omega
    times = np.linspace(0., 8*period, 16001)
    pulse = drive(times, omega)
    np.testing.assert_array_equal(pulse, drive(times, omega, plateau_periods=4.))
    np.testing.assert_allclose(drive(32*period, omega, plateau_periods=32.), 1.2)
    assert drive(34*period, omega, plateau_periods=32.) == 0.
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


def test_failed_run_preserves_only_completed_history_and_original_error(tmp_path, monkeypatch):
    from .experiments.free_surface import time_pml_strip as module
    strip = make_strip(2., np.pi/4, WIDTH, spacing=3., layers=4, quiet_wavelengths=.1)
    output = tmp_path/'new'/'failed.json'
    failure = FloatingPointError('synthetic failed midpoint')

    class FakeMidpoint:
        def __init__(self, supplied, dt, *, plateau_periods=4.):
            assert plateau_periods == 32.
            assert output.parent.is_dir()  # Created before any numerical work.
            self.dt = dt;self.time = 0.
            self.free = np.unique(supplied.top);self.left = np.unique(supplied.left)
            self.psi = np.zeros(len(self.free));self.v = self.psi.copy()
            self.chi = np.zeros((1, 1));self.z = self.chi.copy();self.q = self.chi.copy()

        def step(self):
            if self.time:
                raise failure
            self.time += self.dt
            self.chi = np.full((1, 1), .3)
            return np.ones(len(strip.points))

        def diagnostics(self):
            return np.array([1., .3, 0., 0.])

    monkeypatch.setattr(module, 'make_strip', lambda *args, **kwargs: strip)
    monkeypatch.setattr(module, 'MidpointStrip', FakeMidpoint)
    with pytest.raises(FloatingPointError) as caught:
        module.run(np.pi/100, output, plateau_periods=32.)
    assert caught.value is failure
    report = json.loads(output.read_text())
    assert report['complete'] is False and report['valid_midpoints'] == 1
    assert report['duration_periods'] == 40. and report['fit_period_range'] == [31., 33.]
    assert report['pulse_plateau_periods'] == 32.
    assert report['error']['type'] == 'FloatingPointError'
    assert 'insufficient' in report['harmonic_fit']
    assert report['runtime']['git_revision'] and len(report['source_sha256']) == 4
    assert len(report['parameters_sha256']) == 64
    with np.load(output.with_suffix('.npz')) as arrays:
        assert arrays['potential'].shape == (1, len(strip.points))
        np.testing.assert_array_equal(arrays['endpoint_times'], [0., np.pi/100])
        np.testing.assert_array_equal(arrays['final_chi'], [[.3]])
