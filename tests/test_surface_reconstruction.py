"""Analytical field-convention fixtures; no Capytaine or MATLAB jobs.

Execution still belongs to the coordinator's serial numerical allocation.
"""
import numpy as np
import pytest
from types import SimpleNamespace

from topwater._capytaine_surface import (
    SurfaceKernel, convolve_radiation, recover_transfer, synthesize_spectrum,
    require_surface_readiness, cosine_kernel, exterior_field_points,
    exterior_trace_points, incident_spectrum,
)


def exponential_kernel(*, accepted=True):
    t = np.arange(20001)*.001
    values = np.zeros((len(t), 1, 6))
    values[:, 0, 0] = 3*np.exp(-2*t)
    direct = np.zeros((1, 6)); direct[0, 0] = .4
    return SurfaceKernel(t, values, direct, accepted)


def test_exponential_plus_direct_has_correct_complex_phase():
    kernel = exponential_kernel()
    omega = np.array([.07, .31, 1.13, 3.71, 7.19])
    actual = recover_transfer(kernel, omega)[:, 0, 0]
    expected = .4 + 3/(2-1j*omega)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-6)
    assert np.all(actual.imag > 0), 'exp(-iwt) uses the positive-imaginary causal transform'


def test_cosine_inversion_has_bounded_finite_band_error_and_refines():
    # h=A exp(-a*t), H=D+A/(a-i*w).  The infinite-frequency tail is bounded
    # uniformly in lag by (2*A/pi)*(pi/2-atan(W/a)).  The trapezoid bound follows
    # from |(R cos(w*t))''| <= 10*A/a^3+4*A*t/a^2+A*t^2/a,
    # where R=A*a/(a^2+w^2).  Thus the test budgets actual truncation and
    # quadrature error instead of comparing a finite band to an exact identity.
    A, a, D = 3., 2., .4
    lag = np.linspace(0., .25, 21)
    errors = []
    for maximum, spacing in [(64., .02), (256., .01)]:
        omega = np.arange(round(maximum/spacing)+1)*spacing
        transfer = np.zeros((len(omega), 1, 6), dtype=complex)
        transfer[:, 0, 0] = D + A/(a-1j*omega)
        direct = np.zeros((1, 6)); direct[0, 0] = D
        candidate = cosine_kernel(omega, transfer, lag, direct_velocity=direct)
        error = float(np.max(abs(candidate.values[:, 0, 0]-A*np.exp(-a*lag))))
        tail = 2*A/np.pi*(np.pi/2-np.arctan(maximum/a))
        derivative_bound = 10*A/a**3 + 4*A*lag[-1]/a**2 + A*lag[-1]**2/a
        quadrature = 2/np.pi*maximum*spacing**2/12*derivative_bound
        assert error <= tail+quadrature+1e-12
        np.testing.assert_array_equal(candidate.direct_velocity, direct)
        assert not candidate.accepted
        errors.append(error)
    assert errors[1] < errors[0]/3


def test_mean_footprint_guard_rejects_interior_and_on_waterline():
    vertices = np.array([[-2., -1., 0.], [2., -1., 0.], [2., 1., 0.],
                         [-2., 1., 0.], [0., 0., -1.]])
    hydro = SimpleNamespace(body=SimpleNamespace(mesh=SimpleNamespace(vertices=vertices)))
    exterior = np.array([[2.1, 0.], [-3., .5], [0., 1.1]])
    np.testing.assert_array_equal(exterior_field_points(hydro, exterior), exterior)
    for point in ([[0., 0.]], [[2., 0.]], [[-2., -1.]]):
        with pytest.raises(ValueError, match='strictly outside'):
            exterior_field_points(hydro, point)


def test_exterior_trace_maps_only_disturbance_queries():
    vertices = np.array([[-2., -1., 0.], [2., -1., 0.], [2., 1., 0.],
                         [-2., 1., 0.], [0., 0., -1.]])
    hydro = SimpleNamespace(body=SimpleNamespace(mesh=SimpleNamespace(vertices=vertices)))
    original = np.array([[1.9, 0.], [3., 0.]])
    mapped, continued = exterior_trace_points(hydro, original, .025)
    np.testing.assert_allclose(mapped, [[2.025, 0.], [3., 0.]], atol=1e-14)
    np.testing.assert_array_equal(continued, [True, False])
    train = SimpleNamespace(omega=np.array([0., 1.]), k=np.array([0., .2]),
                            theta=np.array([0.]), C=np.array([[[0., 1.], [0., 0.]]]))
    np.testing.assert_allclose(incident_spectrum(train, original)[1], np.exp(1j*.2*original[:, 0]))
    assert abs(incident_spectrum(train, mapped)[1, 0]-incident_spectrum(train, original)[1, 0]) > .01


def test_identified_acceleration_feedthrough_is_explicit_and_uses_actual_acceleration():
    kernel = exponential_kernel()
    kernel.direct_acceleration = np.zeros((1, 6)); kernel.direct_acceleration[0, 0] = .07
    omega = np.array([.17, 1.13, 3.71])
    np.testing.assert_allclose(recover_transfer(kernel, omega)[:, 0, 0],
                               .4+3/(2-1j*omega)-1j*omega*.07, atol=2e-6)
    t = kernel.lag_s[:101]
    velocity = np.zeros((len(t), 6))
    acceleration = np.zeros_like(velocity); acceleration[:, 0] = 2.
    with pytest.raises(ValueError, match='Actual modal accelerations'):
        convolve_radiation(kernel, t, velocity)
    np.testing.assert_allclose(convolve_radiation(kernel, t, velocity, acceleration)[:, 0], .14, atol=1e-14)


def test_causal_convolution_has_no_future_motion_or_startup_history():
    kernel = exponential_kernel()
    t = kernel.lag_s[:4001]
    w = 1.7
    v = np.zeros((len(t), 6)); v[:, 0] = np.sin(w*t)
    actual = convolve_radiation(kernel, t, v)[:, 0]
    expected = .4*np.sin(w*t) + 3*(2*np.sin(w*t)-w*np.cos(w*t)+w*np.exp(-2*t))/(4+w*w)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-6)
    changed = v.copy(); changed[2001:, 0] += 100
    np.testing.assert_allclose(convolve_radiation(kernel, t, changed)[:2001, 0], actual[:2001], rtol=0, atol=1e-11)
    assert actual[0] == pytest.approx(0, abs=1e-14)


def test_unaccepted_or_short_kernel_cannot_publish_radiation():
    kernel = exponential_kernel(accepted=False)
    t = kernel.lag_s[:11]
    with pytest.raises(ValueError, match='review'):
        convolve_radiation(kernel, t, np.zeros((len(t), 6)))
    kernel.accepted = True
    kernel.lag_s = kernel.lag_s[:5]
    kernel.values = kernel.values[:5]
    with pytest.raises(ValueError, match='entire requested history'):
        convolve_radiation(kernel, t, np.zeros((len(t), 6)))


def test_weighted_spectrum_has_no_hidden_factor_or_conjugation():
    omega = np.array([0., .3, .9])
    coefficients = np.array([[0.], [2+3j], [-.7+.2j]])
    t = np.array([0., .17, 1.2])
    expected = sum((coefficients[i, 0]*np.exp(-1j*w*t)).real for i, w in enumerate(omega))
    np.testing.assert_allclose(synthesize_spectrum(omega, coefficients, t)[:, 0], expected, atol=1e-14)


def test_total_surface_capability_remains_closed_before_physical_pilot():
    with pytest.raises(ValueError, match='pilot'):
        require_surface_readiness()
