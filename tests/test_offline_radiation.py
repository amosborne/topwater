"""Small independent Fourier controls; coordinator must allocate test execution."""
import numpy as np
import pytest

from topwater._capytaine_surface import (
    OfflineRadiationSettings, offline_radiation_response,
    cosine_kernel, convolve_radiation,
)


def compact_case():
    dt, period = 1/32, 8.
    stamp = np.arange(33)*dt
    v = np.zeros((len(stamp), 6)); a = np.zeros_like(v)
    v[:, 0] = stamp**2*(1-stamp)**2
    a[:, 0] = 2*stamp-6*stamp**2+4*stamp**3
    v[:, 4] = -.7*v[:, 0]; a[:, 4] = -.7*a[:, 0]
    w = 2*np.pi*np.arange(129)/period
    H = np.zeros((len(w), 2, 6), dtype=complex)
    H[1:, 0, 0] = .7+.4j
    H[1:, 0, 4] = -.2-.6j
    H[:, 1] = 1.3*H[:, 0]
    settings = OfflineRadiationSettings(period_s=period, fft_dt_s=dt,
        future_s=1., maximum_omega_rad_s=w[-2])
    return stamp, v, a, w, H, settings


def test_full_complex_response_matches_independent_dft_and_retains_lost_phase():
    t, v, a, w, H, settings = compact_case()
    actual, evidence = offline_radiation_response(w, H, t, v, a, settings=settings)
    # Independent direct Fourier SUM, not the implementation FFT. Motion and
    # final v/a are zero outside the saved prefix, so the continuation is zero.
    n, origin = 256, 64
    bins = np.arange(1, 128)
    angle = 2*np.pi*bins[:, None]*(origin+np.arange(len(t)))[None]/n
    Vf = np.exp(-1j*angle) @ v
    Yf = np.einsum('fpm,fm->fp', np.conj(H[bins]), Vf)
    expected = 2/n*(np.exp(1j*angle).T @ Yf).real
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-13)
    # Minimal before/after witness: identical real response, different imaginary
    # response. The retired cosine-only method cannot distinguish these inputs.
    real_only = H.real.astype(complex)
    lost_phase, _ = offline_radiation_response(w, real_only, t, v, a, settings=settings)
    assert np.max(np.abs(actual-lost_phase)) > 1e-3
    kernels = [cosine_kernel(w, data, t, direct_velocity=np.zeros((2, 6)))
               for data in (H, real_only)]
    for kernel in kernels:
        kernel.accepted = True
    np.testing.assert_array_equal(convolve_radiation(kernels[0], t, v),
                                  convolve_radiation(kernels[1], t, v))
    assert evidence['initial_sample_modified'] is False
    assert evidence['universal_causality_claim'] is False
    assert evidence['numerical_accuracy_qualified'] is False


def test_future_hermite_uses_final_acceleration_without_adding_direct_force():
    t, v, a, w, H, settings = compact_case()
    # An exact advance is deliberately noncausal. It makes the supplied future
    # observable and checks that the routine does not replace it by zero or fold
    # its phase into a causal cosine kernel.
    H[:] = 0
    H[:, 0, 0] = np.exp(-1j*w*.25)
    settings = OfflineRadiationSettings(period_s=8., fft_dt_s=1/32,
        future_s=1., maximum_omega_rad_s=w[-1])
    v[:, 0] = t; a[:, 0] = 1.
    response, evidence = offline_radiation_response(w, H, t, v, a, settings=settings)
    u = .25
    expected_future = (2*u**3-3*u**2+1)*v[-1, 0]+(u**3-2*u**2+u)*a[-1, 0]
    np.testing.assert_allclose(response[-1, 0], expected_future, rtol=0, atol=1e-12)
    np.testing.assert_allclose(response[0, 0], .25, rtol=0, atol=1e-12)
    assert evidence['raw_initial_radiation_m'][0] == pytest.approx(.25)
    assert evidence['initial_sample_modified'] is False
    assert evidence['future_policy'].startswith('C1 Hermite')


def test_acceleration_transfer_is_applied_once_as_part_of_full_H():
    t, v, a, w, H, settings = compact_case()
    H[:] = 0
    H[:, 0, 0] = -1j*w*.37
    actual, _ = offline_radiation_response(w, H, t, v, a, settings=settings)
    # Independent derivative of the same finite Fourier representation. This
    # intentionally differs from claiming it equals every saved acceleration.
    bins = np.arange(1, 128); n, origin = 256, 64
    angle = 2*np.pi*bins[:, None]*(origin+np.arange(len(t)))[None]/n
    Vf = np.exp(-1j*angle) @ v[:, 0]
    expected = 2/n*(np.exp(1j*angle).T @ (1j*w[bins]*.37*Vf)).real
    np.testing.assert_allclose(actual[:, 0], expected, rtol=0, atol=1e-12)


def test_offline_policy_rejects_unrepresented_knots_or_short_frequency_bank():
    t, v, a, w, H, settings = compact_case()
    with pytest.raises(ValueError, match='frequency extrapolation'):
        offline_radiation_response(w[:-2], H[:-2], t, v, a, settings=settings)
    with pytest.raises(ValueError, match='exactly represented'):
        offline_radiation_response(w, H, t, v, a,
            settings=OfflineRadiationSettings(period_s=8., fft_dt_s=1/16,
                future_s=1., maximum_omega_rad_s=40.))
