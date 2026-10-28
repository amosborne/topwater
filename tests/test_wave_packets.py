"""Packet construction, independent first-cutoff reference and Fourier identities."""
import numpy as np
import pytest

from topwater.wave_packets import TIME_STEP_S, _diagnostics, minimum_cutoff_packet, quintic_packet


@pytest.mark.parametrize('phase', [0., np.pi/4, np.pi/2])
def test_phase_is_direct_at_midpoint_and_quintic_envelope_is_symmetric(phase):
    packet = quintic_packet(.43, phase, 1., 1., 4., amplitude_m=2.)
    middle = np.flatnonzero(packet.time_s == packet.midpoint_s)[0]
    assert packet.midpoint_s == 4.
    assert packet.envelope[middle] == 1.
    assert packet.elevation_m[middle] == pytest.approx(2*np.cos(phase), abs=1e-14)
    np.testing.assert_allclose(packet.envelope[1:], packet.envelope[:0:-1], atol=1e-14)
    assert np.all(packet.elevation_m[packet.time_s <= 1.] == 0.)


def brute_force(frequency, phase, amplitude, tolerance):
    """Complex full FFT/IFFT, rather than the production incremental real sum."""
    packet = quintic_packet(frequency, phase, 1., 1., 4., amplitude_m=amplitude)
    n = len(packet.time_s)
    c = np.fft.fft(packet.elevation_m)/n
    signed = np.fft.fftfreq(n, packet.period_s/n)
    boundaries = np.array([0., 1., 2., 6., 7., 8.])
    envelope = np.r_[packet.envelope, [0., 0., 1., 1., 0., 0.]]
    gap = float(np.max(np.diff(np.unique(np.r_[packet.time_s, boundaries]))))
    rounding = 128*np.finfo(float).eps*n*(amplitude+2*np.sum(abs(c[:n//2+1])))
    start = np.searchsorted(np.arange(n//2+1)/packet.period_s, frequency)
    prior = []
    for k in range(start, n//2+1):
        retained = c.copy(); retained[abs(signed) > k/packet.period_s] = 0
        wave = np.fft.ifft(retained*n).real
        edge = (np.exp(2j*np.pi*boundaries[:, None]*signed)@retained).real
        peak = float(np.maximum(abs(np.r_[wave, edge])-amplitude*envelope, 0.).max())
        curvature = float(np.sum(abs(retained)*(2*np.pi*signed)**2))
        bound = peak+(curvature+amplitude*10/np.sqrt(3))*gap**2/8+rounding
        if bound <= tolerance:
            return k, wave, prior
        prior.append(bound)
    raise AssertionError('Reference needs a larger numerical grid')


@pytest.mark.parametrize('frequency,phase', [(.125, 0.), (.25, np.pi/2), (.43, .7)])
def test_first_certified_cutoff_matches_independent_full_fft_scan(frequency, phase):
    expected, wave, prior = brute_force(frequency, phase, 1., .02)
    result = minimum_cutoff_packet(frequency, phase, 1., 1., 4.,
        error_tolerance_m=.02)
    assert result.positive_components == expected
    assert result.cutoff_hz == expected/8.
    assert all(value > .02 for value in prior)
    np.testing.assert_allclose(result.elevation_m, wave, rtol=0, atol=1e-12)
    assert result.diagnostics['certified_envelope_exceedance_bound_m'] <= .02


def test_contiguous_lowpass_keeps_below_carrier_bins_dc_and_conjugates():
    result = minimum_cutoff_packet(.43, .7, 1., 1., 4.)
    k, n = result.positive_components, len(result.elevation_m)
    assert result.diagnostics['carrier_start_bin'] == 4
    assert abs(result.full_coefficients[1]) > 1e-6
    np.testing.assert_array_equal(result.coefficients[:k+1], result.full_coefficients[:k+1])
    np.testing.assert_array_equal(result.coefficients[1:k+1], result.coefficients[-1:-k-1:-1].conj())
    assert np.all(result.coefficients[k+1:n-k] == 0.)
    assert result.dc_m == pytest.approx(np.mean(result.target.elevation_m), abs=1e-15)
    np.testing.assert_allclose(np.fft.ifft(result.coefficients*n).imag, 0., atol=1e-14)
    assert result.full_coefficients[n//2].imag == 0.  # Unpaired Nyquist, not duplicated.


def test_odd_sample_count_and_exact_boundaries_are_fourier_evaluated():
    result = minimum_cutoff_packet(.43, .7, .83, 1.07, 3.31)
    t = result.diagnostics['boundary_times_s']
    n = len(result.elevation_m)
    assert n == 711
    assert t[-1] == result.target.period_s
    assert t[-1] > result.target.time_s[-1]
    np.testing.assert_array_equal(result.full_coefficients[1:n//2+1],
                                  result.full_coefficients[-1:n//2:-1].conj())
    np.testing.assert_allclose(np.fft.ifft(result.full_coefficients*n).real,
                               result.target.elevation_m, atol=1e-12)
    f = np.fft.fftfreq(len(result.elevation_m), result.target.period_s/len(result.elevation_m))
    expected = (np.exp(2j*np.pi*t[:, None]*f)@result.coefficients).real
    np.testing.assert_allclose(result.diagnostics['boundary_elevation_m'], expected, atol=1e-12)


def test_tolerance_is_absolute_metres_and_time_step_is_fixed():
    small = minimum_cutoff_packet(.43, .7, 1., 1., 4., amplitude_m=.25)
    large = minimum_cutoff_packet(.43, .7, 1., 1., 4., amplitude_m=2.)
    assert large.positive_components >= small.positive_components
    for result in (small, large):
        assert result.diagnostics['error_tolerance_m'] == .01
        assert result.diagnostics['certified_envelope_exceedance_bound_m'] <= .01
        assert result.diagnostics['time_step_s'] == TIME_STEP_S == .01
        assert result.diagnostics['sample_count'] == 800
        np.testing.assert_array_equal(result.target.time_s, np.arange(800)*.01)


def test_peak_exceedance_is_not_rms_or_dwell_fit():
    t = np.linspace(0., 3., 1001)
    envelope = ((t >= 1.) & (t <= 2.)).astype(float)
    wave = np.zeros(len(t)); wave[0] = .02
    report = _diagnostics(t, wave, envelope, envelope, 1., 1., 2., .5, 3.)
    assert np.sqrt(np.mean(wave**2)) < .01
    assert report['max_envelope_exceedance_m'] == .02
    assert report['dwell_max_error_m'] == 1.  # Diagnostic; not the acceptance criterion.


def test_only_one_final_ifft_and_one_forward_fft_not_one_per_candidate(monkeypatch):
    calls = dict(forward=0, inverse=0)
    forward, inverse = np.fft.rfft, np.fft.irfft
    def counted_forward(*args, **kwargs):
        calls['forward'] += 1
        return forward(*args, **kwargs)
    def counted_inverse(*args, **kwargs):
        calls['inverse'] += 1
        return inverse(*args, **kwargs)
    monkeypatch.setattr(np.fft, 'rfft', counted_forward)
    monkeypatch.setattr(np.fft, 'irfft', counted_inverse)
    result = minimum_cutoff_packet(.43, .7, 1., 1., 4.)
    assert calls == dict(forward=1, inverse=1)
    assert result.diagnostics['checked_cutoffs'] > 1


@pytest.mark.parametrize('argument,value', [('frequency_hz', 0.), ('phase_rad', np.nan),
    ('rest_s', -1.), ('ramp_s', 0.), ('dwell_s', -1.), ('amplitude_m', np.inf)])
def test_invalid_inputs_fail_clearly(argument, value):
    inputs = dict(frequency_hz=.43, phase_rad=.7, rest_s=1., ramp_s=1., dwell_s=4.)
    inputs[argument] = value
    with pytest.raises(ValueError):
        minimum_cutoff_packet(**inputs)


def test_underresolved_and_impossible_certification_fail_clearly():
    with pytest.raises(ValueError, match='Fixed 10ms grid'):
        quintic_packet(10., 0., 1., 1., 4.)
    with pytest.raises(ValueError, match='No cutoff through Nyquist'):
        minimum_cutoff_packet(.43, .7, 1., 1., 4.,
            error_tolerance_m=0.)
    with pytest.raises(ValueError, match='nonnegative'):
        minimum_cutoff_packet(.43, .7, 1., 1., 4., error_tolerance_m=-.01)


@pytest.mark.parametrize('argument,value', [('rest_s', 1.001), ('ramp_s', 1.005), ('dwell_s', 4.001)])
def test_timing_must_align_to_10ms_without_silent_rounding(argument, value):
    inputs = dict(frequency_hz=.43, phase_rad=.7, rest_s=1., ramp_s=1., dwell_s=4.)
    inputs[argument] = value
    with pytest.raises(ValueError, match='align to the fixed 0.01s'):
        minimum_cutoff_packet(**inputs)
