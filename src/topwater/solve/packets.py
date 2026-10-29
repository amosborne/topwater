"""Real quintic packets and ascending physical-cutoff search on declared grids."""
from dataclasses import dataclass

import numpy as np


from topwater.solve.settings import DEFAULTS

TIME_STEP_S = DEFAULTS.time_step_s


@dataclass(frozen=True)
class PacketTarget:
    time_s: np.ndarray
    envelope: np.ndarray
    elevation_m: np.ndarray
    period_s: float
    midpoint_s: float


@dataclass(frozen=True)
class CutoffPacket:
    target: PacketTarget
    cutoff_hz: float
    positive_components: int
    coefficients: np.ndarray
    full_coefficients: np.ndarray
    elevation_m: np.ndarray
    diagnostics: dict

    @property
    def dc_m(self):
        return float(self.coefficients[0].real)


def _values(t, frequency_hz, phase_rad, rest_s, ramp_s, dwell_s, amplitude_m):
    start, end = rest_s+ramp_s, rest_s+ramp_s+dwell_s
    def smoothstep(u):
        u = np.clip(u, 0., 1.)
        return u**3*(10.-15.*u+6.*u*u)
    envelope = smoothstep((t-rest_s)/ramp_s)*smoothstep((end+ramp_s-t)/ramp_s)
    return envelope, amplitude_m*envelope*np.cos(2*np.pi*frequency_hz*(t-(start+end)/2)+phase_rad)


def quintic_packet(frequency_hz, phase_rad, rest_s, ramp_s, dwell_s, *,
                   amplitude_m=1., time_step_s=TIME_STEP_S):
    """Sample rest/ramp/dwell/ramp/rest; phase is DIRECTLY at dwell midpoint."""
    values = np.asarray([frequency_hz, phase_rad, rest_s, ramp_s, dwell_s, amplitude_m], dtype=float)
    if not np.isfinite(values).all() or frequency_hz <= 0 or rest_s < 0 or ramp_s <= 0 or dwell_s < 0 or amplitude_m <= 0:
        raise ValueError('Finite inputs require frequency/amplitude/ramp > 0 and rest/dwell >= 0')
    if not np.isfinite(time_step_s) or time_step_s <= 0:
        raise ValueError('time_step_s must be positive and finite')
    timing = np.asarray([rest_s, ramp_s, dwell_s], dtype=float)
    steps = np.rint(timing/time_step_s)
    if not np.isfinite(steps).all() or np.any(steps > np.iinfo(np.intp).max//5):
        raise ValueError('Timing exceeds the representable timestep grid')
    alignment_roundoff = 8*np.finfo(float).eps*np.maximum(abs(timing), time_step_s)
    if np.any(abs(timing-steps*time_step_s) > alignment_roundoff):
        raise ValueError(f'rest_s, ramp_s and dwell_s must each align to the {time_step_s:g}s timestep')
    rest_steps, ramp_steps, dwell_steps = (int(value) for value in steps)
    n = 2*rest_steps+2*ramp_steps+dwell_steps
    period = n*time_step_s
    if n < 64 or ramp_steps < 16 or time_step_s*frequency_hz > 1/16:
        raise ValueError('Timestep grid requires >=64 total samples and >=16 samples per carrier period and ramp')
    t = np.arange(n)*time_step_s
    envelope, elevation = _values(t, frequency_hz, phase_rad, rest_s, ramp_s, dwell_s, amplitude_m)
    if not np.isfinite(elevation).all():
        raise ValueError('Nonfinite packet target')
    return PacketTarget(t, envelope, elevation, period, period/2)


def _diagnostics(t, wave, envelope, truth, amplitude, start, end, rest, period):
    dwell, quiet = (t >= start) & (t <= end), (t <= rest) | (t >= period-rest)
    ramp = ~(dwell | quiet)
    error = abs(wave-truth)
    exceedance = np.maximum(abs(wave)-amplitude*envelope, 0.)
    worst = int(np.argmax(exceedance))
    dwell_index = np.flatnonzero(dwell)[np.argmax(error[dwell])]
    return dict(max_envelope_exceedance_m=float(exceedance[worst]),
        max_envelope_exceedance_time_s=float(t[worst]),
        dwell_max_error_m=float(error[dwell_index]), dwell_worst_time_s=float(t[dwell_index]),
        rest_max_elevation_m=float(np.max(abs(wave[quiet]))),
        ramp_max_template_error_m=float(np.max(error[ramp])) if ramp.any() else 0.,
        whole_window_max_template_error_m=float(error.max()))


def minimum_cutoff_packet(frequency_hz, phase_rad, rest_s, ramp_s, dwell_s, *,
                          amplitude_m=1., error_tolerance_m=DEFAULTS.packet_tolerance_m, time_step_s=TIME_STEP_S, _harmonics=None):
    """First eligible cutoff whose conservative envelope-exceedance bound passes.

    Eligible cutoffs start at ceil(frequency_hz*period). All lower bins are
    INCLUDED in synthesis. Only |wave| <= amplitude*envelope + tolerance gates
    acceptance; dwell/template mismatch is diagnostic. No continuum minimum is
    claimed: lower cutoffs rejected only by the interpolation guard may pass.
    """
    if not np.isfinite(error_tolerance_m) or error_tolerance_m < 0:
        raise ValueError('error_tolerance_m must be finite and nonnegative metres')
    target = quintic_packet(frequency_hz, phase_rad, rest_s, ramp_s, dwell_s,
                            amplitude_m=amplitude_m, time_step_s=time_step_s)
    n, period = len(target.time_s), target.period_s
    start, end = rest_s+ramp_s, rest_s+ramp_s+dwell_s
    boundaries = np.unique([0., rest_s, start, end, period-rest_s, period])
    boundary_envelope, boundary_truth = _values(boundaries, frequency_hz, phase_rad,
                                               rest_s, ramp_s, dwell_s, amplitude_m)
    t = np.r_[target.time_s, boundaries]
    harmonic = None if _harmonics is None else _harmonics.prepare(t, period)
    envelope = np.r_[target.envelope, boundary_envelope]
    truth = np.r_[target.elevation_m, boundary_truth]
    full = np.fft.rfft(target.elevation_m)/n
    if not np.isfinite(full).all():
        raise ValueError('Nonfinite Fourier coefficients')
    wave = np.full(len(t), full[0].real)
    carrier_start = int(np.searchsorted(np.arange(len(full))/period, frequency_hz, side='left'))
    gap = float(np.max(np.diff(np.unique(t))))
    envelope_curvature = amplitude_m*10/(np.sqrt(3)*ramp_s**2)
    roundoff = 128*np.finfo(float).eps*n*(amplitude_m+2*np.sum(abs(full)))
    curvature, history = 0., []
    for k in range(len(full)):
        if k:
            factor = 1 if n % 2 == 0 and k == n//2 else 2
            basis = np.exp(2j*np.pi*k*t/period) if harmonic is None else harmonic(k)
            wave += factor*(full[k]*basis).real
            curvature += factor*abs(full[k])*(2*np.pi*k/period)**2
        if k < carrier_start:
            continue
        if not np.isfinite(wave).all() or not np.isfinite(curvature):
            raise ValueError('Nonfinite incremental Fourier reconstruction')
        peak = float(np.maximum(abs(wave)-amplitude_m*envelope, 0.).max())
        allowance = float((curvature+envelope_curvature)*gap**2/8)
        bound = peak+allowance+roundoff
        if bound > error_tolerance_m:
            history.append((k, peak, bound))
            continue
        kept = np.zeros_like(full); kept[:k+1] = full[:k+1]
        final = np.fft.irfft(kept*n, n=n)  # Exactly ONE IFFT after search passes.
        identity_error = float(np.max(abs(final-wave[:n])))
        if not np.isfinite(final).all() or identity_error > roundoff:
            raise ValueError('Final IFFT disagrees with incremental real Fourier sum')
        metrics = _diagnostics(t, np.r_[final, wave[n:]], envelope, truth,
                               amplitude_m, start, end, rest_s, period)
        certified = metrics['max_envelope_exceedance_m']+allowance+roundoff
        if not np.isfinite(certified) or certified > error_tolerance_m:
            raise ValueError('Final IFFT validation changed the passing envelope bound')
        metrics.update(error_tolerance_m=float(error_tolerance_m),
            carrier_start_bin=carrier_start, positive_bin_count=k,
            checked_cutoffs=k-carrier_start+1, prior_cutoff_exceedance_and_bound=np.asarray(history).reshape(-1, 3),
            lower_cutoffs_unresolved_by_guard=sum(row[1] <= error_tolerance_m for row in history),
            between_samples_allowance_m=allowance, roundoff_allowance_m=float(roundoff),
            certified_envelope_exceedance_bound_m=float(certified),
            final_ifft_max_difference_m=identity_error,
            boundary_times_s=boundaries.copy(), boundary_elevation_m=wave[n:].copy(),
            validation_grid_max_gap_s=gap, time_step_s=time_step_s, sample_count=n)
        if not all(np.isfinite(metrics[key]) for key in ('dwell_max_error_m', 'rest_max_elevation_m',
                'ramp_max_template_error_m', 'whole_window_max_template_error_m')):
            raise ValueError('Nonfinite error diagnostics')
        mirror = lambda c: np.r_[c, (c[-2:0:-1] if n % 2 == 0 else c[-1:0:-1]).conj()]
        return CutoffPacket(target, k/period, k, mirror(kept), mirror(full), final, metrics)
    raise ValueError(f'No cutoff through Nyquist certifies the envelope tolerance on the {time_step_s:g}s grid; reduce time_step_s or increase error_tolerance_m')
