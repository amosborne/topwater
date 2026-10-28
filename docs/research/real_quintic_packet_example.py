"""Real quintic wave: direct Fourier, DC + topN positive + conjugates; no fit."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

T, CENTER, R, HALF_DWELL, N = 130., 65., 5., 50., 40
F0S, PHASES, BUDGETS = [.01, .1, 1.], [0., np.pi/4, np.pi/2], (10, 20, 40)


def smoothstep(x):
    return x**3*(10-15*x+6*x*x)


def envelope(u):
    return 1-smoothstep(np.clip((abs(u)-HALF_DWELL)/R, 0, 1))


def g_transform(offset, order):
    x, w = np.polynomial.legendre.leggauss(order)
    u = HALF_DWELL+R*(x+1)/2
    delta = np.asarray(offset)
    ramp = R/2*((w*(1-smoothstep((x+1)/2)))@np.cos(2*np.pi*u[:, None]*delta))
    return 2/T*(HALF_DWELL*np.sinc(2*HALF_DWELL*delta)+ramp)


def cosine_coefficients(frequencies, f0, phase, order):
    return (np.exp(1j*phase)*g_transform(frequencies-f0, order)
        + np.exp(-1j*phase)*g_transform(frequencies+f0, order))


def evaluate(u, frequency, coefficient, dc):
    chunk_size = max(1, min(3000, 200000//len(frequency)))
    return dc+np.concatenate([(np.exp(2j*np.pi*chunk[:, None]*frequency)@coefficient).real
        for chunk in np.array_split(u, max(1, int(np.ceil(len(u)/chunk_size))))])


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    out = Path(parser.parse_args().output); out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    t = np.linspace(0, T, 26001); u = t-CENTER; a = envelope(u)
    dwell, quiet = abs(u) <= 50, abs(u) >= 55
    ramp = (abs(u) >= 50) & (abs(u) <= 55)
    fft_size = 32768
    fft_u = np.arange(fft_size)*T/fft_size-CENTER
    arrays, reports = {'time_s': t, 'target_envelope': a}, []
    fig, axes = plt.subplots(3, 3, figsize=(17, 11), layout='constrained')
    specfig, specaxes = plt.subplots(3, 3, figsize=(17, 10), layout='constrained')
    phasefig, phaseaxes = plt.subplots(3, 3, figsize=(17, 10), layout='constrained')
    for i, (f0, phase) in enumerate(zip(F0S, PHASES)):
        jmax = max(80, int(np.ceil((f0+4/R)*T)))
        for expansion in range(6):
            indices = np.arange(1, jmax+1); frequencies = indices/T
            c = cosine_coefficients(frequencies, f0, phase, 256)
            check = cosine_coefficients(frequencies, f0, phase, 128)
            if not np.all(np.isfinite(c)) or not np.all(np.isfinite(check)):
                raise RuntimeError('Nonfinite direct quadrature coefficients')
            quadrature_difference = float(max(abs(c-check)))
            ranked = np.argsort(-abs(c), kind='stable')[:N]
            threshold = float(min(abs(c[ranked])))
            next_f = (jmax+1)/T
            constant = 80/(np.sqrt(3)*R**2*T*(2*np.pi)**3)
            tail_bound = constant*((next_f-f0)**-3+(next_f+f0)**-3)
            if tail_bound+2*quadrature_difference < threshold: break
            jmax *= 2
        else: raise RuntimeError('Candidate expansion did not separate top40 from omitted tail')
        dc = float(g_transform(np.array([f0]), 256)[0]*np.cos(phase))
        target = a*np.cos(2*np.pi*f0*u+phase)
        target_fft_grid = envelope(fft_u)*np.cos(2*np.pi*f0*fft_u+phase)
        target_fft = np.fft.fft(target_fft_grid)/fft_size
        candidate_reference = evaluate(u, frequencies, c, dc)
        full_reference_max_error = float(max(abs(candidate_reference-target)))
        for column, budget in enumerate(BUDGETS):
            chosen = ranked[:budget]; chosen = chosen[np.argsort(frequencies[chosen])]
            sf, sc, sj = frequencies[chosen], c[chosen], indices[chosen]
            if len(np.unique(sf)) != budget or np.any(sf <= 0) or np.any(sj >= fft_size/2):
                raise RuntimeError('Invalid positive component count or Nyquist collision')
            wave = evaluate(u, sf, sc, dc); error = wave-target
            if not np.isfinite(dc) or not np.all(np.isfinite(wave)) or not np.all(np.isfinite(error)):
                raise RuntimeError('Nonfinite real-wave reconstruction')
            spectrum = np.zeros(fft_size, dtype=complex)
            spectrum[0] = dc
            # IFFT origin is t0, whereas coefficients are defined directly at midpoint.
            spectrum[sj] = sc/2*np.exp(-2j*np.pi*sf*CENTER)
            spectrum[-sj] = np.conj(spectrum[sj])
            ifft_wave = np.fft.ifft(spectrum*fft_size)
            direct_fft_grid = evaluate(fft_u, sf, sc, dc)
            fft_reference_c = 2*target_fft[sj]*np.exp(2j*np.pi*sf*CENTER)
            reconstruction_curvature = float(sum(abs(sc)*(2*np.pi*sf)**2))
            omega = 2*np.pi*f0
            target_curvature_bound = 10/(np.sqrt(3)*R**2)+3.75*omega/R+omega**2
            allowance = (reconstruction_curvature+target_curvature_bound)*.005**2/8
            report = dict(carrier_frequency_hz=f0, component_budget=budget,
                phase_at_midpoint_rad=phase, phase_reference_time_s=CENTER,
                positive_component_count=budget, negative_conjugate_count=budget,
                dc_component_count=1, total_stored_fourier_coefficients=2*budget+1,
                dc_mean=dc, reconstructed_minus_target_mean=0., normalization_applied=False,
                frequencies_hz=sf.tolist(), positive_cosine_coefficient_real=sc.real.tolist(),
                positive_cosine_coefficient_imag=sc.imag.tolist(),
                midpoint_real_error=float(abs(evaluate(np.array([0.]), sf, sc, dc)[0]-np.cos(phase))),
                candidate_cutoff_hz=jmax/T, top40_threshold=threshold,
                omitted_positive_coefficient_bound=tail_bound, quadrature_128_256_max_difference=quadrature_difference,
                tail_separation_is_interval_proof=False,
                dwell_max_real_error=float(max(abs(error[dwell]))), quiet_max_real_elevation=float(max(abs(wave[quiet]))),
                ramp_max_real_error=float(max(abs(error[ramp]))), full_max_real_error=float(max(abs(error))),
                full_rms_real_error=float(np.sqrt(np.trapezoid(error**2, t)/T)),
                full_max_error_continuous_upper=float(max(abs(error)))+allowance,
                interpolation_allowance=allowance,
                candidate_full_series_max_error=full_reference_max_error,
                omitted_candidate_positive_energy=float(.5*(sum(abs(c)**2)-sum(abs(sc)**2))),
                retained_series_mean_square=float(dc**2+.5*sum(abs(sc)**2)),
                ifft_imaginary_residue=float(max(abs(ifft_wave.imag))),
                ifft_vs_direct_max_error=float(max(abs(ifft_wave.real-direct_fft_grid))),
                quadrature_vs_sampled_fft_selected_coefficient_error=float(max(abs(sc-fft_reference_c))),
                quadrature_vs_sampled_fft_dc_error=float(abs(dc-target_fft[0].real)),
                sampled_full_fft_ifft_roundtrip_error=float(max(abs(np.fft.ifft(target_fft*fft_size).real-target_fft_grid))))
            identity_keys = ('ifft_imaginary_residue', 'ifft_vs_direct_max_error',
                'sampled_full_fft_ifft_roundtrip_error')
            coefficient_keys = ('quadrature_vs_sampled_fft_selected_coefficient_error',
                'quadrature_vs_sampled_fft_dc_error')
            if any(not np.isfinite(report[key]) or report[key] > 1e-9 for key in identity_keys):
                raise RuntimeError('Numerical Fourier identity validation failed')
            if any(not np.isfinite(report[key]) or report[key] > 1e-8 for key in coefficient_keys):
                raise RuntimeError('Direct quadrature versus independent sampled FFT validation failed')
            reports.append(report)
            arrays.update({f'frequency_{i}_{budget}_hz': sf, f'positive_cosine_coefficients_{i}_{budget}': sc,
                f'ifft_full_spectrum_{i}_{budget}': spectrum, f'elevation_{i}_{budget}_m': wave,
                f'target_elevation_{i}_{budget}_m': target})
            axes[i, column].plot(t, wave, lw=.8, label=f'{budget}-positive-component wave')
            axes[i, column].plot(t, target, color='red', ls='--', lw=.65, label='Prescribed real quintic wave')
            axes[i, column].plot(t, a, color='green', ls=':', lw=.8, label='Prescribed envelope guides')
            axes[i, column].plot(t, -a, color='green', ls=':', lw=.8)
            for boundary in (10, 15, 115, 120): axes[i, column].axvline(boundary, color='grey', ls=':', lw=.5)
            axes[i, column].set(title=f'f0={f0:g} Hz; N={budget}; midpoint phase={np.degrees(phase):g} deg', xlabel='Time (s)', ylabel='Real elevation (m)')
            axes[i, column].legend(fontsize=6)
            specaxes[i, column].bar(sf, abs(sc), width=.003)
            specaxes[i, column].set(title=f'f0={f0:g} Hz; N={budget}; DC={dc:.4g} m', xlabel='Positive frequency (Hz)', ylabel='Cosine amplitude (m)')
            phaseaxes[i, column].bar(sf, np.angle(sc, deg=True), width=.003)
            phaseaxes[i, column].set(title=f'f0={f0:g} Hz; N={budget}; phase directly at65 s', xlabel='Positive frequency (Hz)', ylabel='Degrees', ylim=(-180,180))
    for figure, name, title in ((fig,'pulse-real-quintic-time.png','Real quintic wave; largest10/20/40 positive components; DC + conjugates free'),
            (specfig,'pulse-real-quintic-spectrum.png','Discrete positive cosine amplitudes; no normalization'),
            (phasefig,'pulse-real-quintic-phase.png','Positive cosine phases defined directly at midpoint65 s')):
        figure.suptitle(title); figure.savefig(out/name,dpi=160); plt.close(figure)
    np.savez(out/'packets.npz', **arrays)
    report = dict(policy='Ordinary REAL Fourier truncation, topN positive components; retain DC and negative conjugates, no scaling',
        timing_s=dict(rest=10,ramp=5,dwell=100,end=130,phase_reference=65),
        harmonic_spacing_hz=1/T, fft_reference_samples=fft_size, cases=reports,
        elapsed_s=time.perf_counter()-started, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limitations=['Largest coefficients minimize periodic mean-square truncation error, not peak error.',
            'Quadrature agreement and tail separation are floating-point engineering checks, not interval proof.',
            'Full sampled FFT roundtrip is only a sampled-curve identity; selected coefficient comparison is independent quadrature/FFT evidence.',
            'No physical accuracy gate or propagation; construction repeats every130s; mean preserved but finiteN midpoint value/phase need not match exactly.'])
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__ == '__main__': main()
