"""Fixed quintic -> direct Fourier coefficients -> top40 positive lines. No fit."""
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
F0S, PHASES = [.01, .1, 1.], [0., np.pi/4, np.pi/2]


def smoothstep(x):
    return x**3*(10-15*x+6*x*x)


def envelope(u):
    return 1-smoothstep(np.clip((abs(u)-HALF_DWELL)/R, 0, 1))


def coefficients(frequencies, f0, order):
    x, w = np.polynomial.legendre.leggauss(order)
    u = HALF_DWELL+R*(x+1)/2
    delta = np.asarray(frequencies)-f0
    ramp = R/2*((w*(1-smoothstep((x+1)/2)))@np.cos(2*np.pi*u[:, None]*delta))
    flat = HALF_DWELL*np.sinc(2*HALF_DWELL*delta)
    return 2/T*(flat+ramp)


def evaluate(u, frequency, coefficient):
    return np.concatenate([np.exp(2j*np.pi*chunk[:, None]*frequency)@coefficient
        for chunk in np.array_split(u, max(1, int(np.ceil(len(u)/3000))))])


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    out = Path(parser.parse_args().output); out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    t = np.linspace(0, T, 26001); u = t-CENTER; target_envelope = envelope(u)
    dwell, quiet = abs(u) <= 50, abs(u) >= 55
    ramp = (abs(u) >= 50) & (abs(u) <= 55)
    x, w = np.polynomial.legendre.leggauss(256)
    total_energy = 2/T*(50+R/2*np.sum(w*(1-smoothstep((x+1)/2))**2))
    arrays, reports = {'time_s': t, 'target_envelope': target_envelope}, []
    fig, axes = plt.subplots(3, 3, figsize=(17, 11), layout='constrained')
    specfig, specaxes = plt.subplots(3, 3, figsize=(17, 10), layout='constrained')
    phasefig, phaseaxes = plt.subplots(3, 3, figsize=(17, 10), layout='constrained')
    for i, (f0, phase) in enumerate(zip(F0S, PHASES)):
        jmax = max(80, int(np.ceil((f0+4/R)*T)))
        for expansion in range(6):
            frequencies = np.arange(1, jmax+1)/T
            raw = coefficients(frequencies, f0, 256)
            check = coefficients(frequencies, f0, 128)
            quadrature_difference = float(max(abs(raw-check)))
            chosen = np.argsort(-abs(raw), kind='stable')[:N]
            threshold = float(min(abs(raw[chosen])))
            detuning = (jmax+1)/T-f0
            tail_bound = 80/(np.sqrt(3)*R**2*T*(2*np.pi*detuning)**3)
            if tail_bound+2*quadrature_difference < threshold:
                break
            jmax *= 2
        else:
            raise RuntimeError('Candidate expansion could not separate top40 from omitted-frequency bound')
        ranked = chosen.copy()
        for column, budget in enumerate((10, 20, 40)):
            chosen = ranked[:budget]
            chosen = chosen[np.argsort(frequencies[chosen])]
            selected_f, selected_raw = frequencies[chosen], raw[chosen]
            raw_anchor = float(sum(selected_raw))
            if not np.isfinite(raw_anchor) or abs(raw_anchor) < 1e-10:
                raise RuntimeError('Top40 midpoint sum cannot be safely normalized')
            scale = 1/raw_anchor
            midpoint_coefficients = selected_raw*scale*np.exp(1j*phase)
            z = evaluate(u, selected_f, midpoint_coefficients)
            g = z*np.exp(-1j*phase-2j*np.pi*f0*u)
            target_z = target_envelope*np.exp(1j*phase+2j*np.pi*f0*u)
            curvature = float(sum(abs(midpoint_coefficients)*(2*np.pi*(selected_f-f0))**2))
            allowance = curvature*.005**2/8
            dc = float(coefficients(np.array([0.]), f0, 256)[0])
            negative = coefficients(-frequencies, f0, 256)
            positive_energy = float(sum(raw**2)); negative_energy = float(sum(negative**2))
            selected_energy = float(sum(selected_raw**2))
            dwell_error, quiet_error = float(max(abs(g[dwell]-1))), float(max(abs(g[quiet])))
            ramp_max = float(max(abs(g[ramp])))
            report = dict(carrier_frequency_hz=f0, component_budget=budget, phase_at_midpoint_rad=phase,
                phase_reference_time_s=CENTER, candidate_positive_components=jmax,
                selected_positive_components=len(chosen), nonzero_components_above_1e_12=int(sum(abs(midpoint_coefficients)>1e-12)),
                component_frequencies_hz=selected_f.tolist(), raw_fourier_coefficients=selected_raw.tolist(),
                midpoint_coefficients_real=midpoint_coefficients.real.tolist(), midpoint_coefficients_imag=midpoint_coefficients.imag.tolist(),
                raw_midpoint_sum=raw_anchor, normalization_scale=scale,
                midpoint_anchor_error=float(abs(sum(midpoint_coefficients)-np.exp(1j*phase))),
                candidate_cutoff_hz=jmax/T, omitted_positive_coefficient_bound=tail_bound,
                smallest_selected_raw_magnitude=float(min(abs(selected_raw))), top40_threshold=threshold, quadrature_128_256_max_difference=quadrature_difference,
                top40_tail_separated=True, tail_separation_is_interval_proof=False,
                analytic_target_mean_square=total_energy, selected_raw_positive_energy=selected_energy,
                sampled_dictionary_positive_energy=positive_energy, sampled_dictionary_negative_energy=negative_energy,
                dc_energy=dc**2, signed_dc_template_coefficient=dc,
                target_real_mean=dc*float(np.cos(phase)), reconstructed_real_mean=0.,
                reconstructed_minus_target_real_mean=-dc*float(np.cos(phase)),
                total_unselected_raw_energy_parseval=total_energy-selected_energy,
                unaccounted_energy_beyond_signed_dictionary=total_energy-positive_energy-negative_energy-dc**2,
                normalized_selected_approximation_mean_square_error=total_energy+(scale**2-2*scale)*selected_energy,
                dwell_max_complex_error_observed=dwell_error, quiet_max_amplitude_observed=quiet_error,
                ramp_max_amplitude_observed=ramp_max, curvature_allowance=allowance,
                worst_score_continuous_upper=max(dwell_error, quiet_error, ramp_max-1)+allowance,
                real_target_max_error_observed=float(max(abs(z.real-target_z.real))),
                full_template_complex_max_error_observed=float(max(abs(z-target_z))),
                envelope_magnitude_mirror_residual=float(max(abs(abs(g)-abs(g[::-1])))),
                envelope_conjugate_mirror_residual=float(max(abs(g-np.conj(g[::-1])))))
            report.update(raw_dwell_max_complex_error=float(max(abs(g[dwell]/scale-1))),
                raw_quiet_max_amplitude=float(max(abs(g[quiet]/scale))),
                raw_full_template_complex_max_error=float(max(abs(z/scale-target_z))),
                raw_positive_truncation_energy=positive_energy-selected_energy)
            reports.append(report)
            arrays.update({f'frequency_{i}_{budget}_hz': selected_f, f'raw_coefficients_{i}_{budget}': selected_raw,
                f'midpoint_coefficients_{i}_{budget}': midpoint_coefficients, f'analytic_wave_{i}_{budget}': z,
                f'complex_envelope_{i}_{budget}': g, f'target_analytic_wave_{i}_{budget}': target_z})
            axes[i, column].plot(t, z.real, lw=.7, label=f'{budget}-line real elevation')
            axes[i, column].plot(t, abs(g), color='black', lw=1, label='Reconstructed envelope magnitude')
            axes[i, column].plot(t, target_envelope, color='green', ls=':', lw=1, label='Fixed quintic envelope')
            axes[i, column].plot(t[dwell], target_z.real[dwell], color='red', ls='--', lw=.65, label='Desired carrier during dwell')
            for boundary in (10, 15, 115, 120): axes[i, column].axvline(boundary, color='grey', ls=':', lw=.5)
            axes[i, column].set(title=f'f0={f0:g} Hz; N={budget}; midpoint phase={np.degrees(phase):g} deg', xlabel='Time (s)', ylabel='Elevation (m) / envelope')
            axes[i, column].legend(fontsize=6)
            specaxes[i, column].bar(selected_f, abs(midpoint_coefficients), width=.003)
            specaxes[i, column].axvline(f0, color='black', ls='--', lw=.7)
            specaxes[i, column].set(title=f'f0={f0:g} Hz; N={budget}; normalized amplitudes', xlabel='Physical frequency (Hz)', ylabel='Cosine amplitude (m)')
            phaseaxes[i, column].bar(selected_f, np.angle(midpoint_coefficients, deg=True), width=.003)
            phaseaxes[i, column].set(title=f'f0={f0:g} Hz; N={budget}; phase at65 s', xlabel='Physical frequency (Hz)', ylabel='Degrees', ylim=(-180, 180))
    fig.suptitle('Fixed quintic + direct Fourier + largest N POSITIVE components; no optimizer')
    fig.savefig(out/'pulse40-three-quintic.png', dpi=160); plt.close(fig)
    specfig.suptitle('Discrete component amplitudes; N=10 /20 /40 columns')
    specfig.savefig(out/'pulse40-three-quintic-spectrum.png', dpi=160); plt.close(specfig)
    phasefig.suptitle('Component phases defined directly at midpoint65 s')
    phasefig.savefig(out/'pulse40-three-quintic-phase.png', dpi=160); plt.close(phasefig)
    np.savez(out/'packets.npz', **arrays)
    report = dict(policy='Fixed quintic, direct Fourier, nested top10/20/40 positive magnitudes, separate midpoint normalization',
        timing_s=dict(rest=10, ramp=5, dwell=100, end=130, phase_reference=65),
        harmonic_spacing_hz=1/T, cases=reports, elapsed_s=time.perf_counter()-started,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limitations=['No DC or negative physical frequencies retained or folded.',
            'Positive-frequency projection of COMPLEX target a*exp(i*omega0*u+i*phi), not Fourier truncation of its real part; real-wave truncation would combine signed-frequency template coefficients.',
            'TopN unscaled coefficients minimize periodic mean-square error within positive Fourier basis, not peak error.',
            'Midpoint normalization changes that unscaled projection; peak error need not improve monotonically withN.',
            '128/256 quadrature agreement is diagnostic, not rigorous interval certification.',
            'No physical accuracy gate, propagation or hydrodynamics; only0–130s observation; Fourier construction repeats every130s.'])
    (out/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')


if __name__ == '__main__': main()
