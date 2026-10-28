"""Real quintic target -> FFT -> physical-frequency cutoff -> IFFT. No fit."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

T, CENTER, M = 130., 65., 32768
F0S, PHASES, CUTOFFS = [.01, .1, 1.], [0., np.pi/4, np.pi/2], [1.05, 1.20, 1.50]


def envelope(u):
    x = np.clip((abs(u)-50)/5, 0, 1)
    return 1-x**3*(10-15*x+6*x*x)


def mask_for(size, cutoff):
    jmax = int(np.floor(T*cutoff+1e-10))
    if jmax >= size//2: raise RuntimeError('Requested cutoff reaches numerical Nyquist')
    mask = np.zeros(size, dtype=bool)
    mask[:jmax+1] = True; mask[-jmax:] = True
    return mask, jmax


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    out = Path(parser.parse_args().output); out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    t, fine_t = np.arange(M)*T/M, np.arange(2*M)*T/(2*M)
    u, fine_u = t-CENTER, fine_t-CENTER
    a, fine_a = envelope(u), envelope(fine_u)
    quiet, dwell = abs(u) >= 55, abs(u) <= 50
    arrays, reports = {'time_s':t,'target_envelope':a}, []
    fig, axes = plt.subplots(3,3,figsize=(17,11),layout='constrained')
    specfig, specaxes = plt.subplots(3,3,figsize=(17,10),layout='constrained')
    errfig, erraxes = plt.subplots(3,3,figsize=(17,10),layout='constrained')
    for i,(f0,phase) in enumerate(zip(F0S,PHASES)):
        target = a*np.cos(2*np.pi*f0*u+phase)
        fine_target = fine_a*np.cos(2*np.pi*f0*fine_u+phase)
        full_fft, fine_fft = np.fft.fft(target)/M, np.fft.fft(fine_target)/(2*M)
        roundtrip_error = float(max(abs(np.fft.ifft(full_fft*M).real-target)))
        if roundtrip_error > 1e-10: raise RuntimeError('Target Fourier roundtrip failed')
        for column,cutoff in enumerate(CUTOFFS):
            mask,jmax = mask_for(M,cutoff); fine_mask,fine_jmax = mask_for(2*M,cutoff)
            filtered = np.where(mask,full_fft,0.)
            wave_complex = np.fft.ifft(filtered*M)
            fine_filtered = np.where(fine_mask,fine_fft,0.)
            fine_wave = np.fft.ifft(fine_filtered*(2*M))
            wave, error = wave_complex.real, wave_complex.real-target
            pos = np.arange(1,jmax+1); sf = pos/T
            midpoint_cosine_coeff = 2*filtered[pos]*np.exp(2j*np.pi*sf*CENTER)
            imag_residue = float(max(abs(wave_complex.imag)))
            hermitian_residue = float(max(abs(filtered[pos]-np.conj(filtered[-pos]))))
            mean_error = float(abs(np.mean(wave)-np.mean(target)))
            if (not np.all(np.isfinite(wave)) or imag_residue > 1e-10
                    or hermitian_residue > 1e-10 or mean_error > 1e-10):
                raise RuntimeError('Finite/Hermitian/realness/mean identity check failed')
            if fine_jmax != jmax or int(sum(mask)) != 2*jmax+1:
                raise RuntimeError('Cutoff mask count mismatch')
            report = dict(carrier_frequency_hz=f0,phase_at_midpoint_rad=phase,
                phase_reference_time_s=CENTER,maximum_physical_frequency_hz=cutoff,
                actual_maximum_retained_frequency_hz=jmax/T,positive_bin_count=jmax,
                negative_conjugate_bin_count=jmax,dc_bin_count=1,total_retained_bins=2*jmax+1,
                no_holes=True,normalization_applied=False,target_mean=float(np.mean(target)),
                reconstructed_mean=float(np.mean(wave)),mean_error=mean_error,
                imaginary_residue=imag_residue,hermitian_coefficient_residue=hermitian_residue,
                target_fft_roundtrip_error=roundtrip_error,
                peak_real_error_observed=float(max(abs(error))),rms_real_error=float(np.sqrt(np.mean(error**2))),
                dwell_peak_real_error_observed=float(max(abs(error[dwell]))),
                quiet_peak_real_elevation_observed=float(max(abs(wave[quiet]))),
                finer_fft_retained_coefficient_max_difference=float(max(abs(full_fft[mask]-fine_fft[np.r_[np.arange(jmax+1),np.arange(2*M-jmax,2*M)]]))),
                finer_ifft_common_grid_max_difference=float(max(abs(wave-fine_wave.real[::2]))),
                finer_ifft_peak_real_error_observed=float(max(abs(fine_wave.real-fine_target))),
                discarded_sampled_fourier_energy=float(sum(abs(full_fft[~mask])**2)))
            reports.append(report)
            arrays.update({f'target_{i}_{column}_m':target,f'elevation_{i}_{column}_m':wave,
                f'filtered_fft_{i}_{column}':filtered,f'mask_{i}_{column}':mask,
                f'positive_frequency_{i}_{column}_hz':sf,f'midpoint_cosine_coefficients_{i}_{column}':midpoint_cosine_coeff})
            axes[i,column].plot(t,wave,lw=.8,label='Cutoff reconstruction')
            axes[i,column].plot(t,target,color='red',ls='--',lw=.65,label='Real quintic target')
            axes[i,column].plot(t,a,color='green',ls=':',lw=.7,label='Prescribed envelope guides')
            axes[i,column].plot(t,-a,color='green',ls=':',lw=.7)
            for boundary in (10,15,115,120): axes[i,column].axvline(boundary,color='grey',ls=':',lw=.5)
            title=f'f0={f0:g} Hz; fmax={cutoff:g} Hz; phase={np.degrees(phase):g} deg'
            axes[i,column].set(title=title,xlabel='Time (s)',ylabel='Real elevation (m)')
            axes[i,column].legend(fontsize=6)
            specaxes[i,column].bar(sf,abs(midpoint_cosine_coeff),width=.004)
            specaxes[i,column].set(title=title,xlabel='Positive physical frequency (Hz)',ylabel='Cosine amplitude (m)')
            erraxes[i,column].plot(t,error,lw=.8)
            erraxes[i,column].set(title=title,xlabel='Time (s)',ylabel='Reconstruction minus target (m)')
    for figure,name,title in ((fig,'pulse-fmax-time.png','Fixed quintic real waves; physical fmax cutoff; no optimizer or rescaling'),
            (specfig,'pulse-fmax-spectrum.png','All positive bins through fmax retained; DC and negative conjugates also retained'),
            (errfig,'pulse-fmax-error.png','Pointwise truncation error; rows=carrier; columns=fmax')):
        figure.suptitle(title); figure.savefig(out/name,dpi=160); plt.close(figure)
    np.savez(out/'packets.npz',**arrays)
    report=dict(policy='Sample real quintic -> FFT -> abs(physical frequency)<=fmax -> IFFT',
        timing_s=dict(rest=10,ramp=5,dwell=100,end=130,phase_reference=65),
        samples=M,reference_samples=2*M,period_s=T,frequency_spacing_hz=1/T,cases=reports,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_s=time.perf_counter()-started,
        limitations=['No physical accuracy acceptance; errors are sampled numerical truncation diagnostics.',
            '32768/65536 comparison measures sampling sensitivity, not an interval or continuum proof.',
            'The Fourier representation repeats every130s; no statement about later isolated calm.',
            'Finite cutoff need not preserve target midpoint value exactly; input phase defines target directly at65s.'])
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__ == '__main__': main()
