"""Harmonic ranking and finite-kernel diagnostics only; no time-domain solver."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.interpolate import interp1d
from scipy.io import loadmat, savemat
from threadpoolctl import threadpool_limits


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    started = time.monotonic()
    native = Path('artifacts/reference/generated-bank-native-20261024/attempt1')
    consumer = Path('artifacts/reference/generated-bank-consumer-20261024/bemio-attempt1')
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    manifest = read(native/'manifest.json')
    assert manifest['complete']
    for name, digest in manifest['output_sha256'].items():
        assert sha(native/name) == digest, name
    n = loadmat(native/'native.mat')
    rad = loadmat(native/'radiation.mat')
    h = loadmat(consumer/'bemio.mat', simplify_cells=True)['h']
    w = rad['omega_rad_s'].ravel()
    rho = n['rho'].item()
    K = n['K']
    M = n['M']
    Ainf = n['Ainf']
    h5check = read(consumer/'h5-readback.json')
    assert h5check['complete']
    kernel = np.moveaxis(h['ra_K'], -1, 0)*rho
    lags = np.asarray(h['ra_t']).ravel()
    assert kernel.shape == (2001, 6, 6)
    assert np.allclose(lags, np.arange(2001)*.02, atol=1e-12, rtol=0)
    times = np.arange(3001)*.02
    D = np.diag([2.7432]*3+[1.]*3)
    records = []
    synthetic_errors = []
    x = np.array([0., np.pi/2, np.pi, 3*np.pi/2])
    scales = np.arange(1, 7)
    for coefficient, phase, expected in ((1j, 0., np.sin(x)),
                                        (1., np.pi/2, np.sin(x)),
                                        (1j, np.pi/2, -np.cos(x))):
        actual = np.real(np.exp(1j*(phase-x))[:, None]*(coefficient*scales))
        error = np.max(abs(actual-expected[:, None]*scales))
        assert error < 1e-13
        synthetic_errors.append(float(error))
    with threadpool_limits(1):
        for pair_path in sorted(native.glob('*-pairs.mat')):
            values = loadmat(pair_path)
            omega = values['omega_rad_s'].ravel()
            a = values['amplitude_m'].ravel()
            phase = values['phase_local_rad'].ravel()
            H = values['excitation_complex']
            assert np.array_equal(values['component_index'].ravel(), np.arange(1, 385))
            assert np.array_equal(H, values['fk_complex']+values['diffraction_complex'])
            assert omega.min() >= w.min() and omega.max() <= w.max()
            A = interp1d(w, rad['added_mass'], axis=0, bounds_error=True)(omega)
            B = interp1d(w, rad['radiation_damping'], axis=0, bounds_error=True)(omega)
            harmonic = a[:, None]*np.exp(1j*phase[:, None])*H
            Q = np.empty_like(harmonic)
            condition = np.empty(len(omega))
            for j in range(len(omega)):
                Z = K-omega[j]**2*(M+A[j])-1j*omega[j]*B[j]
                condition[j] = np.linalg.cond(D@Z@D)
                assert np.isfinite(condition[j]) and condition[j] < 1e12, 'Ill-conditioned harmonic solve'
                Q[j] = D@np.linalg.solve(D@Z@D, D@harmonic[j])
            exponent = np.exp(-1j*times[:, None]*omega)
            force = np.real(exponent@harmonic)
            response = np.real(exponent@Q)
            scalar = np.array([[sum(a[j]*(H[j, dof].real*np.cos(phase[j]-omega[j]*t)
                       -H[j, dof].imag*np.sin(phase[j]-omega[j]*t)) for j in range(384))
                       for dof in range(6)] for t in times[[0, 19, 157, -1]]])
            scalar_error = float(np.max(abs(force[[0, 19, 157, -1]]-scalar)))
            np.testing.assert_allclose(force[[0, 19, 157, -1]], scalar, atol=1e-8, rtol=2e-12)
            j = int(np.flatnonzero(a > 0)[0])
            quarter = np.arange(4)*np.pi/(2*omega[j])
            phi = phase[j]-omega[j]*quarter
            independent = a[j]*(np.cos(phi)[:, None]*H[j].real-np.sin(phi)[:, None]*H[j].imag)
            direct = np.real(a[j]*np.exp(1j*phi[:, None])*H[j])
            np.testing.assert_allclose(direct, independent, atol=1e-9, rtol=2e-12)
            # Finite causal kernel convention: B=integral K cos; A=Ainf-integral K sin/omega.
            Brecon = np.empty_like(B)
            Arecon = np.empty_like(A)
            for j, frequency in enumerate(omega):
                Brecon[j] = np.trapezoid(kernel*np.cos(frequency*lags)[:, None, None], lags, axis=0)
                Arecon[j] = Ainf-np.trapezoid(kernel*np.sin(frequency*lags)[:, None, None], lags, axis=0)/frequency
            Zerror = -omega[:, None, None]**2*(Arecon-A)-1j*omega[:, None, None]*(Brecon-B)
            record_id = pair_path.name.removesuffix('-pairs.mat')
            # This file is a prepared native harmonic force table, not a trajectory.
            savemat(out/(record_id+'-harmonic.mat'), dict(time=times[:, None], force=force,
                stationary_response=response, component_response_complex=Q,
                component_force_complex=harmonic, A_screen_interpolated=A, B_screen_interpolated=B,
                condition_scaled_impedance=condition[:, None]))
            records.append(dict(id=record_id, native_pair_sha256=sha(pair_path),
                harmonic_abs_peak_m_rad=np.max(abs(response), axis=0).tolist(),
                harmonic_peak_to_peak_m_rad=np.ptp(response, axis=0).tolist(),
                harmonic_rms_m_rad=np.sqrt(np.mean(response**2, axis=0)).tolist(),
                harmonic_roll_peak_deg=float(np.rad2deg(np.max(abs(response[:, 3])))),
                force_abs_peak_N_Nm=np.max(abs(force), axis=0).tolist(),
                scaled_impedance_condition_max=float(condition.max()),
                scalar_vector_force_max_error=scalar_error,
                actual_harmonic_quarter_cycle_max_error=float(np.max(abs(direct-independent))),
                reconstructed_A_abs_error_per_entry=np.max(abs(Arecon-A), axis=0).tolist(),
                reconstructed_B_abs_error_per_entry=np.max(abs(Brecon-B), axis=0).tolist(),
                reconstructed_Z_abs_error_per_entry=np.max(abs(Zerror), axis=0).tolist()))
    total = np.trapezoid(abs(kernel), lags, axis=0)
    tail = np.trapezoid(abs(kernel[-251:]), lags[-251:], axis=0)
    ratio = np.divide(tail, total, out=np.zeros_like(total), where=total != 0)
    report = dict(complete=True, elapsed_seconds=time.monotonic()-started,
        scope='Stationary harmonic ranking ONLY; no ODE, CIC simulation, initialization, drift or linearity qualification',
        convention='Native exp(-i omega t); no extra conjugation or CG spatial phase',
        radiation_memory_s=40, dt_s=.02, harmonic_window_s=[0, 60],
        radiation_interpolation='Componentwise linear omega interpolation of all36 native A/B entries; screen only',
        synthetic_sign_max_errors=synthetic_errors,
        kernel_abs_last5s_fraction_per_entry=ratio.tolist(),
        kernel_endpoint_abs_per_entry=abs(kernel[-1]).tolist(),
        radiation_B_low_endpoint=rad['radiation_damping'][0].tolist(),
        radiation_B_high_endpoint=rad['radiation_damping'][-1].tolist(), records=records,
        ranking_descending_harmonic_roll=[r['id'] for r in sorted(records, key=lambda r:-r['harmonic_roll_peak_deg'])],
        limitations=['Zero-restoring surge/sway/yaw homogeneous startup terms and drift unresolved',
                    'Harmonic response is not the held-release initial-value trajectory',
                    'Initial immersion, relative-wave geometry and curved-surface effects unqualified',
                    '17panel mesh/frequency cutoffs/kernel memory/temporal uncertainty unassessed'],
        source_sha256=sha(__file__), bemio_sha256=sha(consumer/'bemio.mat'),
        h5_readback=h5check,
        output_sha256={path.name:sha(path) for path in out.glob('*.mat')})
    (out/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(complete=True, records=records, ranking=report['ranking_descending_harmonic_roll'])))


if __name__ == '__main__':
    main()
