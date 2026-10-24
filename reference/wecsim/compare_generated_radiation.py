"""Stored-data harmonic radiation comparison; no BEM, ODE or body simulation."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.interpolate import interp1d
from scipy.io import loadmat, savemat
from threadpoolctl import threadpool_limits
from screen_generated_banks import read, sha


def weights(nu, omega, duration):
    minus, plus = nu-omega, nu+omega
    C = .5*duration*(np.sinc(minus*duration/np.pi)+np.sinc(plus*duration/np.pi))
    S = .25*duration**2*(plus*np.sinc(plus*duration/(2*np.pi))**2
                           -minus*np.sinc(minus*duration/(2*np.pi))**2)
    return C, S


def pv_band(nu, B, omega, Bomega, derivative):
    denominator = nu**2-omega**2
    near = abs(nu-omega) < 1e-12
    integrand = np.empty_like(B)
    integrand[~near] = (B[~near]-Bomega)/denominator[~near, None, None]
    integrand[near] = derivative/(2*omega)
    lo, hi = nu[[0, -1]]
    constant = np.log(abs((hi-omega)*(lo+omega)/((hi+omega)*(lo-omega))))/(2*omega)
    return np.trapezoid(integrand, nu, axis=0)+Bomega*constant


def stats(history):
    # Translations remain metres; every angular column is converted to degrees.
    scale = np.array([1.]*3+[180/np.pi]*3)
    x = history*scale
    return dict(abs_peak_m_deg=np.max(abs(x), axis=0).tolist(),
                peak_to_peak_m_deg=np.ptp(x, axis=0).tolist(),
                rms_m_deg=np.sqrt(np.mean(x*x, axis=0)).tolist())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    started = time.monotonic()
    native = Path('artifacts/reference/generated-bank-native-20261024/attempt1')
    consumer = Path('artifacts/reference/generated-bank-consumer-20261024/bemio-attempt1')
    incident = Path('artifacts/reference/incident-input-diagnostics-20261024')
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    n = loadmat(native/'native.mat')
    rad = loadmat(native/'radiation.mat')
    h = loadmat(consumer/'bemio.mat', simplify_cells=True)['h']
    manifest = read(native/'manifest.json')
    for name, digest in manifest['output_sha256'].items():
        assert sha(native/name) == digest
    rho, gravity = n['rho'].item(), n['g'].item()
    w = rad['omega_rad_s'].ravel()
    np.testing.assert_allclose(np.moveaxis(h['A'], -1, 0)*rho, rad['added_mass'], atol=1e-10, rtol=2e-14)
    np.testing.assert_allclose(np.moveaxis(h['B'], -1, 0)*rho*w[:, None, None], rad['radiation_damping'], atol=1e-10, rtol=2e-14)
    np.testing.assert_allclose(h['Ainf']*rho, n['Ainf'], atol=1e-10, rtol=2e-14)
    np.testing.assert_allclose(h['Khs']*rho*gravity, n['K'], atol=1e-10, rtol=2e-14)
    nu = np.asarray(h['ra_w']).ravel()
    t = np.asarray(h['ra_t']).ravel()
    assert abs(t[0]) < 1e-14 and abs(t[-1]-40.) < 1e-12
    np.testing.assert_allclose(t, np.arange(2001)*.02, atol=1e-12, rtol=0)
    assert abs(nu[0]-.05) < 1e-14 and abs(nu[-1]-12.) < 1e-12
    kernel = np.moveaxis(h['ra_K'], -1, 0)*rho
    normalized_B = np.moveaxis(h['B'], -1, 0)
    interp_Bnorm = interp1d(w, normalized_B, axis=0, bounds_error=True)
    Bstar = interp_Bnorm(nu)*rho*nu[:, None, None]
    Blinear = interp1d(w, rad['radiation_damping'], axis=0)(nu)
    # Analytic constant-band PV fixture and removable-limit fixtures.
    test_omega, T = .73, 40.
    constant = np.ones((len(nu), 6, 6))*np.arange(1, 37).reshape(6, 6)
    lo, hi = nu[[0, -1]]
    expected = constant[0]/(2*test_omega)*np.log(abs((hi-test_omega)*(lo+test_omega)/((hi+test_omega)*(lo-test_omega))))
    np.testing.assert_allclose(pv_band(nu, constant, test_omega, constant[0], np.zeros((6, 6))), expected, atol=1e-12, rtol=2e-14)
    C, S = weights(np.array([test_omega]), test_omega, T)
    np.testing.assert_allclose(C[0], .5*(T+np.sin(2*test_omega*T)/(2*test_omega)), atol=1e-13)
    np.testing.assert_allclose(S[0], (1-np.cos(2*test_omega*T))/(4*test_omega), atol=1e-13)
    # Joint plane angles for Ry(pitch)Rx(roll); its third column is the plane normal.
    for sx, sy in [(0., 0.), (.3, -.4), (-.7, .2)]:
        pitch = -np.arctan(sx)
        roll = np.arctan(sy/np.sqrt(1+sx*sx))
        normal = np.array([np.sin(pitch)*np.cos(roll), -np.sin(roll), np.cos(pitch)*np.cos(roll)])
        expected_normal = np.array([-sx, -sy, 1.])/np.sqrt(1+sx*sx+sy*sy)
        np.testing.assert_allclose(normal, expected_normal, atol=3e-15, rtol=0)
    times = np.arange(3001)*.02
    D = np.diag([2.7432]*3+[1.]*3)
    records = []
    transform_cache = {}
    with threadpool_limits(1):
        for pair_path in sorted(native.glob('*-pairs.mat')):
            values = loadmat(pair_path)
            omega = values['omega_rad_s'].ravel()
            a = values['amplitude_m'].ravel()
            phase = values['phase_local_rad'].ravel()
            force = a[:, None]*np.exp(1j*phase[:, None])*values['excitation_complex']
            A = interp1d(w, rad['added_mass'], axis=0)(omega)
            B = interp1d(w, rad['radiation_damping'], axis=0)(omega)
            AK, BK, AT, BT, AB = [np.empty_like(A) for _ in range(5)]
            Qnative, Qkernel = np.empty_like(force), np.empty_like(force)
            condition = np.empty((len(omega), 2))
            residual = np.zeros(2)
            for j, frequency in enumerate(omega):
                if time.monotonic()-started > 120:
                    raise TimeoutError('Prospective120s comparison budget exceeded')
                if frequency not in transform_cache:
                    cosine = np.cos(frequency*t)
                    sine = np.sin(frequency*t)
                    bk = np.trapezoid(kernel*cosine[:, None, None], t, axis=0)
                    ak = n['Ainf']-np.trapezoid(kernel*sine[:, None, None], t, axis=0)/frequency
                    C, S = weights(nu, frequency, T)
                    bt = (2/np.pi)*np.trapezoid(Bstar*C[:, None, None], nu, axis=0)
                    at = n['Ainf']-(2/(np.pi*frequency))*np.trapezoid(Bstar*S[:, None, None], nu, axis=0)
                    i = min(len(w)-2, max(0, np.searchsorted(w, frequency)-1))
                    slope = (normalized_B[i+1]-normalized_B[i])/(w[i+1]-w[i])
                    bomega = interp_Bnorm(frequency)*rho*frequency
                    derivative = rho*(interp_Bnorm(frequency)+frequency*slope)
                    ab = n['Ainf']+(2/np.pi)*pv_band(nu, Bstar, frequency, bomega, derivative)
                    transform_cache[frequency] = ak, bk, at, bt, ab
                AK[j], BK[j], AT[j], BT[j], AB[j] = transform_cache[frequency]
                for model, (am, damping, Q) in enumerate(((A[j], B[j], Qnative), (AK[j], BK[j], Qkernel))):
                    Z = n['K']-frequency**2*(n['M']+am)-1j*frequency*damping
                    condition[j, model] = np.linalg.cond(D@Z@D)
                    assert np.isfinite(condition[j, model]) and condition[j, model] < 1e12
                    Q[j] = D@np.linalg.solve(D@Z@D, D@force[j])
                    residual[model] = max(residual[model], np.max(abs(Z@Q[j]-force[j])))
            exponent = np.exp(-1j*times[:, None]*omega)
            qnative = np.real(exponent@Qnative)
            qkernel = np.real(exponent@Qkernel)
            difference = qkernel-qnative
            native_energy = np.sum(qnative*qnative, axis=0)
            difference_energy = np.sum(difference*difference, axis=0)
            ratio = [float(np.sqrt(d/n)) if n > 0 else None for d, n in zip(difference_energy, native_energy)]
            record_id = pair_path.name.removesuffix('-pairs.mat')
            with np.load(incident/(record_id+'-offset-0deg.npz')) as source:
                np.testing.assert_allclose(source['time_s'], times, atol=1e-12, rtol=0)
                fields = source['fields']
                eta = fields[0, 0]
                rollplane = np.arctan(fields[2, 0]/np.sqrt(1+fields[1, 0]**2))
                pitchplane = -np.arctan(fields[1, 0])
                relative = {}
                for name, q in [('native', qnative), ('kernel40', qkernel)]:
                    rel = np.column_stack((q[:, 2]-eta, q[:, 3]-rollplane, q[:, 4]-pitchplane))
                    relative[name] = dict(heave_roll_pitch_abs_peak_m_deg=(np.max(abs(rel), axis=0)*[1,180/np.pi,180/np.pi]).tolist(),
                        particular_t0_q_m_rad=q[0].tolist(), particular_t0_relative_heave_roll_pitch_m_rad=rel[0].tolist())
                field_summary = dict(meanCG_t0_eta_slopes_Hessian_etat=fields[:, 0, 0].tolist(),
                    held_mean_pose_t0_relative_heave_roll_pitch_m_rad=[float(-eta[0]), float(-rollplane[0]), float(-pitchplane[0])],
                    spatial_departure_abs_t0_m=float(np.max(abs(source['departure_from_cg_plane_m'][:, 0]))),
                    spatial_departure_abs_peak_m=float(np.max(abs(source['departure_from_cg_plane_m']))),
                    cg_hessian_abs_peak_per_m=np.max(abs(fields[3:6, 0]), axis=1).tolist(),
                    caveat='FixedmeanCG/footprint plane proxies only; Euler angle subtraction is a screening proxy, not exact finite relative orientation; no advected/rotated hull or curvedsurface wetset test; no actualIC selected')
            savemat(out/(record_id+'-comparison.mat'), dict(time=times[:, None], native_particular=qnative,
                kernel40_particular=qkernel, difference=difference, A_native_interpolated=A,
                B_native_interpolated=B, A_kernel40=AK, B_kernel40=BK, A_exacttime40=AT,
                B_exacttime40=BT, A_infinite_time_finiteband=AB,
                Ainf_effective_native_plus_kernel_sine=A+n['Ainf']-AK,
                omega_rad_s=omega[:, None], Q_native=Qnative, Q_kernel40=Qkernel))
            records.append(dict(id=record_id, native=stats(qnative), kernel40=stats(qkernel),
                difference=stats(difference), rms_difference_over_native_rms=ratio,
                zero_native_energy_axes=np.flatnonzero(native_energy == 0).tolist(),
                scaled_condition_max=np.max(condition, axis=0).tolist(),
                ZQ_minus_force_abs_max=residual.tolist(),
                exacttime_vs_discrete_A_abs_error_per_entry=np.max(abs(AT-AK), axis=0).tolist(),
                exacttime_vs_discrete_B_abs_error_per_entry=np.max(abs(BT-BK), axis=0).tolist(),
                finiteband_vs_native_A_abs_error_per_entry=np.max(abs(AB-A), axis=0).tolist(),
                fixedmean_wave_relative_screen=relative, incident_fields=field_summary))
    report = dict(complete=True, elapsed_seconds=time.monotonic()-started, records=records,
        A_at_12_minus_direct_Ainf=(rad['added_mass'][-1]-n['Ainf']).tolist(),
        Bstar_vs_linear_dimensional_B_abs_error_per_entry=np.max(abs(Bstar-Blinear), axis=0).tolist(),
        scope='Storeddata particularresponse diagnostic ONLY; no actualbody/ODE/IC/drift/linearity qualification',
        normalization_checks='h.A*rho, h.B*rho*w, h.Ainf*rho, h.Khs*rho*g checked againstnative',
        fixtures='Constantband PV, exact removable C/S limits and joint RyRx plane-normal closure pass',
        virtualmass_delay='Excluded: actualfeedbackpath needs separate reviewedsource derivation',
        limitations=['FinitebandKK mismatch cannot separate cutofftails fromcoarsemesh/nativeinconsistency',
                    'Stationaryparticularhistories omit startup/free-mode drift',
                    'Fixedmeanwave-relativeproxies do not establish physicalimmersion/wetset/linearity'],
        source_sha256=sha(__file__), output_sha256={p.name:sha(p) for p in out.glob('*.mat')})
    (out/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(complete=True, elapsed_seconds=report['elapsed_seconds'],
                         records=[dict(id=r['id'], native=r['native'], kernel40=r['kernel40'], difference=r['difference']) for r in records])))


if __name__ == '__main__':
    main()
