"""One fixed physical train preflight: initial sixforce screen BEFORE efficient hull bounds."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.io import savemat
from scipy.signal import periodogram
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits
from calm_buoy_train import CalmBuoyTrain, window, window_spectrum


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    started = time.monotonic()
    controls_path = Path('reference/wecsim/calm-buoy-train-controls.json')
    controls = json.loads(controls_path.read_text())
    source = Path(controls['carrier_export'])
    case_path = source/'case.json'
    carriers_path = source/'components.npz'
    observation_path = source/'observation.json'
    case = json.loads(case_path.read_text())
    observed = json.loads(observation_path.read_text())
    physical_path = Path('reference/capytaine/canonical-boat.json')
    physical = json.loads(physical_path.read_text())
    g = physical['water']['gravity_m_s2']
    beam = 9*.3048
    with np.load(carriers_path) as archive:
        carriers = {key:archive[key].copy() for key in archive.files}
    assert case['observation_id'] == controls['observation_id'] and case['seed'] == 42
    assert np.array_equal(carriers['component_id'],np.arange(384))
    assert np.max(abs(window_spectrum(np.array([0.]))-20/np.pi)) < 1e-12
    assert np.array_equal(window(np.array([0.,15.,25.,30.])),np.array([1.,1.,0.,0.]))
    data = CalmBuoyTrain.construct(carriers,case,controls,g,beam)
    args.output.mkdir(parents=True,exist_ok=False)
    np.savez(args.output/'train.npz',**data)
    train = CalmBuoyTrain(args.output/'train.npz')
    table_path = Path('artifacts/reference/calm-buoy-20261024/excitation-attempt1/excitation.npz')
    with np.load(table_path) as archive:
        table = {key:archive[key].copy() for key in archive.files}
    assert str(table['bank_case_sha256'].item()) == sha(case_path)
    assert str(table['bank_components_sha256'].item()) == sha(carriers_path)
    assert np.array_equal(table['global_from_native_R'],np.asarray(case['global_from_boat_rotation']))
    assert np.array_equal(table['global_from_native_t'],np.asarray(case['global_from_boat_translation_m']))
    times = np.arange(3001)*.02
    force,scale,spectral_force = train.force(table,times)
    H0,H0_error = train.static_excitation(table)
    assert np.all(np.isfinite(force))
    limit = np.minimum(.1,.001*scale)
    force_pass = bool(np.all(abs(force[0]) <= limit))
    coefficients = train.C.sum(axis=(0,1))
    eta = (np.exp(-1j*times[:,None]*train.omega)@coefficients).real
    active = np.flatnonzero(carriers['amplitudes'])
    component_complex = train.C[active].sum(axis=1)@np.exp(-1j*times[:,None]*train.omega).T
    per_carrier_arrival = np.zeros(384)
    per_carrier_arrival[active] = abs(component_complex).max(axis=1)/carriers['amplitudes'][active]
    combined_arrival = np.sqrt(np.sum(abs(component_complex)**2,axis=0)/np.sum(carriers['amplitudes']**2))
    arrival_pass = bool(combined_arrival.max() >= controls['minimum_combined_arrival_fraction'])
    arrival = dict(combined_peak_fraction=float(combined_arrival.max()),
        combined_peak_time_s=float(times[combined_arrival.argmax()]),
        per_carrier_peak_fraction=per_carrier_arrival.tolist(),minimum_combined_fraction=.25,
        label='variance-weighted packet envelope arrival diagnostic, not stationaryvariance or Hm0 proof')
    frequencies,PSD = periodogram(eta,fs=50.,window='boxcar',detrend='constant',scaling='density')
    df = frequencies[1]-frequencies[0]
    edges = carriers['frequency_edges_hz']
    band_variance = [float(PSD[(frequencies >= l)&(frequencies < r)].sum()*df) for l,r in zip(edges[:-1],edges[1:])]
    full = data['full_component_spectral_norm']
    retained = data['retained_component_spectral_norm']
    norm = data['spectral_norm_weights']
    actual_headings = train.theta[:,None]+np.array([0.,np.pi])[None,:]
    direction_norm = norm.sum(axis=2)
    moments = [np.sum(direction_norm*np.exp(1j*m*actual_headings))/direction_norm.sum() for m in [1,2]]
    support = dict(omega_range_rad_s=[float(train.omega.min()),float(train.omega.max())],
        removed_component_spectral_norm=(full-retained).tolist(),
        removed_summed_spectral_norm_fraction=float((full.sum()-retained.sum())/full.sum()),
        reverse_retained_norm_fraction=float(data['represented_reverse_component_norm'].sum()/retained.sum()),
        norm_interpretation='componentwise directional train Parseval norm per crestlength; NOT sum global2D energy or stationaryvariance',
        actual_directional_norm_moments=[[float(v.real),float(v.imag)] for v in moments],
        frequency_pruning='No lowcut; Omega>lambda>=2beam cap physical; no renormalization',
        H0_constant_pressure=H0.tolist(),H0_Kheave_max_error=H0_error,
        low_frequency_bridge='H0+(Omega/.1)^2*(H(.1,theta)-H0); approximation, no radiation/Ainf/H5 change')
    spectral = dict(source_observed_band_variance_m2=(np.asarray(observed['band_amplitudes_m'])**2/2).tolist(),
        source_observed_Hm0_m=float(4*np.sqrt(np.sum(np.asarray(observed['band_amplitudes_m'])**2/2))),
        actual60s_centered_variance_m2=float(np.var(eta)),actual60s_mean_eta_m=float(eta.mean()),
        actual60s_periodogram_band_variance_m2=band_variance,
        actual60s_periodogram_total_variance_m2=float(PSD.sum()*df),
        label='actual finite60s hullpoint periodogram, not unchanged stationary observed sea; no energy/phase/amplitude fit')
    initial = dict(initial_force_N_Nm=force[0].tolist(),force_limits=limit.tolist(),
        carrier_force_scale=scale.tolist(),degenerate_axes=(scale == 0).tolist(),force_pass=force_pass)
    savemat(args.output/'body-force.mat',dict(time=times[:,None],force=force))
    np.savez(args.output/'spectral-force.npz',omega=train.omega,spectral_force=spectral_force)
    np.savez(args.output/'finite-field-PSD.npz',time=times,elevation_meanXY=eta,frequency_hz=frequencies,
        PSD_m2_per_hz=PSD,combined_arrival_envelope=combined_arrival)
    write(args.output/'initial-force-screen.json',dict(initial=initial,support=support,spectral=spectral,arrival=arrival,
        elapsed_s=time.monotonic()-started,body_simulated=False))
    gates = dict(initial_force=force_pass,arrival=arrival_pass)
    bounds = None
    if force_pass:
        hull = physical['minimum_hull']; q = hull['equilibrium_quaternion_wxyz']
        R = Rotation.from_quat([*q[1:],q[0]]).as_matrix()
        v = (np.asarray(hull['vertices_material_m'])-hull['cg_material_m'])@R.T+physical['world_cg_m']
        wet = np.asarray(physical['minimum_wetted_triangles_world_m']).reshape(-1,3)
        radius = np.linalg.norm(v[:,:2],axis=1).max()
        bounds = train.initial_bounds(radius,wet[:,2].min(),.05,.05,started+480)
        gates['eta'] = bool(bounds['eta_region_bound_m'] <= .001*carriers['amplitudes'].sum())
        gates['velocity'] = bool(bounds['velocity_region_bound_m_s'] <= .001*np.sum(carriers['amplitudes']*carriers['angular_frequencies']))
    report = dict(complete=True,accepted=bool(all(gates.values())),gates=gates,initial=initial,
        bounds=bounds,initial_eta_limit_m=float(.001*carriers['amplitudes'].sum()),
        initial_speed_limit_m_s=float(.001*np.sum(carriers['amplitudes']*carriers['angular_frequencies'])),
        support=support,spectral=spectral,arrival=arrival,body_simulated=False,elapsed_s=time.monotonic()-started,
        input_sha256={p.as_posix():sha(p) for p in [controls_path,carriers_path,case_path,observation_path,physical_path,table_path]},
        scope='Fixed flat-top spatialtrain, then exact signed-k dispersion; no temporalwave/force ramp or bodyhistoryseed')
    assert time.monotonic()-started <= controls['maximum_elapsed_seconds']
    write(args.output/'diagnostic.json',report)
    print(json.dumps(report),flush=True)
    return 0 if report['accepted'] else 2


if __name__ == '__main__':
    with threadpool_limits(1):
        raise SystemExit(main())
