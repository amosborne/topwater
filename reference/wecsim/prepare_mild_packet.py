"""Incident-only preflight; body execution requires separate allocation after PASS."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.io import loadmat, savemat
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits
import xarray as xr
from incident_packet import IncidentPacket, excitation


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    cp = Path('reference/wecsim/mild-packet-controls.json')
    physical = Path('reference/capytaine/canonical-boat.json')
    control = json.loads(cp.read_text())
    case = json.loads(physical.read_text())
    # Keep the reviewed three-slot evaluator; only slot2 has physical energy.
    case['sea']['amplitudes_m'] = [0., control['amplitude_m'], 0.]
    case['sea']['frequencies_hz'][1] = control['carrier_frequency_hz']
    case['sea']['phases_rad'][1] = control['phase_rad']
    case['sea']['propagation_headings_degrees'][1] = control['heading_deg']
    field = IncidentPacket(case, control)
    assert field.a[0] == field.a[2] == 0
    args.output.mkdir(parents=True, exist_ok=False)
    dt = control['dt_s']
    tf = np.arange(-1000, 4001)*dt
    to = np.arange(3001)*dt
    eta = np.array([field.elevation([[0., 0.]], t)[0, 1] for t in tf])
    np.savetxt(args.output/'wave-physical.txt', np.c_[tf, eta], fmt='%.16e')
    np.savez(args.output/'packet.npz', w=field.w, k=field.k, c=field.c,
             weights=field.weights, direction=field.theta, times=tf, elevation=eta,
             carrier_omega=field.w0, carrier_amplitude=field.a, focus_time=field.tc)
    hull = case['minimum_hull']
    q = hull['equilibrium_quaternion_wxyz']
    R = Rotation.from_quat([*q[1:], q[0]]).as_matrix()
    v = (np.asarray(hull['vertices_material_m'])-hull['cg_material_m'])@R.T+case['world_cg_m']
    wet = np.asarray(case['minimum_wetted_triangles_world_m']).reshape(-1, 3)
    lo = np.r_[v[:, :2].min(axis=0), wet[:, 2].min()]
    hi = np.r_[v[:, :2].max(axis=0), 0.]
    spacing = control['initial_grid_spacing_m']
    axes = [np.linspace(lo[d], hi[d], int(np.ceil((hi[d]-lo[d])/spacing))+1) for d in range(3)]
    xyz = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)
    bounds = field.initial_bounds(xyz, spacing, start+control['maximum_elapsed_seconds'])
    bounds.update(region_min_xyz_m=lo.tolist(), region_max_xyz_m=hi.tolist())
    native = Path('artifacts/reference/minimum-20261024/native-gl2/canonical.nc')
    with xr.open_dataset(native) as ds:
        H = excitation(ds)
        omega = ds.omega.values
        directions = ds.wave_direction.values
    idir = int(np.argmin(abs(directions-field.theta[1])))
    assert abs(directions[idir]-field.theta[1]) < 1e-12
    assert field.w[1].min() >= omega.min() and field.w[1].max() <= omega.max()
    direct = np.zeros((len(to), 6))
    scale = np.zeros(6)
    for axis in range(6):
        h = np.interp(field.w[1], omega, H[axis, idir].real)+1j*np.interp(field.w[1], omega, H[axis, idir].imag)
        direct[:, axis] = (np.exp(-1j*to[:, None]*field.w[1])@(field.c[1]*h)).real
        hc = np.interp(field.w0[1], omega, H[axis, idir].real)+1j*np.interp(field.w0[1], omega, H[axis, idir].imag)
        scale[axis] = field.a[1]*abs(hc)
    bemio_path = Path('artifacts/reference/minimum-20261024/consumer/bemio.mat')
    h = loadmat(bemio_path, simplify_cells=True)['h']
    kt = np.asarray(h['ex_t']).reshape(-1)
    assert np.allclose(kt, np.arange(-1000, 1001)*dt, rtol=0, atol=1e-12)
    kernel = h['ex_K']*case['water']['density_kg_m3']*case['water']['gravity_m_s2']
    force = np.stack([np.convolve(eta, kernel[axis, idir], mode='valid')*dt for axis in range(6)], axis=1)
    assert force.shape == direct.shape == (3001, 6)
    limits = np.minimum(.1, .001*scale)
    arrival = np.max(abs(eta[(tf >= 0) & (tf <= 60)]))/field.a[1]
    eta_limit = .001*field.a.sum()
    speed_limit = .001*np.sum(field.a*field.w0)
    fmax = min(.6, np.sqrt(field.g/(2*np.pi*(2*9*.3048))))
    gates = dict(eta=bool(bounds['eta_region_bound_m'] <= eta_limit),
        velocity=bool(bounds['velocity_region_bound_m_s'] <= speed_limit),
        retained_force=bool(np.all(abs(force[0]) <= limits)),
        direct_force=bool(np.all(abs(direct[0]) <= limits)),
        arrival=bool(arrival >= .25), band=bool(field.w[1].max() <= 2*np.pi*fmax))
    assert time.monotonic()-start <= control['maximum_elapsed_seconds']
    report = dict(complete=True, accepted=bool(all(gates.values())), diagnostic_only=True,
        gates=gates, initial=bounds, initial_eta_limit_m=float(eta_limit),
        initial_speed_limit_m_s=float(speed_limit), initial_retained_force=force[0].tolist(),
        initial_direct_force=direct[0].tolist(), force_limits=limits.tolist(),
        degenerate_axes=(scale == 0).tolist(), carrier_force_scale=scale.tolist(),
        arrival_fraction=float(arrival), spectral_range_hz=(field.w[1,[0,-1]]/(2*np.pi)).tolist(),
        direct_retained_max_difference=np.max(abs(force-direct), axis=0).tolist(),
        elapsed_s=time.monotonic()-start,
        input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [cp,physical,native,bemio_path]})
    savemat(args.output/'body-force.mat',dict(time=to[:,None],force=force,direct=direct))
    (args.output/'diagnostic.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report), flush=True)
    return 0 if report['accepted'] else 2


if __name__ == '__main__':
    with threadpool_limits(1):
        raise SystemExit(main())
