"""Supplement the same frozen packet's calm check over the full hull planform."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits
from incident_packet import IncidentPacket


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    packet = Path('artifacts/reference/minimum-20261024/packet-attempt3')
    prior = json.loads((packet/'diagnostic.json').read_text())
    assert prior['accepted']
    controls_path = Path('reference/wecsim/packet-controls.json')
    controls = json.loads(controls_path.read_text())
    assert hashlib.sha256(controls_path.read_bytes()).hexdigest() == prior['input_sha256'][str(controls_path)]
    case = json.loads(Path('reference/capytaine/canonical-boat.json').read_text())
    field = IncidentPacket(case, controls)
    old = np.load(packet/'packet.npz')
    assert np.array_equal(field.c, old['c']) and np.array_equal(field.w, old['w'])
    hull = case['minimum_hull']
    q = hull['equilibrium_quaternion_wxyz']
    R = Rotation.from_quat([*q[1:], q[0]]).as_matrix()
    v = (np.asarray(hull['vertices_material_m'])-hull['cg_material_m'])@R.T+case['world_cg_m']
    wet = np.asarray(case['minimum_wetted_triangles_world_m']).reshape(-1, 3)
    lo = np.r_[v[:, :2].min(axis=0), wet[:, 2].min()]
    hi = np.r_[v[:, :2].max(axis=0), 0.]
    h = controls['initial_grid_spacing_m']
    axes = [np.linspace(lo[d], hi[d], int(np.ceil((hi[d]-lo[d])/h))+1) for d in range(3)]
    xyz = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)
    bounds = field.initial_bounds(xyz, h, started+480)
    bounds['region_min_xyz_m'] = lo.tolist()
    bounds['region_max_xyz_m'] = hi.tolist()
    gates = dict(eta=bool(bounds['eta_region_bound_m'] <= prior['initial_eta_limit_m']),
                 velocity=bool(bounds['velocity_region_bound_m_s'] <= prior['initial_speed_limit_m_s']))
    report = dict(complete=True, accepted=all(gates.values()), gates=gates, initial=bounds,
                  initial_eta_limit_m=prior['initial_eta_limit_m'], initial_speed_limit_m_s=prior['initial_speed_limit_m_s'],
                  prior_force_arrival_band_report=str(packet/'diagnostic.json'),
                  frozen_packet_sha256=hashlib.sha256((packet/'packet.npz').read_bytes()).hexdigest(),
                  unchanged_spectral_coefficients=True, no_body_simulated=True, elapsed_s=time.monotonic()-started)
    (args.output/'diagnostic.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
    return 0 if report['accepted'] else 2


if __name__ == '__main__':
    with threadpool_limits(1):
        raise SystemExit(main())
