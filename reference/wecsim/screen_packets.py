"""Fixed-order incident-only screening; never select using boat motion."""
import argparse
import copy
import json
from pathlib import Path
import time

import numpy as np
from scipy.io import loadmat
from threadpoolctl import threadpool_limits
import xarray as xr
from incident_packet import IncidentPacket, excitation


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    screen = json.loads(Path('reference/wecsim/packet-screen-controls.json').read_text())
    base = json.loads(Path(screen['base_controls']).read_text())
    case = json.loads(Path('reference/capytaine/canonical-boat.json').read_text())
    ds = xr.open_dataset('artifacts/reference/minimum-20261024/native-gl2/canonical.nc')
    H = excitation(ds)
    omega, directions = ds.omega.values, ds.wave_direction.values
    b = loadmat('artifacts/reference/minimum-20261024/consumer/bemio.mat', simplify_cells=True)
    tau = b['h']['ex_t']
    K = b['h']['ex_K']*case['water']['density_kg_m3']*case['water']['gravity_m_s2']
    to = np.arange(round(10/base['dt_s'])+1)*base['dt_s']
    rows, selected = [], None
    for tc in screen['focus_time_candidates_s_in_order']:
        if time.monotonic()-started > screen['maximum_elapsed_seconds']:
            raise TimeoutError('Declared fixed-order packet screen runtime exceeded')
        controls = copy.deepcopy(base)
        controls['focus_time_s'][screen['leading_component_index']] = tc
        field = IncidentPacket(case, controls)
        per_direct = np.zeros((3, 6))
        per_retained = np.zeros((3, 6))
        scale = np.zeros(6)
        direct_history = np.zeros((len(to), 6))
        for j in range(3):
            idir = int(np.argmin(abs(directions-field.theta[j])))
            eta_history = (np.exp(1j*np.outer(tau, field.w[j]))@field.c[j]).real
            for q in range(6):
                hj = np.interp(field.w[j], omega, H[q, idir].real)+1j*np.interp(field.w[j], omega, H[q, idir].imag)
                per_direct[j, q] = np.sum(field.c[j]*hj).real
                per_retained[j, q] = np.sum(eta_history*K[q, j])*base['dt_s']
                direct_history[:, q] += (np.exp(-1j*np.outer(to, field.w[j]))@(field.c[j]*hj)).real
                hc = np.interp(field.w0[j], omega, H[q, idir].real)+1j*np.interp(field.w0[j], omega, H[q, idir].imag)
                scale[q] += field.a[j]*abs(hc)
        origin = np.column_stack([(np.exp(-1j*np.outer(to, field.w[j]))@field.c[j]).real for j in range(3)])
        arrival = np.max(abs(origin), axis=0)/field.a
        limits = np.minimum(base['maximum_initial_force_N_or_Nm'], base['calm_fraction']*scale)
        direct, retained = per_direct.sum(axis=0), per_retained.sum(axis=0)
        gates = dict(force_direct=bool(np.all(abs(direct) <= limits)),
                     force_retained=bool(np.all(abs(retained) <= limits)),
                     arrival=bool(np.any(arrival >= base['minimum_arrival_fraction'])))
        passed = bool(all(gates.values()))
        rows.append(dict(focus_time_s=tc, gates=gates, passed=passed,
                         direct_per_component=per_direct.tolist(), retained_per_component=per_retained.tolist(),
                         direct_total=direct.tolist(), retained_total=retained.tolist(), limits=limits.tolist(),
                         degenerate_axes=(scale < 1e-12).tolist(),
                         initial_direct_retained_difference=(direct-retained).tolist(),
                         arrival_fraction_by_component=arrival.tolist()))
        if selected is None and passed:
            selected = tc
            controls['acceptance_status'] = 'force/arrival screen passed; full hull-region diagnostic pending'
            (args.output/'selected-controls.json').write_text(json.dumps(controls, indent=2)+'\n')
        # Report all frozen candidates even after selecting the first passing one.
    result = dict(complete=True, selected_first_passing_focus_time_s=selected, candidates=rows,
                  no_body_simulated=True, full_hull_region_check_pending=selected is not None,
                  elapsed_s=time.monotonic()-started)
    (args.output/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    return 0 if selected is not None else 2


if __name__ == '__main__':
    with threadpool_limits(1):
        raise SystemExit(main())
