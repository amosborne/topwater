"""Separate canonical hull geometry error from solver-panel error and cost.

python -m tests.experiments.capytaine_mesh_cost
Requires the validation-only correlation group. This is a static mesh study,
not a time-domain wave-body solver or a real-time qualification.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits

from topwater import SimulationInputs
from topwater.hydrostatics import equilibrate_heave
from topwater.inertia import instantaneous_added_mass
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import benchmark_environment, peak_rss_bytes
from tests.capytaine_reference import MODES, frequency_response, instantaneous_inertia
from tests.reference_hull import input_data


FREQUENCIES = np.array([.8, 1.6, 3.2])  # rad/s, fixed before inspecting results
HEADINGS = np.array([0., np.pi/4, np.pi/2])
CASES = ((2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1))


def subdivide(triangles):
    """Four coplanar children; physical surface and waterline do not change."""
    p, q, r = triangles.transpose(1, 0, 2)
    pq, qr, rp = (p+q)/2, (q+r)/2, (r+p)/2
    return np.concatenate([np.stack(corners, axis=1) for corners in (
        (p, pq, rp), (pq, q, qr), (rp, qr, r), (pq, qr, rp))])


def family_error(actual, reference):
    """Contract tolerance over the complete predeclared nondimensional family."""
    floor = max(1e-10, .005*float(np.max(np.abs(reference))))
    error = np.abs(actual-reference)
    ratio = error/(.05*np.abs(reference)+floor)
    valid = (np.abs(reference) > floor) & (np.abs(actual) > floor)
    phase = np.abs(np.angle(actual[valid]*np.conj(reference[valid])))
    return dict(max_absolute_error=float(error.max()), absolute_floor=floor,
                max_tolerance_ratio=float(ratio.max()), within_tolerance=bool(np.all(ratio <= 1)),
                max_phase_error_radians=float(phase.max()) if phase.size else None,
                phase_compared_entries=int(phase.size))


def encode(array):
    return dict(real=array.real.tolist(), imag=array.imag.tolist()) if np.iscomplexobj(array) else array.tolist()


def decode(value):
    return np.array(value['real'])+1j*np.array(value['imag']) if isinstance(value, dict) else np.array(value)


def mode_errors(actual, reference):
    """Expose weak modes that a family-wide absolute floor can conceal."""
    if actual.shape[-2:] == (6, 6):
        a, b = np.diagonal(actual, axis1=-2, axis2=-1), np.diagonal(reference, axis1=-2, axis2=-1)
        error = np.abs(a-b)
        relative = np.full(b.shape, None, dtype=object)
        nonzero = b != 0
        relative[nonzero] = error[nonzero]/np.abs(b[nonzero])
        return dict(diagonal_reference=b.tolist(), diagonal_absolute_error=error.tolist(),
                    diagonal_relative_error=relative.tolist())
    error = np.max(np.abs(actual-reference), axis=(0, 1))
    peak = np.max(np.abs(reference), axis=(0, 1))
    return dict(mode_peak_reference=peak.tolist(), mode_max_absolute_error=error.tolist(),
        mode_relative_to_peak_error=[float(e/r) if r != 0 else None for e, r in zip(error, peak)])


def run_case(geometry, subdivisions):
    start = perf_counter()
    data = SimulationInputs.from_json(json.dumps(input_data(geometry)))
    rest, wet = equilibrate_heave(data.body, data.initial_state)
    triangles = wet.wetted_triangles
    for _ in range(subdivisions):
        triangles = subdivide(triangles)
    panels = TrianglePanels(triangles, -1)
    geometry_seconds = perf_counter()-start
    with threadpool_limits(1):
        times = []
        for _ in range(6):
            start = perf_counter()
            candidate = instantaneous_added_mass(panels, rest.position)
            times.append(perf_counter()-start)
        candidate_peak = peak_rss_bytes()
        start = perf_counter()
        reference, _ = instantaneous_inertia(triangles, rest.position, data.body.mass)
        infinity_seconds = perf_counter()-start
        start = perf_counter()
        response = frequency_response(triangles, rest.position, data.body.mass, FREQUENCIES, HEADINGS)
        frequency_seconds = perf_counter()-start
        pools = threadpool_info()
    scale = np.array([1., 1., 1., 3.5, 3.5, 3.5])
    inertia_scale = 1000*3.5**3*scale[:, None]*scale[None, :]
    normalized = dict(added_mass=response['added_mass']/inertia_scale,
        radiation_damping=response['radiation_damping']/inertia_scale/FREQUENCIES[:, None, None],
        excitation=response['excitation']/(1000*9.81*3.5**2*scale))
    values = [*response.values(), candidate.added_mass, reference]
    if not all(np.all(np.isfinite(a)) for a in values):
        raise RuntimeError('nonfinite mesh-study response')
    return dict(geometry_refinement=geometry, panel_subdivisions=subdivisions,
        original_wetted_faces=len(wet.wetted_triangles), solver_panels=len(panels),
        center_of_mass=rest.position.tolist(), mass=data.body.mass,
        equilibrium_heave_error_m=float(rest.position[2]-data.initial_state.position[2]),
        displaced_volume=wet.displaced_volume, waterplane_area=wet.waterplane_area,
        wetted_area=wet.wetted_area, geometry_seconds=geometry_seconds,
        candidate_first_seconds=times[0], candidate_warm_seconds=times[1:],
        candidate_warm_median_seconds=float(np.median(times[1:])),
        candidate_peak_rss_bytes=candidate_peak, process_peak_rss_bytes=peak_rss_bytes(),
        candidate_instantaneous_added_mass=candidate.added_mass.tolist(),
        capytaine_instantaneous_added_mass=reference.tolist(),
        same_mesh_instantaneous_correlation=(family_error(candidate.added_mass/inertia_scale, reference/inertia_scale)
            | mode_errors(candidate.added_mass/inertia_scale, reference/inertia_scale)),
        capytaine_instantaneous_seconds=infinity_seconds,
        capytaine_frequency_seconds=frequency_seconds,
        capytaine_response={name: encode(value) for name, value in response.items()},
        normalized_families={name: encode(value) for name, value in normalized.items()},
        environment=benchmark_environment(), threadpools=pools,
        candidate_diagnostics=candidate.diagnostics)


def compare(rows):
    by_case = {(r['geometry_refinement'], r['panel_subdivisions']): r for r in rows}
    pairs = [('solver_panels', (g, s), (g, s+1)) for g, s in CASES if (g, s+1) in by_case]
    # These also contain remaining panel error; adjacent agreement is not truth.
    pairs += [('hull_geometry_and_remaining_panel_error', (2, 2), (3, 2)),
              ('hull_geometry_and_remaining_panel_error', (3, 2), (4, 1)),
              ('hull_geometry_and_remaining_panel_error', (2, 1), (3, 0))]
    result = []
    for kind, coarse, fine in pairs:
        if coarse not in by_case or fine not in by_case:
            continue
        a, b = by_case[coarse], by_case[fine]
        metrics = {}
        for name, value in b['normalized_families'].items():
            actual, reference = decode(a['normalized_families'][name]), decode(value)
            metrics[name] = family_error(actual, reference) | mode_errors(actual, reference)
        result.append(dict(varied=kind, coarse=coarse, fine=fine, families=metrics))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', nargs=2, type=int, help=argparse.SUPPRESS)
    parser.add_argument('--case', nargs=2, type=int, action='append', metavar=('GEOMETRY', 'SUBDIVISIONS'),
                        help='run selected predeclared cases; default: first three fixed-geometry refinements')
    args = parser.parse_args()
    destination = Path('artifacts')
    destination.mkdir(exist_ok=True)
    if args.worker:
        row = run_case(*args.worker)
        path = destination/f'capytaine-mesh-g{args.worker[0]}-s{args.worker[1]}.json'
        path.write_text(json.dumps(row, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(dict(completed=args.worker, panels=row['solver_panels'],
            warm_median_seconds=row['candidate_warm_median_seconds'],
            reference_frequency_seconds=row['capytaine_frequency_seconds'])), flush=True)
        return
    rows = []
    report = dict(qualified=False, capytaine='3.0.0', modes=MODES,
        frequencies_rad_per_s=FREQUENCIES.tolist(), headings_radians=HEADINGS.tolist(),
        water_depth='infinite', density=1000., gravity=9.81, wave_amplitude=1.,
        convention='Re(amplitude*exp(-i*omega*t)); wave elevation phase origin (0,0); world loads about CG',
        scope='Static resolution/cost evidence. Topwater comparison is infinite-frequency only; finite-frequency matrices and total excitation are Capytaine references, not Topwater validation.',
        cases=rows, comparisons=[])
    cases = [tuple(case) for case in args.case] if args.case else CASES[:3]
    if any(case not in CASES for case in cases) or len(set(cases)) != len(cases):
        parser.error(f'choose unique cases from {CASES}')
    for geometry, subdivisions in cases:
        with tempfile.TemporaryDirectory(prefix='topwater-mesh-cost-') as cache:
            env = os.environ.copy()
            env['NUMBA_CACHE_DIR'] = cache
            command = [sys.executable, '-m', 'tests.experiments.capytaine_mesh_cost',
                       '--worker', str(geometry), str(subdivisions)]
            process = subprocess.run(command, env=env, capture_output=True, text=True)
        if process.returncode:
            raise RuntimeError(process.stdout+'\n'+process.stderr)
        row = json.loads((destination/f'capytaine-mesh-g{geometry}-s{subdivisions}.json').read_text(encoding='utf-8'))
        rows.append(row)
        report['comparisons'] = compare(rows)
        (destination/'capytaine-mesh-cost.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        print(process.stdout, end='', flush=True)


if __name__ == '__main__':
    main()
