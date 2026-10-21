"""Fresh-process heave equilibrium -> wetted panels -> six-mode fluid inertia.

Imports and validated input construction are excluded. Every timed repetition
recomputes equilibrium, clipping, numerical geometry, image operators, one LU,
all six auxiliary potentials and the added-mass matrix. The first call includes
JIT compilation; five warm calls follow without retaining operators/factors.
"""

import argparse
import gc
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
from time import perf_counter

from threadpoolctl import threadpool_info, threadpool_limits

from topwater import SimulationInputs
from topwater.hydrostatics import equilibrate_heave
from topwater.inertia import instantaneous_added_mass
from topwater.rankine import TrianglePanels
from .benchmark_rankine import benchmark_environment, peak_rss_bytes
from .reference_hull import input_data


def worker(refinement, threads):
    data = SimulationInputs.from_json(json.dumps(input_data(refinement)))
    samples = []
    baseline = peak_rss_bytes()
    with threadpool_limits(threads or None):
        libraries = threadpool_info()
        for repeat in range(6):
            gc.collect()
            start = perf_counter()
            rest, wet = equilibrate_heave(data.body, data.initial_state)
            after_clipping = perf_counter()
            panels = TrianglePanels(wet.wetted_triangles, -1)
            after_geometry = perf_counter()
            result = instantaneous_added_mass(panels, rest.position)
            elapsed = perf_counter()-start
            diagnostics = result.diagnostics
            timing = dict(total=elapsed, equilibrium_and_clipping=after_clipping-start,
                numerical_geometry=after_geometry-after_clipping,
                assembly=diagnostics['assembly']['assembly_seconds'],
                factorization=diagnostics['factorization']['setup_seconds'],
                six_rhs=diagnostics['solve']['solve_seconds'])
            if repeat:
                samples.append(timing)
            else:
                first = timing
            report = dict(refinement=refinement, complete_hull_faces=len(data.body.mesh.faces),
                wetted_panel_count=len(panels), operator_bytes=diagnostics['assembly']['operator_bytes'],
                factor_bytes=diagnostics['factorization']['factor_bytes'],
                rcond=diagnostics['factorization']['rcond'], added_mass=result.added_mass.tolist(),
                normalized_reciprocity_defect=diagnostics['normalized_reciprocity_defect'],
                normalization_length=diagnostics['normalization_length'],
                backward_error=diagnostics['solve']['backward_error'])
            del result, panels, wet
    report.update(first_call_seconds=first, warm_seconds=samples,
        warm_median_seconds={key: statistics.median(sample[key] for sample in samples) for key in first},
        requested_threads=threads or None, libraries=libraries,
        baseline_peak_rss_bytes=baseline, process_peak_rss_bytes=peak_rss_bytes())
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--threads', type=int, default=1, help='0 preserves ambient limits')
    parser.add_argument('--refinements', nargs='+', type=int, default=[3, 4, 5])
    parser.add_argument('--output', type=Path, default=Path('artifacts/inertia-benchmark.json'))
    args = parser.parse_args()
    if args.threads < 0:
        parser.error('threads must be nonnegative')
    if args.worker is not None:
        print(json.dumps(worker(args.worker, args.threads)))
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        report = dict(benchmark_environment(), results=[])
        for refinement in args.refinements:
            with tempfile.TemporaryDirectory(dir=args.output.parent) as cache:
                output = subprocess.check_output([sys.executable, '-m', __spec__.name,
                    '--worker', str(refinement), '--threads', str(args.threads)], text=True,
                    env=dict(os.environ, NUMBA_CACHE_DIR=str(Path(cache).resolve())))
                row = json.loads(output)
            report['results'].append(row)
            print(f"N={row['wetted_panel_count']}: total {row['warm_median_seconds']['total']:.4f}s, "
                f"peak {row['process_peak_rss_bytes']/2**20:.1f}MiB", flush=True)
        args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
