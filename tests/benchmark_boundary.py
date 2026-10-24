"""Fresh-process end-to-end boundary benchmarks, including reusable LU costs."""
import argparse
import json
from pathlib import Path
import statistics
import subprocess
import sys
from time import perf_counter
import os
import tempfile

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits

from topwater import HullMesh
from topwater.bem import assemble_laplace
from topwater.boundary import BoundarySystem
from topwater.rankine import TrianglePanels
from .benchmark_rankine import benchmark_environment, peak_rss_bytes
from .analytic_geometry import sphere_mesh


def worker(args):
    v, f = sphere_mesh(args.refinement)
    p = TrianglePanels.from_mesh(HullMesh(v, f))
    x, y, z = p.centroids.T
    phi = np.column_stack((x, y, z, x*x-y*y, x*y, x*z))
    gradient = np.zeros((len(p), 6, 3))
    gradient[:, 0, 0] = 1
    gradient[:, 1, 1] = 1
    gradient[:, 2, 2] = 1
    gradient[:, 3, 0], gradient[:, 3, 1] = 2*x, -2*y
    gradient[:, 4, 0], gradient[:, 4, 1] = y, x
    gradient[:, 5, 0], gradient[:, 5, 2] = z, x
    q = np.einsum('nrj,nj->nr', gradient, p.normals)
    mask = (np.ones(len(p), bool) if args.kind == 'dirichlet' else
            np.zeros(len(p), bool) if args.kind == 'neumann' else x > 0)
    setup, single, batch = [], [], []
    baseline = peak_rss_bytes()
    with threadpool_limits(None if args.threads == 0 else args.threads):
        libraries = threadpool_info()
        start = perf_counter()
        op = assemble_laplace(p)
        cold_assembly = perf_counter()-start
        for repeat in range(6):
            start = perf_counter()
            system = BoundarySystem(op, dirichlet=mask, areas=p.areas,
                                    length_scale=3.5, bounded=True,
                                    factorization_threads=args.factorization_threads or None)
            elapsed = perf_counter()-start
            if repeat:
                setup.append(elapsed)
            else:
                cold_setup = elapsed
            start = perf_counter()
            solution = system.solve(potential=phi[:, 0], normal_velocity=q[:, 0])
            elapsed = perf_counter()-start
            if repeat:
                single.append(elapsed)
            start = perf_counter()
            solution = system.solve(potential=phi, normal_velocity=q)
            elapsed = perf_counter()-start
            if repeat:
                batch.append(elapsed)
            factor_bytes = system.diagnostics['factor_bytes']
            factor_pools = system.diagnostics['factorization_blas']
            rcond = system.diagnostics['rcond']
            del system
    return dict(panel_count=len(p), kind=args.kind, requested_blas_threads=args.threads,
                requested_factorization_threads=args.factorization_threads or None,
                factorization_blas=factor_pools,
                libraries=libraries, cold_assembly_seconds=cold_assembly,
                cold_setup_seconds=cold_setup,
                setup_seconds=setup, single_rhs_seconds=single, six_rhs_seconds=batch,
                setup_median_seconds=statistics.median(setup),
                single_rhs_median_seconds=statistics.median(single),
                six_rhs_median_seconds=statistics.median(batch),
                operator_bytes=op.single_layer.nbytes+op.double_layer.nbytes,
                factor_bytes=factor_bytes, rcond=rcond,
                backward_error=solution.diagnostics['backward_error'],
                original_bie_residual=solution.diagnostics['original_bie_residual'],
                baseline_peak_rss_bytes=baseline,
                process_peak_rss_bytes=peak_rss_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--refinement', type=int, default=4, help=argparse.SUPPRESS)
    parser.add_argument('--kind', choices=['dirichlet', 'neumann', 'mixed'], default='mixed')
    parser.add_argument('--threads', type=int, default=1, help='0 preserves the host BLAS default')
    parser.add_argument('--factorization-threads', type=int, default=1,
                        help='0 preserves ambient settings; default 1 matches the solver')
    parser.add_argument('--refinements', type=int, nargs='+', default=[2, 3, 4])
    parser.add_argument('--output', type=Path, default=Path('artifacts/boundary-benchmark.json'))
    args = parser.parse_args()
    if args.threads < 0 or args.factorization_threads < 0:
        parser.error('threads must be nonnegative')
    if args.worker:
        print(json.dumps(worker(args), allow_nan=False))
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results = []
    for refinement in args.refinements:
        with tempfile.TemporaryDirectory(dir=args.output.parent) as cache:
            command = [sys.executable, '-m', 'tests.benchmark_boundary', '--worker',
                       '--refinement', str(refinement), '--kind', args.kind,
                       '--threads', str(args.threads),
                       '--factorization-threads', str(args.factorization_threads)]
            completed = subprocess.run(command, check=True, capture_output=True, text=True,
                                       env=dict(os.environ, NUMBA_CACHE_DIR=str(Path(cache).resolve())))
        row = json.loads(completed.stdout)
        results.append(row)
        print(f"N={row['panel_count']}: setup {row['setup_median_seconds']:.4f}s, "
              f"single {row['single_rhs_median_seconds']:.4f}s, "
              f"six {row['six_rhs_median_seconds']:.4f}s, "
              f"peak {row['process_peak_rss_bytes']/2**20:.1f}MiB", flush=True)
    args.output.write_text(json.dumps(dict(benchmark_environment(), results=results), indent=2)
                           +"\n", encoding='utf-8')


if __name__ == '__main__':
    main()
