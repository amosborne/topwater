"""Fresh-worker flat verification or canoe auxiliary solve; not timesteps."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from time import perf_counter
import numpy as np
from threadpoolctl import threadpool_limits, threadpool_info
from tests.benchmark_rankine import benchmark_environment, peak_rss_bytes
from topwater import SimulationInputs
from topwater.hydrostatics import equilibrate_heave
from tests.reference_hull import input_data
from .canoe import fluid_domain, solve_inertia
from .verify import run


def canoe_run(data, resolution):
    start = perf_counter()
    rest, wet = equilibrate_heave(data.body, data.initial_state)
    mesh, body, free = fluid_domain(wet, resolution)
    added_mass, affine, diagnostics = solve_inertia(mesh, body, free, rest.position)
    return dict(
        seconds=perf_counter() - start,
        peak_rss_bytes=peak_rss_bytes(),
        faces=len(mesh.faces),
        vertices=len(mesh.vertices),
        wetted_faces=len(body),
        added_mass=added_mass.tolist(),
        affine=affine,
        diagnostics=diagnostics,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--canoe",
        action="store_true",
        help="equilibrium through ten-RHS coupling solve",
    )
    parser.add_argument(
        "--quadratic",
        action="store_true",
        help="matched P1/P2 five-period flat-wave comparisons",
    )
    parser.add_argument(
        "--fine",
        action="store_true",
        help="include expensive fine pair for --quadratic",
    )
    parser.add_argument(
        "--degree", type=int, choices=(1, 2), default=2, help=argparse.SUPPRESS
    )
    parser.add_argument("--worker", nargs=3, type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.canoe and args.quadratic:
        parser.error("choose either --canoe or --quadratic")
    case = "quadratic" if args.quadratic else ("canoe" if args.canoe else "flat")
    if args.worker:
        resolution = args.worker
        data = (
            SimulationInputs.from_json(json.dumps(input_data(4)))
            if args.canoe
            else None
        )

        def run_case():
            if args.quadratic:
                from .quadratic import wave

                row = wave(resolution, degree=args.degree)
                row["flat_verification_seconds"] = row["seconds"]
                row["seconds"] = row.pop("total_seconds")
                row["peak_rss_bytes"] = peak_rss_bytes()
                return row
            return canoe_run(data, resolution) if args.canoe else run(resolution)

        with threadpool_limits(1):
            rows = [run_case() for repeat in range(6)]
            threadpools = threadpool_info()
        times = [row["seconds"] for row in rows[1:]]
        print(
            json.dumps(
                dict(
                    case=case,
                    degree=args.degree if args.quadratic else None,
                    resolution=resolution,
                    first=rows[0],
                    warm_seconds=times,
                    warm_median_seconds=float(np.median(times)),
                    warm_min_seconds=min(times),
                    peak_rss_bytes=max(row["peak_rss_bytes"] for row in rows),
                    environment=benchmark_environment(),
                    blas_threads=1,
                    threadpools=threadpools,
                )
            )
        )
    else:
        rows = []
        cases = [((8, 4, 4), 2), ((16, 8, 8), 2)]
        if args.quadratic:
            cases = [((16, 8, 8), 1), ((8, 4, 4), 2)]
            if args.fine:
                cases.extend([((32, 16, 16), 1), ((16, 8, 8), 2)])
        for resolution, degree in cases:
            with tempfile.TemporaryDirectory(prefix="topwater-galerkin-") as cache:
                env = os.environ.copy()
                env["NUMBA_CACHE_DIR"] = cache
                command = [
                    sys.executable,
                    "-m",
                    "tests.experiments.free_surface.benchmark",
                    "--worker",
                    *map(str, resolution),
                ]
                if args.canoe:
                    command.append("--canoe")
                if args.quadratic:
                    command.extend(["--quadratic", "--degree", str(degree)])
                process = subprocess.run(
                    command, env=env, capture_output=True, text=True, check=True
                )
            row = json.loads(process.stdout)
            rows.append(row)
            print(json.dumps(row), flush=True)
        Path("artifacts").mkdir(exist_ok=True)
        suffix = "-quadratic" if args.quadratic else ("-canoe" if args.canoe else "")
        Path(f"artifacts/free-surface-benchmark{suffix}.json").write_text(
            json.dumps(rows, indent=2) + "\n",
            encoding="utf-8",
        )
