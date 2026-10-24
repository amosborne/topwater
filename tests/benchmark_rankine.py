"""Fresh-process benchmarks: python -m tests.benchmark_rankine --help.

Run on an idle machine. Timings include allocation and assembly, exclude mesh
preparation/imports, and retain no operator from a preceding repetition. Process
peak RSS includes imports, geometry, operators and scratch, not just matrix bytes.
"""

import argparse
import ctypes
import gc
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import tempfile
from time import perf_counter

import numpy as np
import scipy
from threadpoolctl import threadpool_info, threadpool_limits

from topwater import HullMesh
from topwater.bem import assemble_laplace
from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels, panel_influence
from .analytic_geometry import sphere_mesh


def benchmark_environment():
    return dict(format_version=1, python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                numba=version("numba"), platform=platform.platform(), processor=platform.processor(),
                logical_cpu_count=os.cpu_count())


def peak_rss_bytes():
    if sys.platform == "win32":
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
                (name, ctypes.c_size_t) for name in (
                    "PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage",
                    "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage",
                    "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]

        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        counters = Counters()
        counters.cb = ctypes.sizeof(counters)
        if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
            raise ctypes.WinError(ctypes.get_last_error())
        return int(counters.PeakWorkingSetSize)
    import resource
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else 1024*value)


def worker(args):
    vertices, faces = sphere_mesh(args.refinement)
    panels = TrianglePanels.from_mesh(HullMesh(vertices, faces))
    controls = LaplaceControls(target_block_size=args.block)
    # Every target lies outside the closed hull; K integrates to zero.
    targets = np.column_stack((np.linspace(-2, 2, 128), np.full(128, 1.), np.full(128, .5)))
    baseline = peak_rss_bytes()
    samples = []
    error = 0.0
    fallbacks = 0
    with threadpool_limits(limits=args.threads):
        libraries = threadpool_info()
        for repeat in range(args.repeats+1):
            gc.collect()
            start = perf_counter()
            if args.operation == "assembly":
                result = assemble_laplace(panels, controls=controls)
            else:
                # Fixed independent bulk targets; geometry/targets prepared outside timing.
                fields = [panel_influence(targets[first:first+args.block], panels,
                                           controls=controls, gradients=True)
                          for first in range(0, len(targets), args.block)]
            elapsed = perf_counter()-start
            if repeat:
                samples.append(elapsed)
            else:
                first_call_seconds = elapsed
            if args.operation == "assembly":
                error = max(error, float(np.max(np.abs(result.double_layer.sum(axis=1)+0.5))))
                fallbacks += result.diagnostics["fallback_pairs"]
            else:
                for field in fields:
                    error = max(error, float(np.max(np.abs(field.double_layer.sum(axis=1)))))
                    fallbacks += field.fallback_pairs
                del field, fields
            if args.operation == "assembly":
                del result
    return dict(operation=args.operation, panel_count=len(panels), block_rows=args.block,
                kernel_threads=1,
                repeats=args.repeats, warmups=1, seconds=samples,
                first_call_seconds=first_call_seconds, first_call_includes_compilation=True,
                target_count=len(panels) if args.operation == "assembly" else len(targets),
                pair_count=len(panels)*(len(panels) if args.operation == "assembly" else len(targets)),
                median_seconds=statistics.median(samples), minimum_seconds=min(samples),
                process_peak_rss_bytes=peak_rss_bytes(), baseline_peak_rss_bytes=baseline,
                operator_bytes=16*len(panels)**2 if args.operation == "assembly" else 0,
                closed_identity_error=error, fallback_pairs_including_warmup=fallbacks,
                libraries=libraries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--refinement", type=int, default=4, help=argparse.SUPPRESS)
    parser.add_argument("--block", type=int, default=128, help=argparse.SUPPRESS)
    parser.add_argument("--operation", choices=["assembly", "field"], default="assembly")
    parser.add_argument("--refinements", type=int, nargs="+", default=[2, 3, 4])
    parser.add_argument("--blocks", type=int, nargs="+", default=[16, 32, 64, 128])
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", type=Path, default=Path("artifacts/rankine-benchmark.json"))
    args = parser.parse_args()
    if args.repeats < 5 or args.threads < 1:
        parser.error("use at least five measured repeats and one thread")
    if args.worker:
        print(json.dumps(worker(args), allow_nan=False))
        return
    results = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for refinement in args.refinements:
        for block in args.blocks:
            command = [sys.executable, "-m", "tests.benchmark_rankine", "--worker",
                       "--refinement", str(refinement), "--block", str(block),
                       "--operation", args.operation, "--repeats", str(args.repeats),
                       "--threads", str(args.threads)]
            # Isolate compilation cache, not just the process. First-call cost is
            # reproducible cold compilation rather than an accidental disk hit.
            with tempfile.TemporaryDirectory(dir=args.output.parent) as cache:
                env = dict(os.environ, NUMBA_CACHE_DIR=str(Path(cache).resolve()))
                completed = subprocess.run(command, check=True, capture_output=True, text=True, env=env)
            result = json.loads(completed.stdout)
            results.append(result)
            print(f"N={result['panel_count']:4} block={block:3}: "
                  f"{result['median_seconds']:.4f} s, peak RSS "
                  f"{result['process_peak_rss_bytes']/2**20:.1f} MiB", flush=True)
    report = dict(benchmark_environment(), results=results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
