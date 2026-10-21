"""Fresh-process canoe timings: python -m tests.benchmark_hydrostatics.

Includes transforms, clipping, allocation and pressure integration; equilibrium
includes root iterations and its final pressure evaluation. Imports and input
validation are outside timings. No BEM work is included.
The first hydrostatics call includes JIT compilation; equilibrium follows that
warmup. Neither first-call measurement isolates compilation time alone.
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

from topwater import SimulationInputs
from topwater.hydrostatics import equilibrate_heave, hydrostatics
from .benchmark_rankine import benchmark_environment, peak_rss_bytes
from .reference_hull import input_data


def worker(refinement):
    data = SimulationInputs.from_json(json.dumps(input_data(refinement)))
    report = dict(refinement=refinement, faces=len(data.body.mesh.faces))
    baseline = peak_rss_bytes()
    for function in (hydrostatics, equilibrate_heave):
        samples = []
        for repeat in range(6):
            gc.collect()
            start = perf_counter()
            value = function(data.body, data.initial_state)
            elapsed = perf_counter()-start
            if repeat:
                samples.append(elapsed)
            else:
                first = elapsed
            result = value[1] if isinstance(value, tuple) else value
            if function is equilibrate_heave:
                assert abs(result.force[2]/(data.body.mass*9.81)-1) < 1e-9
            del result, value
        report[function.__name__] = dict(first_seconds=first,
            warm_seconds=samples, median_seconds=statistics.median(samples))
    report.update(baseline_rss_bytes=baseline, peak_rss_bytes=peak_rss_bytes())
    print(json.dumps(report))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', type=int)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.worker is not None:
        worker(args.worker)
    else:
        report = benchmark_environment()
        report['cases'] = []
        for refinement in (2, 3, 4, 5):
            with tempfile.TemporaryDirectory(prefix='topwater-hydro-') as cache:
                environment = dict(os.environ, NUMBA_CACHE_DIR=cache)
                output = subprocess.check_output([sys.executable, '-m', __spec__.name,
                    '--worker', str(refinement)], env=environment, text=True)
                report['cases'].append(json.loads(output))
        serialized = json.dumps(report, indent=2)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(serialized+'\n', encoding='utf-8')
        print(serialized)
