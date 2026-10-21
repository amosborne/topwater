"""Fresh-worker complete flat verification path, not timestep performance."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
from threadpoolctl import threadpool_limits, threadpool_info
from tests.benchmark_rankine import benchmark_environment
from .verify import run


if __name__ == '__main__':
    if len(sys.argv) > 1:
        resolution = json.loads(sys.argv[1])
        with threadpool_limits(1):
            rows = [run(resolution) for repeat in range(6)]
            threadpools = threadpool_info()
        times = [row['seconds'] for row in rows[1:]]
        print(json.dumps(dict(
            resolution=resolution, first=rows[0], warm_seconds=times,
            warm_median_seconds=float(np.median(times)),
            warm_min_seconds=min(times),
            peak_rss_bytes=max(row['peak_rss_bytes'] for row in rows),
            environment=benchmark_environment(), blas_threads=1,
            threadpools=threadpools,
        )))
    else:
        rows = []
        for resolution in ((8, 4, 4), (16, 8, 8)):
            with tempfile.TemporaryDirectory(prefix='topwater-galerkin-') as cache:
                env = os.environ.copy()
                env['NUMBA_CACHE_DIR'] = cache
                process = subprocess.run([
                    sys.executable, '-m', 'tests.experiments.free_surface.benchmark',
                    json.dumps(resolution),
                ], env=env, capture_output=True, text=True, check=True)
            row = json.loads(process.stdout)
            rows.append(row)
            print(json.dumps(row), flush=True)
        Path('artifacts').mkdir(exist_ok=True)
        Path('artifacts/free-surface-benchmark.json').write_text(
            json.dumps(rows, indent=2)+'\n', encoding='utf-8',
        )
