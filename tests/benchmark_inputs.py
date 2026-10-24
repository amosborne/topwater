"""Time initial input validation, not a timestep or hydrodynamic solve.

Run: uv run python -m tests.benchmark_inputs
The same canonical hull is refined, with generation excluded from validation
timings. Results are JSON on stdout; redirect to artifacts/ for local records.
"""

import json
import platform
import statistics
from time import perf_counter

import numpy as np
import scipy

from topwater import HullMesh
from .reference_boat import mesh_arrays


def main():
    results = []
    for refinement in (2, 4, 6):
        start = perf_counter()
        vertices, faces, _, _ = mesh_arrays(refinement)
        generation_seconds = perf_counter() - start
        HullMesh(vertices, faces)  # Warm up imports and numerical libraries.
        times = []
        for _ in range(5):
            start = perf_counter()
            HullMesh(vertices, faces)
            times.append(perf_counter() - start)
        results.append(dict(
            refinement=refinement, panels=len(faces), vertices=len(vertices),
            generation_seconds=generation_seconds,
            validation_median_seconds=statistics.median(times),
            validation_min_seconds=min(times),
            mesh_bytes=vertices.nbytes + faces.nbytes,
        ))
    print(json.dumps(dict(
        case="canonical 24 ft boat; initial geometry validation only",
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        platform=platform.platform(), processor=platform.processor(),
        results=results,
    ), indent=2))


if __name__ == "__main__":
    main()
