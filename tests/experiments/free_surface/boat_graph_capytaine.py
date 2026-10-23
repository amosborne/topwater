"""Independent A-infinity control for the tilted snapshot's infinite tangent plane.

The graph BVP uses a tapered surface, so this is a nearby physical control,
not an identical nonflat boundary problem or time-domain validation.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np

from topwater import SimulationInputs
from topwater.constants import WATER_DENSITY, GRAVITY
from tests.benchmark_rankine import peak_rss_bytes
from tests.capytaine_reference import radiation_columns, MODES
from .boat_hos_continuation import tilted_geometry
from .boat_spatial_work import subdivide


def run(work_path, output):
    work = json.loads(work_path.read_text(encoding='utf-8'))
    if not work['complete']:
        raise ValueError('complete parent BVP report required')
    for name, digest in work['producer_sha256'].items():
        if sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'parent producer changed: {name}')
    inputs_path = work_path.with_suffix('.source')/'inputs.json'
    inputs = SimulationInputs.from_json(inputs_path.read_text(encoding='utf-8'))
    if sha256(inputs.to_json().encode()).hexdigest() != work['input_sha256']:
        raise ValueError('parent physical input identity changed')
    triangles, level, slope = tilted_geometry(inputs)
    for _ in range(work['body_refinement']-1):
        triangles = subdivide(triangles)
    if sha256(triangles.tobytes()).hexdigest() != work['panel_sha256']:
        raise ValueError('parent wetted panels differ')
    archive = output.with_suffix('.source')
    shutil.copytree(work_path.with_suffix('.source'), archive)
    hashes = dict(work['producer_sha256'])
    for source in [Path(__file__), Path('tests/capytaine_reference.py')]:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(work_path, archive/'parent-report.json')
    report = dict(complete=False, qualified=False, scope=__doc__, producer_sha256=hashes,
        input_sha256=work['input_sha256'], panel_sha256=work['panel_sha256'],
        parent_report_sha256=sha256(work_path.read_bytes()).hexdigest())
    started = perf_counter()
    try:
        # Independent coordinate construction; all six DOFs must rotate together.
        n = np.r_[-slope, 1.]
        n /= np.linalg.norm(n)
        t = np.array([1., 0., slope[0]])
        t /= np.linalg.norm(t)
        frame = np.vstack((t, np.cross(n, t), n))
        origin = np.r_[inputs.initial_state.position[:2], level]
        rotated = (triangles-origin)@frame.T
        bound = 64*np.finfo(float).eps*max(1., float(np.max(abs(triangles-origin))))
        maximum = float(rotated[..., 2].max())
        if maximum > bound:
            raise ValueError('resolved above-plane reference vertex')
        rotated[..., 2] = np.minimum(rotated[..., 2], 0.)
        center = frame@(inputs.initial_state.position-origin)
        result = radiation_columns(rotated, center, inputs.body.mass, np.inf, MODES,
                                   density=WATER_DENSITY, gravity=GRAVITY)
        transform = np.zeros((6, 6))
        transform[:3, :3] = transform[3:, 3:] = frame
        world = transform.T@result['added_mass']@transform
        comparisons = []
        for row, column in zip(work['rows'][2:], [2, 3], strict=True):
            graph = np.asarray(row['instantaneous_added_mass_column'])
            plane = np.asarray(row['tangent_plane_added_mass_column'])
            reference = world[:, column]
            comparisons.append(dict(mode=row['name'], capytaine_added_mass_column=reference.tolist(),
                graph_minus_capytaine=(graph-reference).tolist(),
                tangent_plane_minus_capytaine=(plane-reference).tolist(),
                graph_diagonal_relative_error=float(graph[column]/reference[column]-1),
                tangent_plane_diagonal_relative_error=float(plane[column]/reference[column]-1)))
        report.update(complete=True, reference_angular_frequency='infinite',
            density=WATER_DENSITY, gravity=GRAVITY, water_depth='infinite',
            reference_world_added_mass=world.tolist(), comparisons=comparisons,
            maximum_vertex_height_before_roundtrip_repair_m=maximum,
            roundtrip_roundoff_bound_m=bound, green_function_settings=result['green_function_settings'])
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        report['peak_process_rss_bytes'] = peak_rss_bytes()
        output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-report', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists():
        parser.error('fresh report and source archive required')
    run(args.work_report, args.output)
