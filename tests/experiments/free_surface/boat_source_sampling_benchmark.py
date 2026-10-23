"""Warm paired serial/parallel sampling and full-RHS controls at an owned evolved state."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
from unittest.mock import patch
import numpy as np

from topwater import SimulationInputs
from topwater.rankine import panel_influence
from tests.forced_graph_case import boundary_data
from .boat_forced_motion import ForcedBody, state_from_vector
from .boat_graph_geometry import constant_twist_pose, wet_snapshot
from .boat_image_plane import elevated_potential


def serial_potential(points, panels, height):
    """Independent serial scheduling reference; shared public accuracy guards."""
    result = np.empty((len(points), len(panels)))
    for first in range(0, len(points), 128):
        last = min(first+128, len(points))
        reflected = points[first:last].copy()
        reflected[:, 2] = 2*height-reflected[:, 2]
        result[first:last] = (panel_influence(points[first:last], panels).single_layer
                              -panel_influence(reflected, panels).single_layer)
    return result


def run(inputs, parent, boundary_path, report):
    assert parent['complete'] and parent['input_sha256'] == sha256(inputs.to_json().encode()).hexdigest()
    assert sha256(boundary_path.read_bytes()).hexdigest() == parent['boundary_sha256']
    controls = parent['controls']
    boundary = boundary_data(inputs, controls['grid'], controls['length'])
    with np.load(boundary_path) as saved:
        assert set(saved.files) == set(boundary)
        assert all(np.array_equal(saved[k], boundary[k]) for k in saved.files)
    model = ForcedBody(inputs, boundary, controls['derivative_step'], controls['body_refinement'], controls['order'])
    assert model.setup == parent['setup']
    vector = np.array(parent['heldout']['state'])
    state, time = state_from_vector(vector), parent['heldout']['time_seconds']
    step = model.step
    panels = [wet_snapshot(inputs.body, pose, boundary['plane_origin'], boundary['frame'],
        boundary['grid_xy'], boundary['length'], model.refinement)[0] for pose in
        [state, constant_twist_pose(state, step), constant_twist_pose(state, -step)]]
    points = np.column_stack((boundary['grid_xy'], boundary['eta'].ravel()))
    reference = [serial_potential(points, p, model.image) for p in panels]
    candidate = [elevated_potential(points, p, model.image) for p in panels]
    assert all(np.array_equal(a, b) for a,b in zip(reference, candidate, strict=True))
    report['source_rows'] = []
    for name in ['serial', 'parallel', 'parallel', 'serial']:
        function = serial_potential if name == 'serial' else elevated_potential
        started = perf_counter()
        values = [function(points, p, model.image) for p in panels]
        assert all(np.array_equal(a, b) for a,b in zip(reference, values, strict=True))
        report['source_rows'].append(dict(method=name, seconds=perf_counter()-started))

    def evaluate(function):
        # Verification-only substitution; neither moving driver nor source owner
        # has a runtime backend/method switch. Worker threads never run CuPy.
        with patch('tests.experiments.free_surface.boat_graph_bvp.elevated_potential', function), \
             patch('tests.experiments.free_surface.boat_frozen_field.elevated_potential', function):
            return model.evaluate(time, vector)

    rhs_reference = evaluate(serial_potential)
    evaluate(elevated_potential)
    report['rhs_rows'] = []
    for name in ['serial', 'parallel', 'parallel', 'serial']:
        result = evaluate(serial_potential if name == 'serial' else elevated_potential)
        pressure_error = float(np.max(abs(result.pressure-rhs_reference.pressure)))
        acceleration_error = float(np.max(abs(result.acceleration-rhs_reference.acceleration)))
        assert pressure_error < 1e-8 and acceleration_error < 1e-11
        report['rhs_rows'].append(dict(method=name, diagnostics=result.diagnostics,
            pressure_max_difference_pa=pressure_error, acceleration_max_difference=acceleration_error))
    report.update(complete=True, exact_source_matrix_parity=True,
                  scope='Warm same-process ABBA at evolved body pose; includes all current RHS geometry/field/coupling work. Not an evolving-sea or trajectory qualification.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--parent-motion', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resource-context', default='Uncharacterized shared-machine load')
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists():
        parser.error('fresh output and archive required')
    archive = args.output.with_suffix('.source'); archive.mkdir(parents=True)
    hashes = {}
    for source in [*Path('src/topwater').glob('*.py'), *Path('tests').glob('*.py'),
                   *Path(__file__).parent.glob('*.py'), *Path(__file__).parent.glob('*.cu'),
                   *Path(__file__).parent.glob('*.cuh')]:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target); hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    shutil.copyfile(args.parent_motion, archive/'parent-motion.json')
    parent = json.loads(args.parent_motion.read_text())
    boundary = args.parent_motion.with_name(args.parent_motion.stem+'-boundary.npz')
    shutil.copyfile(boundary, archive/'boundary.npz')
    inputs = SimulationInputs.from_json(args.inputs.read_text())
    report = dict(complete=False, qualified=False, producer_sha256=hashes,
        resource_context=args.resource_context,
        parent_motion_sha256=sha256(args.parent_motion.read_bytes()).hexdigest(),
        input_sha256=sha256(inputs.to_json().encode()).hexdigest())
    started = perf_counter()
    try:
        run(inputs, parent, boundary, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
