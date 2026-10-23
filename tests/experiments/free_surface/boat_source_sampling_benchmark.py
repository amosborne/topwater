"""Warm paired CPU/GPU sampling and full-RHS controls at an owned evolved state."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from argparse import Namespace
from time import perf_counter
from unittest.mock import patch
import numpy as np

from topwater import SimulationInputs
from tests.benchmark_rankine import peak_rss_bytes
from tests.forced_graph_case import boundary_data
from .boat_forced_motion import ForcedBody, state_from_vector, run as run_motion
from .boat_benchmark_contract import performance_comparison
from .boat_motion_metrics import body_poses, motion_errors, attitude_degrees
from .boat_integration_controls import RELATIVE_TOLERANCE
from .boat_graph_geometry import constant_twist_pose, wet_snapshot
from .boat_image_plane import elevated_potential as serial_potential
from .cuda_source_potential import paired_potential, SourcePotential
from .cuda_graph import cp


def cpu_reference(points, panels, height):
    return SourcePotential(cp.asarray(serial_potential(points, panels, height)),
        dict(implementation='independent guarded CPU reference'))


def run(inputs, parent, boundary_path, report, args):
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
    warmup_started=perf_counter()
    reference = [serial_potential(points, p, model.image) for p in panels]
    candidate = [paired_potential(points, p, model.image) for p in panels]
    report['source_warmup_seconds']=perf_counter()-warmup_started
    errors=[]
    for a,b in zip(reference,candidate,strict=True):
        np.testing.assert_allclose(cp.asnumpy(b.values), a, atol=3e-14, rtol=3e-12)
        errors.append(float(np.max(abs(cp.asnumpy(b.values)-a))))
    report['source_matrix_max_errors']=errors
    report['gpu_guard_diagnostics']=[b.diagnostics for b in candidate]
    report['source_rows'] = []
    for name in ['cpu', 'gpu', 'gpu', 'cpu']:
        function = cpu_reference if name == 'cpu' else paired_potential
        started = perf_counter()
        values = [function(points, p, model.image) for p in panels]
        for a,b in zip(reference,values,strict=True):
            np.testing.assert_allclose(cp.asnumpy(b.values), a, atol=3e-14, rtol=3e-12)
        report['source_rows'].append(dict(method=name, seconds=perf_counter()-started))

    def evaluate(function):
        # Verification-only substitution; neither moving driver nor source owner
        # has a runtime backend/method switch.
        with patch('tests.experiments.free_surface.boat_graph_bvp.paired_potential', function), \
             patch('tests.experiments.free_surface.boat_frozen_field.paired_potential', function):
            return model.evaluate(time, vector)

    warmup_started=perf_counter()
    rhs_reference = evaluate(cpu_reference)
    evaluate(paired_potential)
    report['rhs_warmup_seconds']=perf_counter()-warmup_started
    report['rhs_rows'] = []
    for name in ['cpu', 'gpu', 'gpu', 'cpu']:
        result = evaluate(cpu_reference if name == 'cpu' else paired_potential)
        pressure_error = float(np.max(abs(result.pressure-rhs_reference.pressure)))
        acceleration_error = float(np.max(abs(result.acceleration-rhs_reference.acceleration)))
        assert pressure_error < 1e-6 and acceleration_error < 1e-9
        report['rhs_rows'].append(dict(method=name, diagnostics=result.diagnostics,
            pressure_max_difference_pa=pressure_error, acceleration_max_difference=acceleration_error))
    report.update(complete=True, source_matrix_parity_within_declared_tolerance=True,
                  scope='Warm same-process ABBA at evolved body pose; includes all current RHS geometry/field/coupling work. Not an evolving-sea or trajectory qualification.')
    if args.moving_case:
        assert controls['rtol'] == RELATIVE_TOLERANCE and controls['absolute_tolerance_scale'] == 1.
        report['complete'] = False
        reports, paths = [], []
        for number, name in enumerate(['cpu', 'gpu', 'gpu', 'cpu']):
            path = args.output.with_name(args.output.stem+f'-moving-{number}-{name}.json')
            assert not path.exists() and not path.with_name(path.stem+'-states.npz').exists()
            settings = Namespace(**{**controls, 'output':path, 'resource_status':args.resource_status,
                'resource_context':args.resource_context})
            case = dict(complete=False, integration_complete=False, qualified=False, rhs=[],
                input_sha256=report['input_sha256'], producer_sha256=report['producer_sha256'],
                implementation=name, controls={k:str(v) if isinstance(v,Path) else v for k,v in vars(settings).items()},
                producer_archive=str(args.output.with_suffix('.source')))
            function = cpu_reference if name == 'cpu' else paired_potential
            try:
                with patch('tests.experiments.free_surface.boat_graph_bvp.paired_potential', function), \
                     patch('tests.experiments.free_surface.boat_frozen_field.paired_potential', function):
                    run_motion(inputs, settings, case)
            finally:
                path.write_text(json.dumps(case,indent=2,allow_nan=False),encoding='utf-8')
            # Both exact current implementations were warmed above on the owned
            # evolved state, before timing these fresh-model integration cases.
            case['measurement_environment'].update(cache_policy='prewarmed current producers',
                startup_treatment='Same process; both current implementations warmed with complete owned-state RHS before timed ABBA cases')
            assert case['heldout']['full_context_oracle']['pressure_max_difference_pa'] < 1e-5
            path.write_text(json.dumps(case,indent=2,allow_nan=False),encoding='utf-8')
            reports.append(case);paths.append(path)
        times = np.linspace(0., controls['duration'], 101)
        comparisons=[]
        for reference_index,candidate_index in [(0,1),(3,2)]:
            comparison=performance_comparison(reports[reference_index],reports[candidate_index])
            with np.load(paths[reference_index].with_name(paths[reference_index].stem+'-states.npz')) as z:
                reference_pose=body_poses(z['times'],z['states'],times,dense_q=z['dense_q'])
                heading=attitude_degrees(z['states'][0,3:7])[2]
            with np.load(paths[candidate_index].with_name(paths[candidate_index].stem+'-states.npz')) as z:
                pose=body_poses(z['times'],z['states'],times,dense_q=z['dense_q'])
            translation,rotation=motion_errors(reference_pose,pose,heading)
            assert np.max(abs(translation))<1e-10 and np.max(abs(rotation))<1e-8
            comparisons.append(dict(reference=str(paths[reference_index]),candidate=str(paths[candidate_index]),
                performance=comparison,max_translation_difference_m=np.max(abs(translation),axis=0).tolist(),
                max_rotation_difference_degrees=np.max(abs(rotation),axis=0).tolist()))
        report.update(moving_comparisons=comparisons,moving_reports=[str(path) for path in paths],complete=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--parent-motion', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resource-context', default='Uncharacterized shared-machine load')
    parser.add_argument('--resource-status', choices=['uncontrolled','contended','uncontended'],default='uncontrolled')
    parser.add_argument('--moving-case', action='store_true',help='also run matched fresh-model moving-body ABBA cases after warmup')
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists():
        parser.error('fresh output and archive required')
    archive_started=perf_counter()
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
        archive_seconds=perf_counter()-archive_started,
        resource_context=args.resource_context,
        resource_status=args.resource_status,
        parent_motion_sha256=sha256(args.parent_motion.read_bytes()).hexdigest(),
        input_sha256=sha256(inputs.to_json().encode()).hexdigest())
    started = perf_counter()
    try:
        run(inputs, parent, boundary, report, args)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['numerical_seconds'] = perf_counter()-started
        report['total_seconds'] = report['archive_seconds']+report['numerical_seconds']
        report['benchmark_process_peak_rss_bytes']=peak_rss_bytes()
        report['ending_cuda_pool_reserved_bytes']=cp.get_default_memory_pool().total_bytes()
        report['memory_scope']='Whole verification process includes CPU/GPU warmups, reference matrices and ABBA; not per-case RAM or peak VRAM.'
        report['timing_scope']='Aggregate includes archive, warmups and verification. Matched moving comparisons separately report warmed model setup/integration only.'
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
