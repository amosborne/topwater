"""Compare accepted RK body motion on a common time/heading frame."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .boat_motion_metrics import (attitude_degrees, body_poses, motion_errors,
    ANGLE_ERROR_TARGET_DEGREES, POSITION_ERROR_TARGET_METERS, ACCURACY_DURATION_SECONDS)
from .p4_fluid_ordering import fresh
from .boat_moving import METHODS, SOURCES
from .boat_benchmark_contract import performance_comparison


def run(reference, candidate, output, driver_change_note=None, integration_change_note=None,
        derivative_change_note=None):
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    inputs, trajectories = [], []
    for path in map(Path, (reference, candidate)):
        report = json.loads(path.read_text())
        if (report['format_version'] != 1 or report['methods'] != METHODS
                or set(report['source_sha256']) != set(SOURCES)):
            raise ValueError('current single-path baseline report required')
        physical_path = path.with_name(path.stem+'-inputs.json')
        if sha256(physical_path.read_bytes()).hexdigest() != report['input_sha256']:
            raise ValueError('owned physical inputs changed')
        with np.load(path.with_suffix('.npz'), allow_pickle=False) as data:
            times, states, dense = data['accepted_times'], data['accepted_states'], data['accepted_dense_q']
        if len(times) < 2 or dense.shape[-1] != 3:
            raise ValueError('accepted RK23 trajectory and its owned cubic coefficients required')
        inputs.append(dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest(),
            arrays_sha256=sha256(path.with_suffix('.npz').read_bytes()).hexdigest(),
            input_sha256=report['input_sha256'], complete=report['complete'],
            accepted_seconds=float(times[-1]-times[0]), seconds=report['seconds'],
            operator_counts=report['operator_counts'], discretization=report['discretization']))
        trajectories.append((times, states, report, dense))
    ta, ya, ra, da = trajectories[0]
    tb, yb, rb, db = trajectories[1]
    if (ra['input_sha256'] != rb['input_sha256'] or not np.array_equal(ya[0], yb[0]) or ta[0] != tb[0]):
        raise ValueError('same physical inputs and initial state/time required')
    if ra['methods'] != rb['methods'] or ra['discretization'] != rb['discretization']:
        raise ValueError('same baseline methods and discretization required')
    if ra['physical_source_sha256'] != rb['physical_source_sha256']:
        raise ValueError('physical input/constants implementation differs')
    if set(ra['source_sha256']) != set(rb['source_sha256']):
        raise ValueError('producer graph differs; use a separately reviewed refactor audit')
    derivative_owners = ('boat_action_stencil.py', 'boat_momentum_tangent.py',
        'boat_incident_tangent.py', 'cuda_directional.py', 'cuda_directional.cu',
        'cuda_reverse_gradient.py', 'cuda_scalar_reverse.cu')
    derivative_changes = {}
    for name in sorted(ra['source_sha256'].keys()-{'boat_moving.py'}):
        hashes = tuple(report['source_sha256'][name] for report in (ra, rb))
        if hashes[0] != hashes[1]:
            if name not in derivative_owners:
                raise ValueError('physical helper differs between pilots: '+name)
            derivative_changes[name] = dict(reference=hashes[0], candidate=hashes[1])
    if derivative_changes and not derivative_change_note:
        raise ValueError('derivative helper changed; explicit reviewed implementation change note required')
    driver_changed = ra['source_sha256']['boat_moving.py'] != rb['source_sha256']['boat_moving.py']
    if driver_changed and not driver_change_note:
        raise ValueError('driver source changed; explicit reviewed change note required')
    integration_changes = {}
    for name in ('rtol', 'absolute_tolerance_by_component', 'max_step_seconds', 'first_step_seconds'):
        if ra['controls'][name] != rb['controls'][name]:
            integration_changes[name] = dict(reference=ra['controls'][name], candidate=rb['controls'][name])
    if integration_changes and not integration_change_note:
        raise ValueError('integration controls differ; explicit reviewed change note required')
    end = min(ta[-1], tb[-1])
    exact = np.intersect1d(ta[ta <= end], tb[tb <= end])
    # Evaluate both owned accepted-step polynomials, never a trial stage.
    grid = np.unique(np.r_[ta[ta <= end], tb[tb <= end], np.linspace(ta[0], end, 501)])
    heading = attitude_degrees(ya[0, 3:7])[2]
    result = dict(qualified=False,
        scope='Shared accepted motion comparison using owned RK dense polynomials; not an independently converged or open-water reference.',
        inputs=inputs, heading_frame='Fixed shared initial boat heading; surge forward, sway port, heave world up.',
        heading_degrees=float(heading), overlap_seconds=float(end-ta[0]),
        required_accuracy_duration_seconds=ACCURACY_DURATION_SECONDS,
        covers_required_duration=bool(end-ta[0] >= ACCURACY_DURATION_SECONDS),
        driver_changed=driver_changed, reviewed_driver_change_note=driver_change_note,
        integration_control_changes=integration_changes, reviewed_integration_change_note=integration_change_note,
        derivative_changes=derivative_changes, reviewed_derivative_change_note=derivative_change_note,
        performance_comparison=performance_comparison(ra, rb),
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_motion_compare.py', 'boat_motion_metrics.py', 'boat_benchmark_contract.py')},
        samples=[], limitations=[
            'Matching this unqualified finer run does not establish physical accuracy.',
            'Subminute comparisons are screening evidence only; requested motion accuracy applies throughout 60 seconds.',
            'Short weakly excited motion cannot qualify behavior over representative wave/roll periods.',
            'Accepted-step polynomials are integrator dense output, not extra temporal accuracy.',
            'Sampling them is a numerical comparison, not a rigorous continuous-time error bound.',
            'Failed-run totals and differently achieved horizons are not a matched full-run speedup.'])
    saved = {}
    for label, times in (('exact_common_knots', exact), ('dense_common_grid', grid)):
        a, b = body_poses(ta, ya, times, dense_q=da), body_poses(tb, yb, times, dense_q=db)
        translation, rotation = motion_errors(a, b, heading)
        maximum_translation = np.max(np.abs(translation), axis=0)
        maximum_rotation = np.max(np.abs(rotation), axis=0)
        result['samples'].append(dict(kind=label, count=len(times),
            first_time_seconds=float(times[0]), last_time_seconds=float(times[-1]),
            maximum_translation_error_m=maximum_translation.tolist(),
            maximum_roll_pitch_yaw_error_degrees=maximum_rotation.tolist(),
            within_translation_target=(maximum_translation <= POSITION_ERROR_TARGET_METERS).tolist(),
            within_angular_target=(maximum_rotation <= ANGLE_ERROR_TARGET_DEGREES).tolist()))
        saved[label+'_times'] = times
        saved[label+'_translation_error_m'] = translation
        saved[label+'_roll_pitch_yaw_error_degrees'] = rotation
    result['reference_maximum_absolute_roll_degrees'] = float(max(
        abs(attitude_degrees(q)[0]) for q in body_poses(ta, ya, grid, dense_q=da)[:, 3:7]))
    np.savez(output.with_suffix('.npz'), **saved)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--candidate', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--driver-change-note')
    parser.add_argument('--integration-change-note')
    parser.add_argument('--derivative-change-note')
    run(**vars(parser.parse_args()))
