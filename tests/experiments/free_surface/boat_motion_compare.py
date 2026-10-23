"""Compare accepted RK body motion on a common time/heading frame."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .boat_motion_metrics import (attitude_degrees, body_poses, motion_errors,
    ANGLE_ERROR_TARGET_DEGREES, POSITION_ERROR_TARGET_METERS)
from .p4_fluid_ordering import fresh


def run(reference, candidate, output, driver_change_note=None, integration_change_note=None):
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    inputs, trajectories = [], []
    for path in map(Path, (reference, candidate)):
        report = json.loads(path.read_text())
        with np.load(path.with_suffix('.npz'), allow_pickle=False) as data:
            times, states = data['accepted_times'], data['accepted_states']
            dense = data['accepted_dense_q'] if 'accepted_dense_q' in data else None
        if report['controls'].get('integrator', 'RK23') != 'RK23' and dense is None:
            raise ValueError('non-RK23 trajectory requires its owned dense coefficients')
        method = report['controls'].get('integrator', 'RK23')
        if method not in ('RK23', 'RK45') or (dense is not None and dense.shape[-1] != {'RK23': 3, 'RK45': 4}[method]):
            raise ValueError('dense coefficient degree must match the declared integrator')
        if len(times) < 2:
            raise ValueError('at least two accepted states per trajectory required')
        metadata_path = Path(report['geometry'][0]['path'])
        metadata = json.loads(metadata_path.read_text())
        inputs.append(dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest(),
            arrays_sha256=sha256(path.with_suffix('.npz').read_bytes()).hexdigest(),
            initial_geometry_metadata_sha256=sha256(metadata_path.read_bytes()).hexdigest(),
            complete=report['complete'], accepted_seconds=float(times[-1]-times[0]),
            seconds=report['seconds'], operator_counts=report['operator_counts'],
            refinement=report['boundary_refinement'],
            quadrature_order=report['operator_policy']['quadrature_order']))
        inputs[-1]['integrator'] = report['controls'].get('integrator', 'RK23')
        trajectories.append((times, states, report, metadata, dense))
    ta, ya, ra, ma, da = trajectories[0]
    tb, yb, rb, mb, db = trajectories[1]
    if (ra['fixture_sha256'] != rb['fixture_sha256'] or not np.array_equal(ya[0], yb[0])
            or ma['actual_sea_realization'] != mb['actual_sea_realization'] or ta[0] != tb[0]):
        raise ValueError('same physical fixture, initial state/time and sea realization required')
    physics = ('boat_world_boundary.py', 'boat_parent_clip_proof.py', 'boat_world_field.py',
               'boat_geometry_template.py', 'body_pressure.py', 'incident_baseline_benchmark.py',
               'boat_ale_pressure.py', 'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh',
               'canoe.py')
    for name in ('action_coupling', 'action_tangent'):
        if ra['controls'].get(name, False) != rb['controls'].get(name, False):
            raise ValueError('action policy differs: '+name)
    if ra['controls'].get('action_coupling', False):
        physics += ('boat_action_stencil.py', 'boat_body_action.py', 'boat_parent_action.py',
                    'cuda_shape_gradient.py', 'cuda_shape.cu', 'cuda_dual.cuh', 'cuda_wide_dual.cuh')
    if ra['controls'].get('action_tangent', False):
        physics += ('boat_chart_tangent.py', 'boat_action_tangent.py')
    for name in physics:
        if ra['source_sha256'][name] != rb['source_sha256'][name]:
            raise ValueError('physical helper differs between pilots: '+name)
    driver_changed = ra['source_sha256']['boat_moving.py'] != rb['source_sha256']['boat_moving.py']
    if driver_changed and not driver_change_note:
        raise ValueError('driver source changed; explicit reviewed change note required')
    integration_changes = {}
    methods = tuple(report['controls'].get('integrator', 'RK23') for report in (ra, rb))
    if methods[0] != methods[1]:
        integration_changes['integrator'] = dict(reference=methods[0], candidate=methods[1])
    for name in ('rtol', 'absolute_tolerance_by_component', 'max_step_seconds', 'first_step_seconds'):
        if ra['controls'][name] != rb['controls'][name]:
            integration_changes[name] = dict(reference=ra['controls'][name], candidate=rb['controls'][name])
    if integration_changes and not integration_change_note:
        raise ValueError('integration controls differ; explicit reviewed change note required')
    for name in ('derivative_step_seconds', 'pressure_derivative', 'minimum_chart_derivative_step_seconds'):
        if ra['controls'][name] != rb['controls'][name]:
            raise ValueError('integration control differs: '+name)
    if ra['controls'].get('common_pressure_chart', False) != rb['controls'].get('common_pressure_chart', False):
        raise ValueError('common pressure chart policy differs')
    if ra['boundary_refinement']['region'] != rb['boundary_refinement']['region']:
        raise ValueError('refinement region differs')
    for name in ('enabled', 'kind', 'limits'):
        if ra['operator_policy'][name] != rb['operator_policy'][name]:
            raise ValueError('operator reuse policy differs: '+name)
    end = min(ta[-1], tb[-1])
    exact = np.intersect1d(ta[ta <= end], tb[tb <= end])
    # Evaluate both owned accepted-step polynomials, never a trial stage.
    grid = np.unique(np.r_[ta[ta <= end], tb[tb <= end], np.linspace(ta[0], end, 501)])
    heading = attitude_degrees(ya[0, 3:7])[2]
    result = dict(qualified=False,
        scope='Shared accepted motion comparison using owned RK dense polynomials; not an independently converged or open-water reference.',
        inputs=inputs, heading_frame='Fixed shared initial boat heading; surge forward, sway port, heave world up.',
        heading_degrees=float(heading), overlap_seconds=float(end-ta[0]),
        driver_changed=driver_changed, reviewed_driver_change_note=driver_change_note,
        integration_control_changes=integration_changes, reviewed_integration_change_note=integration_change_note,
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_motion_compare.py', 'boat_motion_metrics.py')},
        samples=[], limitations=[
            'Matching this unqualified finer run does not establish physical accuracy.',
            'Short weakly excited motion cannot qualify behavior over representative wave/roll periods.',
            'Accepted-step polynomials are integrator dense output, not extra temporal accuracy.',
            'Sampling them is a numerical comparison, not a rigorous continuous-time error bound.',
            'Combined body/quadrature change does not isolate either source of sensitivity.',
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
    run(**vars(parser.parse_args()))
