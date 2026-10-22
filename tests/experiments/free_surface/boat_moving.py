"""Bounded coupled nonlinear six-DOF boat evolution in the diagnostic tank.

Actual observed sea, evolving WORLD surface state and Bernoulli pressure.
This is not a qualified open-water model or a simulation throughput gate.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import scipy
from scipy.integrate import RK23
from scipy.spatial.transform import Rotation

from topwater import InitialState, sea_states
from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from tests.reference_boat import GRAVITY, WATER_DENSITY, reference_body
from .boat_geometry_template import evaluate_template, InvalidTemplate
from .boat_world_boundary import construct, body_material_samples, interpolate_body_trace
from .boat_world_field import current_field, surface_rates, refine_boundary
from .body_pressure import evaluate_body_pressure, quadrature
from .boat_ale_pressure import differentiate
from .boat_motion_metrics import (attitude_degrees, ANGLE_ERROR_TARGET_DEGREES,
                                 POSITION_ERROR_TARGET_METERS)
from .boat_operator_reuse import LocalOperator, InvalidReuse, LIMITS
from .boat_past_operator import PastLocalOperator, LIMITS as PAST_LIMITS
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(output, duration=.02, refinement=1, max_rhs=13, budget_seconds=90., max_step=.01,
        quadrature_order=3, operator_reuse=None, charted_pressure=False):
    from .cuda_operator import prepare_energy, cp
    if isinstance(operator_reuse, bool):
        operator_reuse = 'one-direction' if operator_reuse else None
    if (not np.isfinite(duration) or duration <= 0 or refinement not in (0, 1, 2)
            or max_rhs < 4 or not np.isfinite(budget_seconds) or budget_seconds <= 0
            or not np.isfinite(max_step) or max_step <= 0 or quadrature_order not in (3, 6, 12)
            or operator_reuse not in (None, 'one-direction', 'past-linear', 'past-curved')):
        raise ValueError('positive bounded duration/budget, >=4 RHS and refinement 0..2 required')
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    report = dict(complete=False, qualified=False,
        scope='Coupled nonlinear six-DOF and free-surface state in a 12x6x2 m diagnostic tank; no open-water or temporal/spatial qualification',
        scipy_version=scipy.__version__, requested_duration_seconds=duration,
        motion_requirement=dict(roll_pitch_yaw_error_target_degrees=ANGLE_ERROR_TARGET_DEGREES,
            surge_sway_heave_error_target_meters=POSITION_ERROR_TARGET_METERS,
            translation_comparison_frame='Common reference heading for horizontal surge/sway; world up for heave.',
            operating_scope='Ordinary motion tolerable for a person aboard; no fixed roll-angle cutoff.',
            status='Requested numerical motion convergence and operating scope; not yet validated.'),
        state_layout=['position_world_m[3]', 'orientation_wxyz[4]',
            'velocity_world_m_per_s[3]', 'angular_velocity_world_rad_per_s[3]',
            'parent_eta_m[N]', 'parent_psi_m2_per_s[N]'],
        boundary_refinement=dict(region='body', levels=refinement), rhs=[], geometry=[],
        operator_policy=dict(enabled=operator_reuse is not None, kind=operator_reuse, quadrature_order=quadrature_order,
            limits=(LIMITS if operator_reuse == 'one-direction' else PAST_LIMITS) if operator_reuse else None,
            scope='Past-only local model; atomic whole-RHS fallback. Geometry limits are heuristics, not force-error bounds. Stage history can affect refresh decisions.'),
        controls=dict(max_rhs=max_rhs, soft_elapsed_budget_seconds=budget_seconds,
            pressure_derivative=('charted nodal with explicit material fallback on chart failure'
                                 if charted_pressure else 'material quadrature'),
            derivative_step_seconds=.001, rtol=.001, first_step_seconds=min(max_step, duration),
            minimum_chart_derivative_step_seconds=.001/64,
            max_step_seconds=max_step, topology_templates_retained=4),
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_moving.py', 'boat_geometry_template.py', 'boat_world_boundary.py',
                         'boat_parent_clip_proof.py', 'boat_world_field.py', 'body_pressure.py',
                         'canoe.py', 'incident_baseline_benchmark.py',
                         'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh', 'boat_operator_reuse.py',
                         'boat_past_operator.py', 'boat_ale_pressure.py', 'boat_motion_metrics.py')})
    started = perf_counter()
    times, states, solver = [], [], None
    templates, cold_count, operator_model = [], 0, None
    previous_stencil = None
    saved = {}
    body = reference_body()

    def elapsed_check():
        if perf_counter()-started > budget_seconds:
            raise RuntimeError('bounded moving-case budget exhausted between operations')

    def cold(**inputs):
        nonlocal cold_count
        elapsed_check()
        path = output.with_name(output.stem+f'-geometry-{cold_count}.json')
        cold_count += 1
        before = perf_counter()
        _, arrays, metadata = construct(path, export_template=True, **inputs)
        report['geometry'].append(dict(path=str(path), seconds=perf_counter()-before,
            sha256=sha256(path.with_suffix('.npz').read_bytes()).hexdigest()))
        templates.insert(0, arrays)
        del templates[4:]
        return arrays, metadata

    def geometry(pose, eta, psi, label, chart=None):
        elapsed_check()
        before = perf_counter()
        reasons = []
        for index, template in enumerate(templates if chart is None else [chart]):
            try:
                arrays, diagnostics = evaluate_template(template, pose=pose,
                    parent_eta=eta, parent_psi=psi)
            except InvalidTemplate as error:
                reasons.append(str(error))
                continue
            if chart is None:
                templates.insert(0, templates.pop(index))
            mode = 'current guarded chart'
            break
        else:
            if chart is not None:
                raise InvalidTemplate('pressure offsets leave the base chart: '+str(reasons))
            arrays, metadata = cold(body=body, pose=pose,
                parent_xy=initial_geometry['parent_nodes'][:, :2], parents=initial_geometry['parents'],
                parent_eta=eta, parent_psi=psi)
            diagnostics = dict(reference_seconds=metadata['seconds'])
            mode = 'exact chart refresh'
        for _ in range(refinement):
            arrays = refine_boundary(arrays, region='body')
        row['geometry'].append(dict(label=label, mode=mode, seconds=perf_counter()-before,
            rejected_charts=reasons, diagnostics=diagnostics, faces=len(arrays['faces'])))
        return arrays

    def field(arrays, twist, time, acceleration_columns, model):
        elapsed_check()
        before = perf_counter()
        prediction = None
        if model is None:
            assembly = prepare_energy(arrays['faces'], TrianglePanels(arrays['points'][arrays['faces']], 1),
                                      order=quadrature_order)
            if operator_reuse and (not assembly.uniform_quadrature or assembly.quadrature_order != quadrature_order):
                raise ValueError('operator reuse requires the declared uniform quadrature owner')
            energy = assembly.energy
        else:
            energy, prediction = model.predict(arrays, time)
        assembled = perf_counter()
        value = current_field(arrays, energy, values, twist,
                              acceleration_columns=acceleration_columns, time=time)
        row['fields'].append(dict(mode='fresh' if model is None else 'predicted', prediction=prediction,
            energy_seconds=assembled-before,
            solve_seconds=perf_counter()-assembled, relative_residual=value['relative_field_residual']))
        if not np.isfinite(value['relative_field_residual']) or value['relative_field_residual'] > 1e-9:
            if model is not None:
                raise InvalidReuse('predicted field exceeds the algebraic residual budget')
            raise ValueError('current field exceeds the audited algebraic residual budget')
        return value, energy

    def rhs(time, state):
        nonlocal row, operator_model, previous_stencil
        elapsed_check()
        if len(report['rhs']) >= max_rhs:
            raise RuntimeError('bounded moving-case RHS count exhausted')
        before = perf_counter()
        row = dict(time=float(time), complete=False, geometry=[], fields=[])
        report['rhs'].append(row)
        print(json.dumps(dict(event='rhs-start', index=len(report['rhs']), time=float(time))), flush=True)
        if not np.isfinite(state).all() or np.linalg.norm(state[3:7]) == 0:
            raise ValueError('nonfinite state or zero quaternion')
        quaternion = state[3:7]/np.linalg.norm(state[3:7])
        pose = InitialState(state[:3], quaternion)
        twist, eta, psi = state[7:13], state[13:13+size], state[13+size:]
        saved['attempted_rhs_state'], saved['attempted_rhs_time'] = state.copy(), float(time)
        arrays = geometry(pose, eta, psi, 'base')
        selected_chart = templates[0]
        samples, groups = body_material_samples(arrays, quadrature()[0])
        hints = np.broadcast_to(np.arange(len(arrays['body_faces']))[:, None], groups.shape)
        rotation = rotation_matrix(quaternion)
        centers = {}

        def stencil(model, use_chart=charted_pressure, step=.001):
            if id(model) not in centers:
                current, center_energy = field(arrays, twist, time, True, model)
                eta_rate, psi_rate, surface = surface_rates(arrays['points'],
                    arrays['faces'][arrays['free_faces']], arrays['free_nodes'], current['parent_map'],
                    psi, current['flux'])
                centers[id(model)] = current, center_energy, eta_rate, psi_rate, surface
            current, center_energy, eta_rate, psi_rate, surface = centers[id(model)]
            moved = []
            for sign in (1, -1):
                offset = sign*step
                orientation = (Rotation.from_rotvec(offset*twist[3:])
                    *Rotation.from_quat(quaternion, scalar_first=True)).as_quat(scalar_first=True)
                try:
                    moved.append(geometry(InitialState(pose.position+offset*twist[:3], orientation),
                        eta+offset*eta_rate, psi+offset*psi_rate, f'offset-{sign}',
                        chart=selected_chart if use_chart else None))
                except ValueError as error:
                    if model is not None:
                        raise InvalidReuse('predicted direction has invalid offset geometry: '+str(error)) from error
                    raise
            if model is not None:
                row['reuse_direction'] = model.check_direction(*moved)
            traces, energies = [], [center_energy]
            for sign, offset_arrays in zip((1, -1), moved):
                moved_field, energy = field(offset_arrays, twist, time+sign*step, False, model)
                energies.append(energy)
                traces.append(moved_field['potential'] if use_chart else
                    interpolate_body_trace(offset_arrays, moved_field['potential'], samples, groups,
                                           face_hints=hints))
            if use_chart:
                try:
                    rate, mesh_velocity, diagnostics = differentiate(arrays, *moved, *traces,
                        step, pose.position, twist)
                except ValueError as error:
                    if model is not None:
                        raise InvalidReuse(str(error)) from error
                    raise
            else:
                rate = (traces[0]-traces[1])/(2*step)
                mesh_velocity = twist[:3]+np.cross(twist[3:], arrays['points']-pose.position)
                diagnostics = {}
            pressure = evaluate_body_pressure(arrays['points'], arrays['faces'], arrays['body_faces'],
                pose.position, body.mass, rotation@body.inertia_cg@rotation.T, twist, current['potential'],
                rate, mesh_velocity,
                current['acceleration_potentials'], density=WATER_DENSITY, gravity=GRAVITY)
            return (pressure, eta_rate, psi_rate, surface, [arrays, *moved], energies,
                    dict(representation='charted nodal' if use_chart else 'material quadrature',
                         step_seconds=step, diagnostics=diagnostics))

        result, audit_candidate = None, None
        refresh = (operator_reuse in ('past-linear', 'past-curved') and operator_model is not None
                   and abs(time-operator_model.time) >= PAST_LIMITS['refresh_age_seconds']-1e-14)
        if operator_model is not None:
            try:
                candidate = stencil(operator_model,
                    use_chart=operator_model.pressure_derivative_representation == 'charted nodal',
                    step=operator_model.derivative_step_seconds)
                if refresh:
                    audit_candidate = candidate
                    row['refresh_reason'] = 'scheduled age refresh; valid old model audited'
                else:
                    result = candidate
                    row['operator_mode'] = 'reused'
                    row['operator_anchor_rhs_index'] = operator_model.rhs_index
            except (InvalidReuse, np.linalg.LinAlgError, FloatingPointError) as error:
                row['reuse_rejection'] = str(error)
                row['refresh_audit'] = dict(available=False, reason=str(error))
                # Discard ALL candidate values, including their surface direction.
                # A fresh base determines both fresh pressure-offset geometries.
        if result is None:
            for halving in range(7 if charted_pressure else 1):
                try:
                    result = stencil(None, step=.001/2**halving)
                    break
                except InvalidTemplate as error:
                    if not charted_pressure:
                        raise
                    row.setdefault('pressure_chart_step_rejections', []).append(
                        dict(step_seconds=.001/2**halving, reason=str(error)))
            if result is None:
                # No common nodal chart at these steps: explicitly use strictly
                # wet material samples. Never mix the two offset representations.
                row['pressure_derivative_fallback'] = 'no common chart after six step halvings'
                result = stencil(None, use_chart=False)
            row['operator_mode'] = 'fresh'
        pressure, eta_rate, psi_rate, surface, stencil_arrays, stencil_energies, derivative = result
        row['pressure_derivative_representation'] = derivative['representation']
        row['pressure_derivative'] = derivative['diagnostics']
        row['pressure_derivative_step_seconds'] = derivative['step_seconds']
        if audit_candidate is not None:
            old_pressure, old_eta, old_psi = audit_candidate[:3]
            row['refresh_audit'] = dict(available=True, anchor_rhs_index=operator_model.rhs_index,
                candidate_pressure_derivative=audit_candidate[6]['representation'],
                fresh_pressure_derivative=derivative['representation'],
                candidate_derivative_step_seconds=audit_candidate[6]['step_seconds'],
                fresh_derivative_step_seconds=derivative['step_seconds'],
                operator_only_comparison=(audit_candidate[6]['representation'] == derivative['representation']
                    and audit_candidate[6]['step_seconds'] == derivative['step_seconds']),
                hydrodynamic_load_error_world=(old_pressure.hydrodynamic_load-pressure.hydrodynamic_load).tolist(),
                acceleration_error_world=(old_pressure.acceleration-pressure.acceleration).tolist(),
                parent_height_rate_error_m_per_s=(old_eta-eta_rate).tolist(),
                parent_potential_rate_error_m2_per_s2=(old_psi-psi_rate).tolist())
        raw_quaternion = state[3:7]
        quaternion_rate = .5*np.r_[-twist[3:]@raw_quaternion[1:],
            raw_quaternion[0]*twist[3:]+np.cross(twist[3:], raw_quaternion[1:])]
        rate = np.r_[twist[:3], quaternion_rate, pressure.acceleration, eta_rate, psi_rate]
        if not np.isfinite(rate).all():
            raise FloatingPointError('nonfinite coupled boat state rate')
        # Commit a model only after its entire fresh RHS has succeeded. RK-rejected
        # stages may train later predictions; events retain that computational history.
        if operator_reuse and row['operator_mode'] == 'fresh':
            training_started = perf_counter()
            current_stencil = (stencil_arrays, stencil_energies, float(time), len(report['rhs']))
            try:
                if operator_reuse == 'one-direction':
                    operator_model = LocalOperator(stencil_arrays, stencil_energies, time, len(report['rhs']))
                elif previous_stencil is not None:
                    operator_model = PastLocalOperator(previous_stencil, current_stencil, time, len(report['rhs']),
                                                       curved=operator_reuse == 'past-curved')
                else:
                    raise InvalidReuse('no previous completed fresh stencil')
                operator_model.derivative_step_seconds = derivative['step_seconds']
                operator_model.pressure_derivative_representation = derivative['representation']
                row['operator_training'] = dict(operator_model.diagnostics,
                    derivative_step_seconds=derivative['step_seconds'],
                    pressure_derivative_representation=derivative['representation'])
            except (ValueError, np.linalg.LinAlgError, FloatingPointError) as error:
                operator_model = None
                row['operator_training_rejection'] = str(error)
            if operator_reuse in ('past-linear', 'past-curved'):
                previous_stencil = current_stencil
            row['operator_training_seconds'] = perf_counter()-training_started
        row.update(complete=True, seconds=perf_counter()-before, surface=surface,
            acceleration_world=pressure.acceleration.tolist(),
            hydrodynamic_load_world=pressure.hydrodynamic_load.tolist(), pressure=pressure.diagnostics)
        saved['last_rhs_state'], saved['last_rhs_rate'] = state.copy(), rate
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        print(json.dumps(dict(event='rhs-complete', index=len(report['rhs']), time=float(time),
                              seconds=row['seconds'], cold_geometries=cold_count)), flush=True)
        return rate

    row = None
    try:
        before = perf_counter()
        cp.cuda.get_current_stream().synchronize()
        report['gpu_context_seconds'] = perf_counter()-before
        initial_geometry, metadata = cold()
        report['fixture_sha256'] = metadata['fixture_sha256']
        values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
        size = len(initial_geometry['parent_nodes'])
        initial = np.r_[initial_geometry['pose_position'], initial_geometry['pose_orientation'],
            [.03, -.02, .01, .01, -.015, .02], initial_geometry['parent_nodes'][:, 2], initial_geometry['parent_psi']]
        times, states = [0.], [initial.copy()]
        # Exploratory physical-unit integration controls, not acceptance limits.
        # They must be refined before interpreting a trajectory quantitatively.
        atol = np.r_[np.full(3, .001), np.full(4, .0001), np.full(3, .001),
            np.full(3, .0001), np.full(size, .001), np.full(size, .01)]
        report['controls']['absolute_tolerance_by_component'] = atol.tolist()
        solver = RK23(rhs, 0., initial, duration, first_step=min(max_step, duration),
            max_step=max_step, rtol=.001, atol=atol)
        while solver.status == 'running':
            message = solver.step()
            if solver.status == 'failed':
                raise RuntimeError(message)
            times.append(float(solver.t)); states.append(solver.y.copy())
            saved.update(accepted_times=np.array(times), accepted_states=np.array(states))
            np.savez(output.with_suffix('.npz'), **saved)
        if not np.array_equal(saved['last_rhs_state'], solver.y):
            raise RuntimeError('accepted endpoint lacks matching FSAL diagnostics')
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        saved.update(accepted_times=np.array(times), accepted_states=np.array(states))
        report.update(seconds=perf_counter()-started,
            accepted_simulated_seconds=float(times[-1]) if times else 0., accepted_steps=max(0, len(times)-1),
            cold_geometry_count=cold_count, process_peak_rss_bytes=peak_rss_bytes(),
            reserved_gpu_pool_bytes=cp.get_default_memory_pool().total_bytes())
        report['operator_counts'] = dict(
            fresh_rhs=sum(r.get('operator_mode') == 'fresh' for r in report['rhs']),
            reused_rhs=sum(r.get('operator_mode') == 'reused' for r in report['rhs']),
            fresh_assemblies=sum(f['mode'] == 'fresh' for r in report['rhs'] for f in r['fields']),
            predicted_fields=sum(f['mode'] == 'predicted' for r in report['rhs'] for f in r['fields']),
            failed_reuse_attempts=sum('reuse_rejection' in r for r in report['rhs']))
        if states:
            angles = np.array([attitude_degrees(state[3:7]) for state in states])
            rolls = angles[:, 0]
            saved['accepted_roll_pitch_yaw_degrees'] = angles
            report['maximum_absolute_accepted_roll_degrees'] = float(np.max(np.abs(rolls)))
            report['maximum_accepted_quaternion_norm_error'] = float(np.max(abs(
                np.linalg.norm(np.array(states)[:, 3:7], axis=1)-1)))
        np.savez(output.with_suffix('.npz'), **saved)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--duration', type=float, default=.02)
    parser.add_argument('--refinement', type=int, default=1)
    parser.add_argument('--max-rhs', type=int, default=13)
    parser.add_argument('--budget-seconds', type=float, default=90.)
    parser.add_argument('--max-step', type=float, default=.01)
    parser.add_argument('--quadrature-order', type=int, default=3, choices=(3, 6, 12))
    parser.add_argument('--operator-reuse', nargs='?', const='one-direction', choices=('one-direction', 'past-linear', 'past-curved'))
    parser.add_argument('--charted-pressure', action='store_true')
    result = run(**vars(parser.parse_args()))
    print(json.dumps({key: result[key] for key in
        ('complete', 'seconds', 'accepted_simulated_seconds', 'accepted_steps', 'cold_geometry_count')}))
