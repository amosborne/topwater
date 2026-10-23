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
from scipy.integrate import RK23, RK45
from scipy.spatial.transform import Rotation

from topwater import InitialState, sea_states
from topwater.frames import rotation_matrix
from topwater.hydrostatics import FloodingError, LossOfContactError
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from tests.reference_boat import GRAVITY, WATER_DENSITY, reference_body, equilibrium
from .boat_geometry_template import evaluate_template, InvalidTemplate, InvalidFan
from .boat_world_boundary import construct, body_material_samples, interpolate_body_trace
from .boat_world_field import current_field, surface_rates, refine_boundary
from .body_pressure import evaluate_body_pressure, quadrature
from .boat_ale_pressure import differentiate
from .boat_motion_metrics import (attitude_degrees, ANGLE_ERROR_TARGET_DEGREES,
                                 POSITION_ERROR_TARGET_METERS, ACCURACY_DURATION_SECONDS)
from .boat_operator_reuse import LocalOperator, InvalidReuse, LIMITS
from .boat_past_operator import PastLocalOperator, LIMITS as PAST_LIMITS
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(output, duration=.02, refinement=1, max_rhs=13, budget_seconds=90., max_step=.01,
        quadrature_order=3, operator_reuse=None, charted_pressure=False, common_pressure_chart=False,
        action_coupling=False, initial_xy=None, action_tangent=None, integrator='RK23', tolerance_scale=1.,
        momentum_tangent=None, integer_geometry=True, shape_gradient=None):
    from .cuda_operator import prepare_energy, cp
    if action_tangent is None:
        action_tangent = action_coupling
    if momentum_tangent is None:
        momentum_tangent = action_tangent
    if shape_gradient is None:
        shape_gradient = 'reverse' if action_coupling else 'forward'
    if isinstance(operator_reuse, bool):
        operator_reuse = 'one-direction' if operator_reuse else None
    if (not np.isfinite(duration) or duration <= 0 or refinement not in (0, 1, 2)
            or max_rhs < 4 or not np.isfinite(budget_seconds) or budget_seconds <= 0
            or not np.isfinite(max_step) or max_step <= 0 or quadrature_order not in (3, 6, 12)
            or operator_reuse not in (None, 'one-direction', 'past-linear', 'past-curved')):
        raise ValueError('positive bounded duration/budget, >=4 RHS and refinement 0..2 required')
    if common_pressure_chart and not charted_pressure:
        raise ValueError('common pressure chart requires charted pressure')
    if action_coupling and (operator_reuse or common_pressure_chart):
        raise ValueError('experimental action coupling requires fresh operators and one owned chart')
    if initial_xy is not None and (np.shape(initial_xy) != (2,) or not np.isfinite(initial_xy).all()):
        raise ValueError('initial_xy must contain two finite world coordinates')
    if action_tangent and not action_coupling:
        raise ValueError('analytic action tangent requires action coupling')
    if momentum_tangent and not action_tangent:
        raise ValueError('exact momentum tangent requires analytic action tangent')
    if shape_gradient not in ('forward', 'reverse') or (shape_gradient == 'reverse' and not action_coupling):
        raise ValueError('reverse shape gradient requires action coupling')
    if integrator not in ('RK23', 'RK45'):
        raise ValueError('integrator must be RK23 or RK45')
    if not np.isfinite(tolerance_scale) or tolerance_scale <= 0:
        raise ValueError('positive finite tolerance scale required')
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    report = dict(complete=False, qualified=False,
        scope='Coupled nonlinear six-DOF and free-surface state in a 12x6x2 m diagnostic tank; no open-water or temporal/spatial qualification',
        scipy_version=scipy.__version__, requested_duration_seconds=duration,
        motion_requirement=dict(roll_pitch_yaw_error_target_degrees=ANGLE_ERROR_TARGET_DEGREES,
            surge_sway_heave_error_target_meters=POSITION_ERROR_TARGET_METERS,
            accuracy_duration_seconds=ACCURACY_DURATION_SECONDS,
            error_scope='Motion history throughout the full required duration, not endpoint only.',
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
            integrator=integrator, dense_output='Saved SciPy accepted-step RkDenseOutput coefficients',
            action_coupling=action_coupling,
            action_tangent=action_tangent,
            momentum_tangent=momentum_tangent,
            integer_geometry=bool(integer_geometry),
            shape_gradient=shape_gradient,
            initial_horizontal_position_m=None if initial_xy is None else list(initial_xy),
            common_pressure_chart=common_pressure_chart,
            pressure_derivative=('charted nodal with explicit material fallback on chart failure'
                                 if charted_pressure else 'material quadrature'),
            derivative_step_seconds=.001, rtol=.001*tolerance_scale, tolerance_scale=tolerance_scale,
            first_step_seconds=min(max_step, duration),
            minimum_chart_derivative_step_seconds=.001/64,
            max_step_seconds=max_step, topology_templates_retained=4),
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_moving.py', 'boat_geometry_template.py', 'boat_world_boundary.py',
                         'boat_parent_clip_proof.py', 'boat_integer_intersections.py', 'boat_world_field.py', 'body_pressure.py',
                         'canoe.py', 'incident_baseline_benchmark.py',
                         'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh', 'boat_operator_reuse.py',
                         'boat_past_operator.py', 'boat_ale_pressure.py', 'boat_motion_metrics.py')})
    if action_coupling:
        report['source_sha256'].update({name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_action_stencil.py', 'boat_body_action.py', 'boat_parent_action.py',
                'cuda_shape_gradient.py', 'cuda_shape.cu', 'cuda_dual.cuh', 'cuda_wide_dual.cuh')})
        report['controls'].update(pressure_derivative='diagnostic only, matched action offsets',
            derivative_step_seconds=1e-5/2, minimum_chart_derivative_step_seconds=1e-5/2**13,
            action_difference_policy='Compute h and h/2; use h/2 without a claimed error bound.')
    if action_tangent:
        report['source_sha256'].update({name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_chart_tangent.py', 'boat_action_tangent.py')})
        report['controls']['action_difference_policy'] = 'Analytic current shape derivatives; h/h2 time-offset fields, use h/2 without a claimed error bound.'
    if momentum_tangent:
        report['source_sha256'].update({name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_momentum_tangent.py', 'boat_incident_tangent.py',
                         'cuda_directional.py', 'cuda_directional.cu')})
        report['controls'].update(action_difference_policy='Analytic current shape and mixed-field time derivatives.',
            pressure_derivative='diagnostic only, exact current nodal ALE rate',
            derivative_step_seconds=None, minimum_chart_derivative_step_seconds=None)
    if shape_gradient == 'reverse':
        report['source_sha256'].update({name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('cuda_reverse_gradient.py', 'cuda_scalar_reverse.cu')})
    started = perf_counter()
    times, states, dense_q, solver = [], [], [], None
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
        _, arrays, metadata = construct(path, export_template=True, integer_predicates=integer_geometry, **inputs)
        report['geometry'].append(dict(path=str(path), seconds=perf_counter()-before,
            sha256=sha256(path.with_suffix('.npz').read_bytes()).hexdigest()))
        templates.insert(0, arrays)
        del templates[4:]
        return arrays, metadata

    def geometry(pose, eta, psi, label, chart=None):
        elapsed_check()
        before = perf_counter()
        reasons, last_error = [], None
        for index, template in enumerate(templates if chart is None else [chart]):
            try:
                arrays, diagnostics = evaluate_template(template, pose=pose,
                    parent_eta=eta, parent_psi=psi)
            except InvalidTemplate as error:
                reasons.append(str(error))
                last_error = error
                continue
            if chart is None:
                templates.insert(0, templates.pop(index))
            mode = 'current guarded chart'
            break
        else:
            if chart is not None:
                raise type(last_error)('pressure offsets leave the base chart: '+str(reasons)) from last_error
            arrays, metadata = cold(body=body, pose=pose,
                parent_xy=initial_geometry['parent_nodes'][:, :2], parents=initial_geometry['parents'],
                parent_eta=eta, parent_psi=psi)
            diagnostics = dict(reference_seconds=metadata['seconds'])
            mode = 'exact chart refresh'
            if action_tangent:
                # Cold rational intersections and floating chart replay can
                # differ by ulps. Choose ONE replay before building either the
                # primal operator or its tangent; never relax their ownership.
                arrays, replay = evaluate_template(templates[0], pose=pose,
                    parent_eta=eta, parent_psi=psi)
                diagnostics.update(replay)
                mode = 'exact chart refresh with owned replay'
        for _ in range(refinement):
            arrays = refine_boundary(arrays, region='body')
        row['geometry'].append(dict(label=label, mode=mode, seconds=perf_counter()-before,
            rejected_charts=reasons, diagnostics=diagnostics, faces=len(arrays['faces'])))
        return arrays

    def field(arrays, twist, time, acceleration_columns, model, *, retain_assembly=False):
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
        return value, assembly if retain_assembly else energy

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
        base_context = dict(arrays=arrays, chart=templates[0], centers={})
        rotation = rotation_matrix(quaternion)
        failed_offset = None

        def stencil(model, use_chart=charted_pressure, step=.001, context=None):
            nonlocal failed_offset
            context = base_context if context is None else context
            arrays, selected_chart, centers = context['arrays'], context['chart'], context['centers']
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
                moved_pose = InitialState(pose.position+offset*twist[:3], orientation)
                moved_eta, moved_psi = eta+offset*eta_rate, psi+offset*psi_rate
                failed_offset = moved_pose, moved_eta, moved_psi, offset
                try:
                    moved.append(geometry(moved_pose, moved_eta, moved_psi, f'offset-{sign}',
                        chart=selected_chart if use_chart else None))
                except ValueError as error:
                    if model is not None:
                        raise InvalidReuse('predicted direction has invalid offset geometry: '+str(error)) from error
                    raise
            if model is not None:
                row['reuse_direction'] = model.check_direction(*moved)
            if not use_chart:
                samples, groups = body_material_samples(arrays, quadrature()[0])
                hints = np.broadcast_to(np.arange(len(arrays['body_faces']))[:, None], groups.shape)
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
        if action_coupling:
            from .boat_action_stencil import evaluate as action_stencil
            def tangent_geometry():
                from .boat_chart_tangent import chart_tangent
                from .boat_action_tangent import refine_tangent
                coarse, derivative = chart_tangent(base_context['chart'], pose=pose, parent_eta=eta, parent_psi=psi)
                return refine_tangent(coarse, derivative['point_jacobian'], derivative['parent_map_jacobian'], refinement)
            current, assembly = field(arrays, twist, time, True, None, retain_assembly=True)
            result = action_stencil(body, arrays, current, assembly, state, values, time,
                lambda current_pose, height, potential, label: geometry(current_pose, height,
                    potential, label, chart=base_context['chart']),
                lambda moved, current_time: field(moved, twist, current_time, False, None),
                geometry_tangent=tangent_geometry if action_tangent else None,
                momentum_tangent=momentum_tangent, shape_gradient=shape_gradient)
            row['operator_mode'] = 'fresh'
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
                    if (common_pressure_chart and isinstance(error, InvalidFan)
                            and 'common_chart_attempt' not in row):
                        # The physical arrangement passed. Try a different
                        # triangulation, including a wholly NEW CENTER solve.
                        # Never difference fields from different charts.
                        old_templates = list(templates)
                        moved_pose, moved_eta, moved_psi, offset = failed_offset
                        attempt = dict(accepted=False, signed_offset_seconds=offset,
                            reason=str(error), old_chart_key=str(base_context['chart']['node_chart_key']))
                        row['common_chart_attempt'] = attempt
                        repair_started = perf_counter()
                        try:
                            chart, _ = cold(body=body, pose=moved_pose,
                                parent_xy=initial_geometry['parent_nodes'][:, :2], parents=initial_geometry['parents'],
                                parent_eta=moved_eta, parent_psi=moved_psi)
                            rebuilt = geometry(pose, eta, psi, 'common-chart-base', chart=chart)
                            new_context = dict(arrays=rebuilt, chart=chart, centers={})
                            repaired = stencil(None, step=.001/2**halving, context=new_context)
                            if not np.isfinite(np.r_[repaired[1], repaired[2], repaired[0].acceleration]).all():
                                raise FloatingPointError('nonfinite common-chart RHS')
                        except (FloodingError, LossOfContactError):
                            raise
                        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as rejected:
                            attempt['rejection'] = str(rejected)
                            templates[:] = old_templates
                        else:
                            base_context, result = new_context, repaired
                            attempt.update(accepted=True, new_chart_key=str(chart['node_chart_key']),
                                old_faces=len(arrays['faces']), new_faces=len(rebuilt['faces']))
                            operator_model, previous_stencil = None, None
                            if audit_candidate is not None:
                                row['refresh_audit'] = dict(available=False,
                                    reason='common chart changed the discrete trial space')
                                audit_candidate = None
                            break
                        finally:
                            attempt['seconds'] = perf_counter()-repair_started
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
        if initial_xy is None:
            initial_geometry, metadata = cold()
        else:
            rest = equilibrium()[0]
            initial_geometry, metadata = cold(pose=InitialState(np.r_[initial_xy, rest.position[2]], rest.orientation))
        report['fixture_sha256'] = metadata['fixture_sha256']
        values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
        size = len(initial_geometry['parent_nodes'])
        initial = np.r_[initial_geometry['pose_position'], initial_geometry['pose_orientation'],
            [.03, -.02, .01, .01, -.015, .02], initial_geometry['parent_nodes'][:, 2], initial_geometry['parent_psi']]
        times, states = [0.], [initial.copy()]
        # Exploratory physical-unit integration controls, not acceptance limits.
        # They must be refined before interpreting a trajectory quantitatively.
        atol = np.r_[np.full(3, .001), np.full(4, .0001), np.full(3, .001),
            np.full(3, .0001), np.full(size, .001), np.full(size, .01)]*tolerance_scale
        report['controls']['absolute_tolerance_by_component'] = atol.tolist()
        solver = (RK23 if integrator == 'RK23' else RK45)(rhs, 0., initial, duration, first_step=min(max_step, duration),
            max_step=max_step, rtol=.001*tolerance_scale, atol=atol)
        while solver.status == 'running':
            message = solver.step()
            if solver.status == 'failed':
                raise RuntimeError(message)
            times.append(float(solver.t)); states.append(solver.y.copy())
            dense_q.append(solver.dense_output().Q.copy())
            saved.update(accepted_times=np.array(times), accepted_states=np.array(states), accepted_dense_q=np.array(dense_q))
            np.savez(output.with_suffix('.npz'), **saved)
        if not np.array_equal(saved['last_rhs_state'], solver.y):
            raise RuntimeError('accepted endpoint lacks matching FSAL diagnostics')
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        saved.update(accepted_times=np.array(times), accepted_states=np.array(states), accepted_dense_q=np.array(dense_q))
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
        if solver is not None and report['complete']:
            attempts = (solver.nfev-1)//solver.n_stages
            report['integration_stats'] = dict(rhs_calls=solver.nfev, attempted_steps=attempts,
                rejected_steps=attempts-max(0, len(times)-1), stages_per_attempt=solver.n_stages)
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
    parser.add_argument('--common-pressure-chart', action='store_true')
    parser.add_argument('--action-coupling', action='store_true')
    parser.add_argument('--initial-xy', nargs=2, type=float)
    parser.add_argument('--action-tangent', action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument('--integrator', choices=('RK23', 'RK45'), default='RK23')
    parser.add_argument('--tolerance-scale', type=float, default=1.)
    parser.add_argument('--momentum-tangent', action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument('--integer-geometry', action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('--shape-gradient', choices=('forward', 'reverse'), default=None)
    result = run(**vars(parser.parse_args()))
    print(json.dumps({key: result[key] for key in
        ('complete', 'seconds', 'accepted_simulated_seconds', 'accepted_steps', 'cold_geometry_count')}))
