"""Single-path nonlinear six-DOF action baseline in the diagnostic tank.

Boat, mass properties, initial state and actual wave phases are input data.
The finite tank and unresolved force/spatial accuracy are not open-water
qualification. Independent derivative references live in verification harnesses.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import scipy
from scipy.integrate import RK23

from topwater import InitialState, SimulationInputs, sea_states
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from .boat_geometry_template import evaluate_template, InvalidTemplate
from .boat_world_boundary import construct
from .boat_world_field import current_field, refine_boundary
from .boat_chart_tangent import chart_tangent
from .boat_action_tangent import refine_tangent
from .boat_action_stencil import evaluate as action_rhs
from .boat_motion_metrics import (attitude_degrees, ANGLE_ERROR_TARGET_DEGREES,
                                 POSITION_ERROR_TARGET_METERS, ACCURACY_DURATION_SECONDS)
from .incident_baseline_benchmark import coefficients
from .boat_integration_controls import (RELATIVE_TOLERANCE, BODY_ABSOLUTE_TOLERANCE,
    SURFACE_HEIGHT_ABSOLUTE_TOLERANCE_M, SURFACE_POTENTIAL_ABSOLUTE_TOLERANCE_M2_PER_S, ERROR_NORM)
from .boat_benchmark_contract import measurement_environment, array_identity
from .p4_fluid_ordering import fresh


METHODS = dict(integrator='RK23', geometry='exact integer with guarded chart replay',
    geometry_derivative='analytic chart', momentum_derivative='exact current',
    shape_gradient='sample-local reverse', formulation='restricted action')
FIELD_RESIDUAL_LIMIT = 1e-9
SOURCES = ('boat_moving.py', 'boat_geometry_template.py', 'boat_world_boundary.py',
    'boat_parent_clip_proof.py', 'boat_integer_intersections.py', 'boat_world_field.py',
    'body_pressure.py', 'incident_baseline_benchmark.py', 'cuda_operator.py',
    'cuda_energy.cu', 'cuda_pair.cuh', 'boat_ale_pressure.py', 'boat_motion_metrics.py',
    'boat_action_stencil.py', 'boat_parent_action.py', 'boat_chart_tangent.py',
    'boat_action_tangent.py', 'boat_momentum_tangent.py', 'boat_incident_tangent.py',
    'cuda_directional.py', 'cuda_directional.cu', 'cuda_dual.cuh',
    'cuda_reverse_gradient.py', 'cuda_scalar_reverse.cu', 'cuda_panel_geometry.py',
    'boat_integration_controls.py', 'boat_benchmark_contract.py')


def run(inputs, output, duration=.02, refinement=1, max_rhs=13, budget_seconds=90.,
        max_step=.01, quadrature_order=3, tolerance_scale=1.,
        resource_status='uncontrolled', resource_context='Uncharacterized shared-machine load'):
    """Evolve supplied physical inputs using the one adopted solver path."""
    from .cuda_operator import prepare_energy, cp
    if not isinstance(inputs, SimulationInputs):
        raise TypeError('inputs must be SimulationInputs')
    if (not np.isfinite(duration) or duration <= 0 or refinement not in (0, 1, 2)
            or max_rhs < 4 or not np.isfinite(budget_seconds) or budget_seconds <= 0
            or not np.isfinite(max_step) or max_step <= 0 or quadrature_order not in (3, 6, 12)
            or not np.isfinite(tolerance_scale) or tolerance_scale <= 0):
        raise ValueError('positive bounded duration/budget/tolerance, >=4 RHS and refinement 0..2 required')
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    beam = float(np.ptp(inputs.body.mesh.vertices[:, 1]))
    sea_states.validate_envelope(inputs.waves, beam_m=beam, gravity=GRAVITY)
    physical = inputs.to_json()
    output.with_name(output.stem+'-inputs.json').write_text(physical, encoding='utf-8', newline='\n')
    report = dict(complete=False, qualified=False, format_version=1,
        scope='Coupled nonlinear six-DOF and free-surface state in a 12x6x2 m diagnostic tank; no open-water or temporal/spatial qualification',
        scipy_version=scipy.__version__, requested_duration_seconds=duration,
        input_sha256=sha256(physical.encode()).hexdigest(), methods=METHODS,
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
        discretization=dict(body_refinement=refinement, quadrature_order=quadrature_order),
        rhs=[], geometry=[],
        controls=dict(max_rhs=max_rhs, soft_elapsed_budget_seconds=budget_seconds,
            dense_output='Saved SciPy accepted-step RkDenseOutput coefficients',
            rtol=RELATIVE_TOLERANCE*tolerance_scale, tolerance_scale=tolerance_scale,
            first_step_seconds=min(max_step, duration), max_step_seconds=max_step,
            topology_templates_retained=4),
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in SOURCES},
        physical_source_sha256={name:sha256(Path(name).read_bytes()).hexdigest() for name in
            ('src/topwater/constants.py', 'src/topwater/inputs.py', 'src/topwater/sea_states.py',
             'src/topwater/frames.py', 'src/topwater/rankine.py', 'src/topwater/_validation.py',
             'src/topwater/surface.py', 'tests/laplace_reference.py')})
    started = perf_counter()
    times, states, dense_q, solver = [], [], [], None
    templates, cold_count, saved = [], 0, {}
    body = inputs.body

    def elapsed_check():
        if perf_counter()-started > budget_seconds:
            raise RuntimeError('bounded moving-case budget exhausted between operations')

    def cold(**state):
        nonlocal cold_count
        elapsed_check()
        path = output.with_name(output.stem+f'-geometry-{cold_count}.json')
        cold_count += 1
        before = perf_counter()
        _, arrays, metadata = construct(path, body=body, waves=inputs.waves,
            export_template=True, **state)
        report['geometry'].append(dict(path=str(path), seconds=perf_counter()-before,
            sha256=sha256(path.with_suffix('.npz').read_bytes()).hexdigest()))
        templates.insert(0, arrays)
        del templates[4:]
        return arrays, metadata

    def geometry(pose, eta, psi):
        elapsed_check()
        before = perf_counter()
        reasons = []
        for index, template in enumerate(templates):
            try:
                arrays, diagnostics = evaluate_template(template, pose=pose,
                    parent_eta=eta, parent_psi=psi)
            except InvalidTemplate as error:
                reasons.append(str(error))
                continue
            templates.insert(0, templates.pop(index))
            mode = 'current guarded chart'
            break
        else:
            _, metadata = cold(pose=pose,
                parent_xy=initial_geometry['parent_nodes'][:, :2], parents=initial_geometry['parents'],
                parent_eta=eta, parent_psi=psi)
            # One floating replay owns both the primal operator and tangent.
            arrays, diagnostics = evaluate_template(templates[0], pose=pose,
                parent_eta=eta, parent_psi=psi)
            diagnostics['reference_seconds'] = metadata['seconds']
            mode = 'exact chart refresh with owned replay'
        chart = templates[0]
        for _ in range(refinement):
            arrays = refine_boundary(arrays, region='body')
        row['geometry'].append(dict(mode=mode, seconds=perf_counter()-before,
            rejected_charts=reasons, diagnostics=diagnostics, faces=len(arrays['faces'])))
        return arrays, chart

    def rhs(time, state):
        nonlocal row
        elapsed_check()
        if len(report['rhs']) >= max_rhs:
            raise RuntimeError('bounded moving-case RHS count exhausted')
        before = perf_counter()
        row = dict(time=float(time), complete=False, geometry=[], fields=[])
        report['rhs'].append(row)
        print(json.dumps(dict(event='rhs-start', index=len(report['rhs']), time=float(time))), flush=True)
        if not np.isfinite(state).all() or np.linalg.norm(state[3:7]) == 0:
            raise ValueError('nonfinite state or zero quaternion')
        pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
        twist, eta, psi = state[7:13], state[13:13+size], state[13+size:]
        saved['attempted_rhs_state'], saved['attempted_rhs_time'] = state.copy(), float(time)
        arrays, chart = geometry(pose, eta, psi)
        before_field = perf_counter()
        assembly = prepare_energy(arrays['faces'], TrianglePanels(arrays['points'][arrays['faces']], 1),
                                  order=quadrature_order)
        assembled = perf_counter()
        current = current_field(arrays, assembly.energy, values, twist, acceleration_columns=True, time=time)
        residual = current['relative_field_residual']
        row['fields'].append(dict(energy_seconds=assembled-before_field,
            solve_seconds=perf_counter()-assembled, relative_residual=residual))
        if not np.isfinite(residual) or residual > FIELD_RESIDUAL_LIMIT:
            raise ValueError('current field exceeds the audited algebraic residual budget')
        def geometry_tangent():
            coarse, derivative = chart_tangent(chart, pose=pose, parent_eta=eta, parent_psi=psi)
            return refine_tangent(coarse, derivative['point_jacobian'], derivative['parent_map_jacobian'], refinement)
        result = action_rhs(body, arrays, current, assembly, state, values, time, geometry_tangent)
        raw_quaternion = state[3:7]
        quaternion_rate = .5*np.r_[-twist[3:]@raw_quaternion[1:],
            raw_quaternion[0]*twist[3:]+np.cross(twist[3:], raw_quaternion[1:])]
        rate = np.r_[twist[:3], quaternion_rate, result.acceleration, result.eta_rate, result.psi_rate]
        if not np.isfinite(rate).all():
            raise FloatingPointError('nonfinite coupled boat state rate')
        row.update(complete=True, seconds=perf_counter()-before, surface=result.surface,
            acceleration_world=result.acceleration.tolist(),
            hydrodynamic_load_world=result.hydrodynamic_load.tolist(), pressure=result.diagnostics,
            pressure_derivative=result.ale)
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
        initial_geometry, metadata = cold(pose=inputs.initial_state)
        values = coefficients(inputs.waves)
        size = len(initial_geometry['parent_nodes'])
        initial = np.r_[initial_geometry['pose_position'], initial_geometry['pose_orientation'],
            inputs.initial_state.linear_velocity,
            rotation_matrix(inputs.initial_state.orientation)@inputs.initial_state.angular_velocity,
            initial_geometry['parent_nodes'][:, 2], initial_geometry['parent_psi']]
        times, states = [0.], [initial.copy()]
        # Exploratory physical-unit integration controls, not acceptance limits.
        # They must be refined before interpreting a trajectory quantitatively.
        atol = np.r_[BODY_ABSOLUTE_TOLERANCE, np.full(size, SURFACE_HEIGHT_ABSOLUTE_TOLERANCE_M),
            np.full(size, SURFACE_POTENTIAL_ABSOLUTE_TOLERANCE_M2_PER_S)]*tolerance_scale
        report['controls']['absolute_tolerance_by_component'] = atol.tolist()
        solver = RK23(rhs, 0., initial, duration, first_step=min(max_step, duration),
            max_step=max_step, rtol=RELATIVE_TOLERANCE*tolerance_scale, atol=atol)
        report['benchmark_contract'] = dict(physics=report['scope'], input_sha256=report['input_sha256'],
            constants=dict(gravity_m_per_s2=GRAVITY, water_density_kg_per_m3=WATER_DENSITY),
            state_layout=report['state_layout'], state_count=len(initial), error_norm=ERROR_NORM,
            duration_seconds=duration, integrator='RK23',
            solve_controls=dict(method='current Cholesky/direct action coupling', true_field_residual_limit=FIELD_RESIDUAL_LIMIT),
            tolerances=dict(rtol=float(solver.rtol), absolute_by_component=solver.atol.tolist()),
            max_step_seconds=max_step, first_step_seconds=min(max_step, duration),
            spatial=dict(**report['discretization'], initial_discrete_geometry_sha256=array_identity(initial_geometry)))
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
        report['operator_counts'] = dict(rhs=len(report['rhs']),
            assemblies=sum(len(row['fields']) for row in report['rhs']))
        report['achieved_duration_seconds'] = report['accepted_simulated_seconds']
        report['benchmark_seconds'] = report['seconds']
        report['measurement_environment'] = measurement_environment(cp, resource_status, resource_context,
            'Cold geometry exports; RHS reports and accepted-step NPZ checkpoints in measured loop; final output excluded')
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
    parser.add_argument('--inputs', required=True, type=Path)
    parser.add_argument('--output', required=True)
    parser.add_argument('--duration', type=float, default=.02)
    parser.add_argument('--refinement', type=int, default=1)
    parser.add_argument('--max-rhs', type=int, default=13)
    parser.add_argument('--budget-seconds', type=float, default=90.)
    parser.add_argument('--max-step', type=float, default=.01)
    parser.add_argument('--quadrature-order', type=int, default=3, choices=(3, 6, 12))
    parser.add_argument('--tolerance-scale', type=float, default=1.)
    parser.add_argument('--resource-status', choices=['uncontrolled', 'contended', 'uncontended'], default='uncontrolled')
    parser.add_argument('--resource-context', default='Uncharacterized shared-machine load')
    options = vars(parser.parse_args())
    options['inputs'] = SimulationInputs.from_json(options['inputs'].read_text())
    result = run(**options)
    print(json.dumps({key: result[key] for key in
        ('complete', 'seconds', 'accepted_simulated_seconds', 'accepted_steps', 'cold_geometry_count')}))
