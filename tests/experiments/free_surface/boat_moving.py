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
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(output, duration=.02, refinement=1, max_rhs=13, budget_seconds=90., max_step=.01):
    from .cuda_operator import prepare_energy, cp
    if (not np.isfinite(duration) or duration <= 0 or refinement not in (0, 1, 2)
            or max_rhs < 4 or not np.isfinite(budget_seconds) or budget_seconds <= 0
            or not np.isfinite(max_step) or max_step <= 0):
        raise ValueError('positive bounded duration/budget, >=4 RHS and refinement 0..2 required')
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    report = dict(complete=False, qualified=False,
        scope='Coupled nonlinear six-DOF and free-surface state in a 12x6x2 m diagnostic tank; no open-water or temporal/spatial qualification',
        scipy_version=scipy.__version__, requested_duration_seconds=duration,
        state_layout=['position_world_m[3]', 'orientation_wxyz[4]',
            'velocity_world_m_per_s[3]', 'angular_velocity_world_rad_per_s[3]',
            'parent_eta_m[N]', 'parent_psi_m2_per_s[N]'],
        boundary_refinement=dict(region='body', levels=refinement), rhs=[], geometry=[],
        controls=dict(max_rhs=max_rhs, soft_elapsed_budget_seconds=budget_seconds,
            derivative_step_seconds=.001, rtol=.001, first_step_seconds=min(max_step, duration),
            max_step_seconds=max_step, topology_templates_retained=4),
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_moving.py', 'boat_geometry_template.py', 'boat_world_boundary.py',
                         'boat_parent_clip_proof.py', 'boat_world_field.py', 'body_pressure.py',
                         'canoe.py', 'incident_baseline_benchmark.py',
                         'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh')})
    started = perf_counter()
    times, states, solver = [], [], None
    templates, cold_count = [], 0
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

    def geometry(pose, eta, psi, label):
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

    def field(arrays, twist, time, acceleration_columns):
        elapsed_check()
        before = perf_counter()
        energy = prepare_energy(arrays['faces'], TrianglePanels(arrays['points'][arrays['faces']], 1)).energy
        assembled = perf_counter()
        value = current_field(arrays, energy, values, twist,
                              acceleration_columns=acceleration_columns, time=time)
        row['fields'].append(dict(energy_seconds=assembled-before,
            solve_seconds=perf_counter()-assembled, relative_residual=value['relative_field_residual']))
        if not np.isfinite(value['relative_field_residual']) or value['relative_field_residual'] > 1e-9:
            raise ValueError('current field exceeds the audited algebraic residual budget')
        return value

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
        quaternion = state[3:7]/np.linalg.norm(state[3:7])
        pose = InitialState(state[:3], quaternion)
        twist, eta, psi = state[7:13], state[13:13+size], state[13+size:]
        arrays = geometry(pose, eta, psi, 'base')
        current = field(arrays, twist, time, True)
        eta_rate, psi_rate, surface = surface_rates(arrays['points'],
            arrays['faces'][arrays['free_faces']], arrays['free_nodes'], current['parent_map'],
            psi, current['flux'])
        samples, groups = body_material_samples(arrays, quadrature()[0])
        hints = np.broadcast_to(np.arange(len(arrays['body_faces']))[:, None], groups.shape)
        traces = []
        for sign in (1, -1):
            offset = sign*.001
            orientation = (Rotation.from_rotvec(offset*twist[3:])
                *Rotation.from_quat(quaternion, scalar_first=True)).as_quat(scalar_first=True)
            moved = geometry(InitialState(pose.position+offset*twist[:3], orientation),
                eta+offset*eta_rate, psi+offset*psi_rate, f'offset-{sign}')
            moved_field = field(moved, twist, time+offset, False)
            traces.append(interpolate_body_trace(moved, moved_field['potential'], samples, groups,
                                                  face_hints=hints))
        rotation = rotation_matrix(quaternion)
        pressure = evaluate_body_pressure(arrays['points'], arrays['faces'], arrays['body_faces'],
            pose.position, body.mass, rotation@body.inertia_cg@rotation.T, twist, current['potential'],
            (traces[0]-traces[1])/.002,
            twist[:3]+np.cross(twist[3:], arrays['points']-pose.position),
            current['acceleration_potentials'], density=WATER_DENSITY, gravity=GRAVITY)
        raw_quaternion = state[3:7]
        quaternion_rate = .5*np.r_[-twist[3:]@raw_quaternion[1:],
            raw_quaternion[0]*twist[3:]+np.cross(twist[3:], raw_quaternion[1:])]
        rate = np.r_[twist[:3], quaternion_rate, pressure.acceleration, eta_rate, psi_rate]
        if not np.isfinite(rate).all():
            raise FloatingPointError('nonfinite coupled boat state rate')
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
        if states:
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
    result = run(**vars(parser.parse_args()))
    print(json.dumps({key: result[key] for key in
        ('complete', 'seconds', 'accepted_simulated_seconds', 'accepted_steps', 'cold_geometry_count')}))
