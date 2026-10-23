"""Compare candidate action and Bernoulli body loads on identical frozen fields."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.spatial.transform import Rotation

from topwater import InitialState, sea_states
from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from tests.reference_boat import GRAVITY, WATER_DENSITY, reference_body
from .boat_ale_pressure import differentiate
from .boat_body_action import body_derivatives, known_load
from .boat_contact_reaction import load_snapshot, refine_snapshot
from .boat_geometry_template import evaluate_template, InvalidTemplate
from .boat_world_field import current_field
from .body_pressure import evaluate_body_pressure
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(reference, output, scalar_poses=()):
    from .cuda_operator import prepare_energy
    reference, output = Path(reference), fresh(output)
    control = json.loads(reference.read_text())
    if not control['complete']:
        raise ValueError('complete parent-action audit required')
    for name, expected in control['source_sha256'].items():
        if sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() != expected:
            raise ValueError('current helper differs from parent audit: '+name)
    for path, expected in control['input_sha256'].items():
        if sha256(Path(path).read_bytes()).hexdigest() != expected:
            raise ValueError('parent audit input changed: '+path)
    vertex = Path(next(path for path in control['input_sha256'] if path.endswith('.json') and '-geometry-' not in path))
    vertex_control = json.loads(vertex.read_text())
    snapshot = Path(next(path for path in vertex_control['input_sha256'] if path.endswith('.json') and '-geometry-' not in path))
    source, data, arrays, pilot, chart = load_snapshot(snapshot)
    with np.load(reference.with_suffix('.npz'), allow_pickle=False) as saved:
        state, energy, gradient = (saved[key] for key in ('state', 'energy', 'gradient'))
        eta_rate, psi_rate = saved['eta_rate'], saved['psi_rate'][1]
        np.testing.assert_array_equal(state, data['state'])
        owned_potential, owned_right = saved['potential'], saved['right']
        gauge_gradient = saved['gauge_gradient']
    metadata_path = Path(source['chart_identity']['path']).with_suffix('.json')
    metadata = json.loads(metadata_path.read_text())
    values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    time, order = source['time_seconds'], control['order']
    count, twist = len(eta_rate), state[7:13]
    eta, psi = state[13:13+count], state[13+count:]
    pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
    body = reference_body()
    rotation = rotation_matrix(pose.orientation)
    inertia = rotation@body.inertia_cg@rotation.T
    def geometry(current_pose, height, potential=psi):
        return refine_snapshot(evaluate_template(chart, pose=current_pose, parent_eta=height, parent_psi=potential)[0], source, pilot)
    report = dict(complete=False, qualified=False, order=order,
        scope='Frozen-state body action versus Bernoulli pressure on identical candidate parent rates and field stencils; no trajectory qualification.',
        input_sha256={str(path): sha256(path.read_bytes()).hexdigest() for path in (reference, reference.with_suffix('.npz'), vertex, snapshot, metadata_path)},
        source_sha256={path.name: sha256(path.read_bytes()).hexdigest() for path in
            [Path(__file__), Path(__file__).with_name('boat_body_action.py'), Path(__file__).with_name('boat_ale_pressure.py'), Path(__file__).with_name('body_pressure.py')]},
        parent_source_sha256=control['source_sha256'])
    output.parent.mkdir(parents=True, exist_ok=True)
    saved, started = {}, perf_counter()
    try:
        field = current_field(arrays, energy, values, twist, time=time)
        np.testing.assert_array_equal(field['potential'], owned_potential)
        np.testing.assert_array_equal(field['right'], owned_right)
        derivatives = body_derivatives(arrays, field, gradient, values, twist, time, eta_rate, geometry)
        gauge_arrays = dict(arrays, parent_psi=psi+1.)
        gauge_field = current_field(gauge_arrays, energy, values, twist, time=time)
        gauge_derivatives = body_derivatives(gauge_arrays, gauge_field, gauge_gradient, values, twist, time,
            eta_rate, lambda current_pose, height: geometry(current_pose, height, psi+1.))
        report['body_geometry_evaluations'] = derivatives['geometry_evaluations']
        report['pose_steps'] = derivatives['steps'].tolist()
        report['pose_step_units'] = ['m']*3+['rad']*3
        report['pose_terms'] = {name: value.tolist() for name, value in derivatives['terms'].items()}
        report['scalar_pose_controls'] = []
        for index in scalar_poses:
            if index not in range(6):
                raise ValueError('scalar pose index must be 0..5')
            step = derivatives['steps'][index]
            functionals = []
            for offset in (step, -step, step/2, -step/2):
                direction = offset*np.eye(3)[index % 3]
                moved_pose = (InitialState(pose.position+direction, pose.orientation) if index < 3 else
                    InitialState(pose.position, (Rotation.from_rotvec(direction)*Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)))
                moved = geometry(moved_pose, eta)
                assembly = prepare_energy(moved['faces'], TrianglePanels(moved['points'][moved['faces']], 1), order=order)
                moved_field = current_field(moved, assembly.energy, values, twist, time=time, acceleration_columns=False)
                phi = moved_field['potential']
                functionals.append(float(WATER_DENSITY*(phi@assembly.energy@phi/2-phi@moved_field['right'])))
            scalar = np.array([(functionals[0]-functionals[1])/(2*step), (functionals[2]-functionals[3])/step])
            chain = -sum(derivatives['terms'][key][:, index] for key in ('kinetic', 'boundary_load', 'moving_trace'))
            report['scalar_pose_controls'].append(dict(index=index, functional_values_joules=functionals,
                scalar_derivative=scalar.tolist(), chain_derivative=chain.tolist(), difference=(chain-scalar).tolist()))
        hydrostatic = evaluate_body_pressure(arrays['points'], arrays['faces'], arrays['body_faces'], pose.position,
            body.mass, inertia, np.zeros(6), np.zeros(len(energy)), np.zeros(len(energy)), np.zeros_like(arrays['points']),
            np.zeros((len(energy), 6)), density=WATER_DENSITY, gravity=GRAVITY)
        report['hydrostatic_control'] = dict(action=derivatives['terms']['gravity'].tolist(),
            pressure=hydrostatic.known_pressure_load.tolist(),
            difference=(derivatives['terms']['gravity']-hydrostatic.known_pressure_load).tolist())
        for halving in range(13):
            h = 1e-5/2**halving
            offsets = {}
            try:
                for dt in (h, -h, h/2, -h/2):
                    q = (Rotation.from_rotvec(dt*twist[3:])*Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
                    moved = geometry(InitialState(pose.position+dt*twist[:3], q), eta+dt*eta_rate, psi+dt*psi_rate)
                    for key in ('faces', 'free_nodes', 'free_faces', 'body_faces', 'wall_faces', 'bottom_faces', 'node_chart_key'):
                        if not np.array_equal(moved[key], arrays[key]):
                            raise ValueError('body time derivative leaves owned chart')
                    offsets[dt] = moved
            except InvalidTemplate:
                continue
            break
        else:
            raise ValueError('no common body time stencil')
        fields, gauge_fields = {}, {}
        report['offset_assemblies'] = []
        for dt, moved in offsets.items():
            before = perf_counter()
            assembly = prepare_energy(moved['faces'], TrianglePanels(moved['points'][moved['faces']], 1), order=order)
            fields[dt] = current_field(moved, assembly.energy, values, twist, time=time+dt, acceleration_columns=False)
            gauge_fields[dt] = current_field(dict(moved, parent_psi=moved['parent_psi']+1.), assembly.energy,
                values, twist, time=time+dt, acceleration_columns=False)
            report['offset_assemblies'].append(dict(offset=dt, seconds=perf_counter()-before))
            saved[f'potential_{dt}'] = fields[dt]['potential']
        added_mass = WATER_DENSITY*field['body_load'].T@field['acceleration_potentials']
        report['added_mass_action'] = added_mass.tolist()
        report['added_mass_symmetry_relative'] = float(np.linalg.norm(added_mass-added_mass.T)/np.linalg.norm(added_mass))
        rigid = np.zeros((6, 6)); rigid[:3, :3] = body.mass*np.eye(3); rigid[3:, 3:] = inertia
        weight_gyro = np.r_[0., 0., -body.mass*GRAVITY, -np.cross(twist[3:], inertia@twist[3:])]
        report['comparisons'] = []
        for level, step in enumerate((h, h/2)):
            rate, mesh_velocity, diagnostics = differentiate(arrays, offsets[step], offsets[-step], fields[step]['potential'], fields[-step]['potential'], step, pose.position, twist)
            pressure = evaluate_body_pressure(arrays['points'], arrays['faces'], arrays['body_faces'], pose.position,
                body.mass, inertia, twist, field['potential'], rate, mesh_velocity, field['acceleration_potentials'],
                density=WATER_DENSITY, gravity=GRAVITY)
            action, transport, angular = known_load(derivatives, field, twist, fields[step], fields[-step], step, level)
            # Both routes contain the identical -rho B^T phi_dot. Remove it
            # before interpreting a dynamic mismatch as a time-rate problem.
            geometric_pressure = evaluate_body_pressure(arrays['points'], arrays['faces'], arrays['body_faces'], pose.position,
                body.mass, inertia, twist, field['potential'], np.zeros_like(rate), mesh_velocity,
                field['acceleration_potentials'], density=WATER_DENSITY, gravity=GRAVITY)
            phi_dot = (fields[step]['potential']-fields[-step]['potential'])/(2*step)
            common = -WATER_DENSITY*field['body_load'].T@phi_dot
            b_dot = (fields[step]['body_load']-fields[-step]['body_load'])/(2*step)
            geometric_action = derivatives['pose_load'][level]+angular-WATER_DENSITY*b_dot.T@field['potential']
            acceleration = np.linalg.solve(rigid+added_mass, action+weight_gyro)
            report['comparisons'].append(dict(step=step, action_known_load=action.tolist(), pressure_known_load=pressure.known_pressure_load.tolist(),
                load_difference=(action-pressure.known_pressure_load).tolist(), momentum_time_derivative=transport.tolist(), angular_term=angular.tolist(),
                added_mass_difference=(added_mass-pressure.added_mass).tolist(), action_acceleration=acceleration.tolist(), pressure_acceleration=pressure.acceleration.tolist(),
                acceleration_difference=(acceleration-pressure.acceleration).tolist(), pressure_diagnostics=pressure.diagnostics, ale_diagnostics=diagnostics))
            report['comparisons'][-1]['geometric_identity'] = dict(
                common_potential_rate_load=common.tolist(),
                pressure_temporal_identity_error=(pressure.known_pressure_load-geometric_pressure.known_pressure_load-common).tolist(),
                momentum_product_remainder=(transport-common+WATER_DENSITY*b_dot.T@field['potential']).tolist(),
                geometric_action_load=geometric_action.tolist(), geometric_pressure_load=geometric_pressure.known_pressure_load.tolist(),
                geometric_load_difference=(geometric_action-geometric_pressure.known_pressure_load).tolist())
            gauge_action, _, _ = known_load(gauge_derivatives, gauge_field, twist, gauge_fields[step], gauge_fields[-step], step, level)
            report['comparisons'][-1]['constant_gauge_action_load_change'] = (gauge_action-action).tolist()
        saved.update(state=state, energy=energy, gradient=gradient, eta_rate=eta_rate, psi_rate=psi_rate)
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['seconds'] = perf_counter()-started
        np.savez(output.with_suffix('.npz'), **saved)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--scalar-poses', nargs='*', type=int, default=[])
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))
