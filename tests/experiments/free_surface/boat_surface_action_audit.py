"""Two-direction frozen-state test of the restricted surface Euler equation."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.spatial.transform import Rotation

from topwater import InitialState, sea_states
from topwater.rankine import TrianglePanels
from tests.reference_boat import GRAVITY, WATER_DENSITY
from tests.reference_hull import surface_integrals
from .boat_contact_reaction import load_snapshot, refine_snapshot
from .boat_geometry_template import evaluate_template, InvalidTemplate
from .boat_world_field import current_field, surface_rates
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def parent_mass(arrays):
    """Integrate fixed world parent basis products over exposed projected faces."""
    mapping = np.zeros((len(arrays['points']), len(arrays['parent_nodes'])))
    for column in range(3):
        ids = arrays['free_parent_indices'][:, column]
        valid = ids >= 0
        mapping[arrays['free_nodes'][valid], ids[valid]] += arrays['free_parent_weights'][valid, column]
    faces = arrays['faces'][arrays['free_faces']]
    triangles = arrays['points'][faces, :2]
    area = np.linalg.det(triangles[:, 1:]-triangles[:, :1])/2
    if np.any(area <= 0):
        raise ValueError('invalid projected free triangle')
    trace = mapping[faces]
    unit = (np.ones((3, 3))+np.eye(3))/12
    return np.einsum('f,fci,cd,fdj->ij', area, trace, unit, trace, optimize=True)


def run(audit, output, order=6):
    from .cuda_operator import prepare_energy
    if order not in (6, 12, 24):
        raise ValueError('bounded action control uses order six, twelve or twenty-four')
    audit, output = Path(audit), fresh(output)
    source, data, arrays, pilot, chart = load_snapshot(audit)
    fixture = Path(__file__).resolve().parents[2]/'reference_boat.py'
    if sha256(fixture.read_bytes()).hexdigest() != pilot['fixture_sha256']:
        raise ValueError('current physical fixture differs from the saved pilot')
    state, time = data['state'], source['time_seconds']
    size = (len(state)-13)//2
    eta, psi, twist = state[13:13+size], state[13+size:], state[7:13]
    pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
    metadata_path = Path(source['chart_identity']['path']).with_suffix('.json')
    metadata = json.loads(metadata_path.read_text())
    values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    report = dict(complete=False, qualified=False, order=order,
        scope='Finite differences of complete stationary mixed functional and moving surface one-form in two parent directions; no trajectory or implementation of replacement rates.',
        input_sha256={str(p): sha256(p.read_bytes()).hexdigest()
            for p in (audit, audit.with_suffix('.npz'), metadata_path, fixture)},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_surface_action_audit.py', 'boat_contact_reaction.py', 'boat_geometry_template.py',
                         'boat_world_field.py', 'boat_world_boundary.py', 'canoe.py', 'incident_baseline_benchmark.py',
                         'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh')},
        directions=[], operations=[])
    saved, started = dict(state=state), perf_counter()
    output.parent.mkdir(parents=True, exist_ok=True)

    def geometry(height, current_pose=pose):
        current, _ = evaluate_template(chart, pose=current_pose, parent_eta=height, parent_psi=psi)
        return refine_snapshot(current, source, pilot)

    def evaluate(current, label):
        if perf_counter()-started > 90:
            raise RuntimeError('90 s action audit budget exhausted')
        before = perf_counter()
        energy = prepare_energy(current['faces'], TrianglePanels(current['points'][current['faces']], 1), order=order).energy
        field = current_field(current, energy, values, twist, time=time, acceleration_columns=False)
        potential = field['potential']
        functional = WATER_DENSITY*(potential@energy@potential/2-potential@field['right'])
        volume, center = surface_integrals(current['points'], current['faces'])
        gravity_energy = WATER_DENSITY*GRAVITY*volume*center[2]
        report['operations'].append(dict(label=label, seconds=perf_counter()-before,
            functional_joules=float(functional), gravity_joules=float(gravity_energy),
            field_relative_residual=field['relative_field_residual']))
        return functional, gravity_energy, parent_mass(current), field

    try:
        _, _, mass, base_field = evaluate(arrays, 'base')
        expected_mass = base_field['parent_map'].T@data['free_mass']@base_field['parent_map']
        np.testing.assert_allclose(mass, expected_mass, rtol=1e-13, atol=1e-13)
        eta_rate, psi_rate, surface = surface_rates(arrays['points'], arrays['faces'][arrays['free_faces']],
            arrays['free_nodes'], base_field['parent_map'], psi, base_field['flux'])
        current_weak_load = WATER_DENSITY*mass@psi_rate
        # Cdot includes body pose and eta motion, not psi evolution.
        for halving in range(13):
            dt = 1e-5/2**halving
            try:
                time_geometries = []
                for offset in (dt, -dt, dt/2, -dt/2):
                    orientation = (Rotation.from_rotvec(offset*twist[3:])
                        *Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
                    time_geometries.append(geometry(eta+offset*eta_rate,
                        InitialState(pose.position+offset*twist[:3], orientation)))
            except InvalidTemplate:
                continue
            break
        else:
            raise ValueError('no same-chart mass transport stencil')
        time_mass = [parent_mass(current) for current in time_geometries]
        cdot = [(time_mass[0]-time_mass[1])/(2*dt), (time_mass[2]-time_mass[3])/dt]
        report['mass_transport'] = dict(step_seconds=dt,
            h_h2_matrix_maximum_difference=float(np.max(abs(cdot[0]-cdot[1]))))
        outer_nodes = np.intersect1d(arrays['free_nodes'],
            np.unique(arrays['faces'][arrays['wall_faces']]))
        outer_rows = np.flatnonzero(np.isin(arrays['free_nodes'], outer_nodes))
        for index, label in ((7, 'central interior contact parent'), (10, 'bow interior contact parent')):
            if np.max(abs(base_field['parent_map'][outer_rows, index])) > 1e-14:
                raise ValueError('height variation must vanish on the external wall boundary')
            direction = np.eye(size)[index]
            for halving in range(13):
                step = 1e-5/2**halving
                try:
                    moved = [geometry(eta+offset*direction) for offset in (step, -step, step/2, -step/2)]
                except InvalidTemplate:
                    continue
                break
            else:
                raise ValueError('no same-chart parent-height stencil')
            evaluations = [evaluate(current, f'eta{index}:{offset:g}')
                for current, offset in zip(moved, (step, -step, step/2, -step/2), strict=True)]
            row = dict(parent_index=index, label=label, step_m=step,
                external_wall_trace_maximum=float(np.max(abs(base_field['parent_map'][outer_rows, index]))),
                current_bernoulli_weak_load_newtons=float(current_weak_load[index]), derivatives=[])
            for level, (plus, minus, h) in enumerate(((evaluations[0], evaluations[1], step),
                                                      (evaluations[2], evaluations[3], step/2))):
                df, du = (plus[0]-minus[0])/(2*h), (plus[1]-minus[1])/(2*h)
                dc = (plus[2]-minus[2])/(2*h)
                transport = WATER_DENSITY*(psi@dc@eta_rate-(cdot[level].T@psi)[index])
                variational = -df-du+transport
                row['derivatives'].append(dict(step_m=h, functional_derivative_newtons=float(df),
                    gravity_derivative_newtons=float(du),
                    independently_integrated_gravity_derivative_newtons=float(WATER_DENSITY*GRAVITY*(mass@eta)[index]),
                    mass_transport_newtons=float(transport), variational_weak_load_newtons=float(variational),
                    discrepancy_newtons=float(variational-current_weak_load[index])))
            report['directions'].append(row)
        saved.update(parent_mass=mass, eta_rate=eta_rate, psi_rate=psi_rate, mass_time_derivative=cdot[-1])
        report.update(complete=True, base_surface=surface)
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
    parser.add_argument('--audit', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--order', type=int, default=6)
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))
