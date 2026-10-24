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
from topwater.surface import surface_integrals
from .boat_contact_reaction import load_snapshot, refine_snapshot
from .boat_geometry_template import evaluate_template, InvalidTemplate
from .boat_world_field import current_field, surface_rates
from .boat_parent_action import parent_mass
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(audit, output, order=6, directional=False, vertex=False, vertex_batch_width=None):
    from .cuda_operator import prepare_energy
    if order not in (3, 6, 12, 24):
        raise ValueError('bounded action control uses order three, six, twelve or twenty-four')
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
        directional_check=directional,
        vertex_gradient_check=vertex,
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_surface_action_audit.py', 'boat_contact_reaction.py', 'boat_geometry_template.py',
                         'boat_world_field.py', 'boat_parent_action.py', 'boat_world_boundary.py', 'canoe.py', 'incident_baseline_benchmark.py',
                         'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh', 'cuda_directional.py',
                         'cuda_directional.cu', 'cuda_dual.cuh', 'cuda_wide_dual.cuh', 'cuda_shape_gradient.py', 'cuda_panel_geometry.py', 'cuda_shape.cu')},
        directions=[], operations=[])
    saved, started = dict(state=state), perf_counter()
    output.parent.mkdir(parents=True, exist_ok=True)

    def geometry(height, current_pose=pose):
        current, _ = evaluate_template(chart, pose=current_pose, parent_eta=height, parent_psi=psi)
        return refine_snapshot(current, source, pilot)

    def evaluate(current, label, *, retain=False):
        if perf_counter()-started > 90:
            raise RuntimeError('90 s action audit budget exhausted')
        before = perf_counter()
        assembly = prepare_energy(current['faces'], TrianglePanels(current['points'][current['faces']], 1), order=order)
        energy = assembly.energy
        field = current_field(current, energy, values, twist, time=time, acceleration_columns=False)
        potential = field['potential']
        functional = WATER_DENSITY*(potential@energy@potential/2-potential@field['right'])
        volume, center = surface_integrals(current['points'], current['faces'])
        gravity_energy = WATER_DENSITY*GRAVITY*volume*center[2]
        report['operations'].append(dict(label=label, seconds=perf_counter()-before,
            functional_joules=float(functional), gravity_joules=float(gravity_energy),
            field_relative_residual=field['relative_field_residual']))
        return functional, gravity_energy, parent_mass(current), field, assembly if retain else None

    try:
        _, _, mass, base_field, base_assembly = evaluate(arrays, 'base', retain=directional or vertex)
        if vertex:
            from .cuda_shape_gradient import vertex_gradient, DEFAULT_BATCH_WIDTH
            width = DEFAULT_BATCH_WIDTH if vertex_batch_width is None else vertex_batch_width
            before = perf_counter()
            gradient = vertex_gradient(arrays['points'], arrays['faces'], base_field['potential'],
                assembly=base_assembly, density=WATER_DENSITY, order=order, batch_width=width)
            report['vertex_gradient'] = dict(seconds=perf_counter()-before,
                quadrature_order=order, batch_width=width, retained_current_factor=True)
            saved['energy_vertex_gradient'] = gradient
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
                derivative_row = dict(step_m=h, functional_derivative_newtons=float(df),
                    gravity_derivative_newtons=float(du),
                    independently_integrated_gravity_derivative_newtons=float(WATER_DENSITY*GRAVITY*(mass@eta)[index]),
                    mass_transport_newtons=float(transport), variational_weak_load_newtons=float(variational),
                    discrepancy_newtons=float(variational-current_weak_load[index]))
                if directional or vertex:
                    positive, negative = moved[2*level:2*level+2]
                    if any(not np.array_equal(current[key], arrays[key])
                           for current in (positive, negative) for key in ('faces', 'free_nodes')):
                        raise ValueError('shape derivative requires identical same-chart numbering')
                    moving = (positive['points']-negative['points'])/(2*h)
                    db = (plus[3]['right']-minus[3]['right'])/(2*h)
                    load_term = -WATER_DENSITY*(base_field['potential']@db)
                    fixed = arrays['free_nodes']
                    trace_change = (plus[3]['potential'][fixed]-minus[3]['potential'][fixed])/(2*h)
                    trace_term = WATER_DENSITY*(base_field['flux']@trace_change)
                    saved[f'eta{index}_h{level}_moving'] = moving
                    def chain(dk):
                        analytic = dk+load_term+trace_term
                        return dict(fixed_trace_energy_derivative_newtons=float(dk),
                            boundary_load_term_newtons=float(load_term), moving_trace_term_newtons=float(trace_term),
                            complete_functional_derivative_newtons=float(analytic),
                            scalar_finite_difference_error_newtons=float(analytic-df))
                    if directional:
                        from .cuda_directional import energy_direction_action
                        before = perf_counter()
                        action, diagnostics = energy_direction_action(arrays['points'], arrays['faces'],
                            moving, base_field['potential'], assembly=base_assembly)
                        derivative_row['directional'] = dict(seconds=perf_counter()-before,
                            **chain(WATER_DENSITY*(base_field['potential']@action)/2), diagnostics=diagnostics)
                        saved[f'eta{index}_h{level}_energy_action'] = action
                    if vertex:
                        derivative_row['vertex_gradient'] = chain(np.sum(gradient*moving))
                row['derivatives'].append(derivative_row)
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
    parser.add_argument('--directional', action='store_true')
    parser.add_argument('--vertex', action='store_true')
    parser.add_argument('--vertex-batch-width', type=int, choices=(1, 3, 6))
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))
