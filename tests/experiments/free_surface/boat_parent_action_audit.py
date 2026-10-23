"""All-parent action loads, open-port work and gauge at one saved boat state."""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from topwater import sea_states
from topwater.rankine import TrianglePanels
from tests.reference_boat import WATER_DENSITY
from . import boat_world_field as world
from .boat_contact_reaction import load_snapshot, refine_snapshot
from .boat_geometry_template import evaluate_template
from .boat_parent_action import action_rates
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(reference, output):
    from .cuda_operator import prepare_energy
    from .cuda_shape_gradient import vertex_gradient
    reference, output = Path(reference), fresh(output)
    control = json.loads(reference.read_text())
    if not control['complete'] or not control.get('vertex_gradient_check'):
        raise ValueError('complete owned vertex/action control required')
    for path, digest in control['input_sha256'].items():
        if sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('control input changed: '+path)
    for name in ('cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh', 'canoe.py', 'incident_baseline_benchmark.py'):
        if sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() != control['source_sha256'][name]:
            raise ValueError('physical primal/load helper changed: '+name)
    snapshot = next(Path(path) for path in control['input_sha256'] if path.endswith('.json') and '-geometry-' not in path)
    source, data, arrays, pilot, chart = load_snapshot(snapshot)
    with np.load(reference.with_suffix('.npz'), allow_pickle=False) as archive:
        gradient = archive['energy_vertex_gradient']
        if not np.array_equal(archive['state'], data['state']):
            raise ValueError('control state differs')
    fixture = Path(__file__).resolve().parents[2]/'reference_boat.py'
    if sha256(fixture.read_bytes()).hexdigest() != pilot['fixture_sha256']:
        raise ValueError('physical fixture changed')
    metadata_path = Path(source['chart_identity']['path']).with_suffix('.json')
    metadata = json.loads(metadata_path.read_text())
    values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    state, time, order = data['state'], source['time_seconds'], control['order']
    size = (len(state)-13)//2
    psi, twist = state[13+size:], state[7:13]
    report = dict(complete=False, qualified=False, order=order,
        scope='Candidate all-parent restricted-action rates with open-port virtual work at one frozen state; no trajectory, body-force or physical accuracy qualification.',
        input_sha256={str(path): sha256(path.read_bytes()).hexdigest()
            for path in (reference, reference.with_suffix('.npz'), fixture, metadata_path)},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_parent_action_audit.py', 'boat_parent_action.py', 'boat_contact_reaction.py',
                         'boat_world_field.py', 'boat_world_boundary.py', 'boat_geometry_template.py', 'canoe.py',
                         'incident_baseline_benchmark.py', 'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh',
                         'cuda_shape_gradient.py', 'cuda_panel_geometry.py', 'cuda_shape.cu', 'cuda_dual.cuh', 'cuda_wide_dual.cuh')},
        geometry_replay_verification=source.get('geometry_replay_verification', {}),
        cached_gradient_scope='Reused ONLY at its identical frozen geometry/physical state; not an operator refresh/cache policy.')
    output.parent.mkdir(parents=True, exist_ok=True)
    saved, started = {}, perf_counter()
    try:
        assembly = prepare_energy(arrays['faces'], TrianglePanels(arrays['points'][arrays['faces']], 1), order=order)
        field = world.current_field(arrays, assembly.energy, values, twist, time=time)
        # This trusted, SHA-verified archived function is a numerical refactor
        # oracle in the audit ONLY. It cannot enter the candidate rate method.
        old_path = snapshot.with_name(snapshot.stem+'-source')/'boat_world_field.py'
        if sha256(old_path.read_bytes()).hexdigest() != control['source_sha256']['boat_world_field.py']:
            raise ValueError('original field oracle is not the control producer')
        old_node = next(node for node in ast.parse(old_path.read_text()).body
                        if isinstance(node, ast.FunctionDef) and node.name == 'current_field')
        namespace = dict(vars(world))
        exec(compile(ast.Module(body=[old_node], type_ignores=[]), str(old_path), 'exec'), namespace)
        old_field = namespace['current_field'](arrays, assembly.energy, values, twist, time=time)
        for name in field:
            np.testing.assert_array_equal(field[name], old_field[name], err_msg=name)
        report['boundary_refactor'] = dict(all_field_outputs_bitwise_equal=True,
            includes='Potential, reaction, physical RHS, trace map, body load and all six acceleration responses.',
            original_source_sha256=sha256(old_path.read_bytes()).hexdigest())
        functional = WATER_DENSITY*(field['potential']@assembly.energy@field['potential']/2-field['potential']@field['right'])
        np.testing.assert_allclose(functional, control['operations'][0]['functional_joules'], rtol=1e-12, atol=1e-8)
        def geometry(current_pose, height, potential_state=psi):
            current = evaluate_template(chart, pose=current_pose, parent_eta=height, parent_psi=potential_state)[0]
            return refine_snapshot(current, source, pilot)
        before = perf_counter()
        rates = action_rates(arrays, field, gradient, values, twist, time, geometry)
        report['geometry_load_rate_seconds'] = perf_counter()-before
        old_eta, old_psi, _ = world.surface_rates(arrays['points'], arrays['faces'][arrays['free_faces']],
            arrays['free_nodes'], field['parent_map'], psi, field['flux'])
        np.testing.assert_allclose(rates['eta_rate'], old_eta, rtol=1e-11, atol=1e-11)
        report['interior_chain_replay'] = []
        for row in control['directions']:
            index = row['parent_index']
            derivative = row['derivatives'][-1]
            # Compare complete functional derivatives before gravity/transport.
            df = -sum(rates['terms'][name][1, index] for name in ('kinetic', 'body_load', 'incident_load', 'moving_trace'))
            report['interior_chain_replay'].append(dict(parent_index=index,
                functional_derivative_difference_newtons=float(df-derivative['vertex_gradient']['complete_functional_derivative_newtons'])))
        gauge_arrays = dict(arrays, parent_psi=psi+1.)
        gauge_field = world.current_field(gauge_arrays, assembly.energy, values, twist, time=time)
        before = perf_counter()
        gauge_gradient = vertex_gradient(arrays['points'], arrays['faces'], gauge_field['potential'],
            assembly=assembly, order=order, density=WATER_DENSITY)
        gauge_rates = action_rates(gauge_arrays, gauge_field, gauge_gradient, values, twist, time,
            lambda current_pose, height: geometry(current_pose, height, psi+1.))
        report['gauge_check'] = dict(shift_m2_per_s=1., seconds=perf_counter()-before,
            maximum_eta_rate_change=float(np.max(abs(gauge_rates['eta_rate']-rates['eta_rate']))),
            maximum_corrected_psi_rate_change=float(np.max(abs(gauge_rates['psi_rate']-rates['psi_rate']))),
            maximum_raw_weak_load_change_newtons=float(np.max(abs(gauge_rates['raw_action_weak']-rates['raw_action_weak']))),
            maximum_corrected_weak_load_change_newtons=float(np.max(abs(gauge_rates['corrected_action_weak']-rates['corrected_action_weak']))))
        report['rates'] = dict(eta_rate=rates['eta_rate'].tolist(),
            psi_rate_h=rates['psi_rate'][0].tolist(), psi_rate_h2=rates['psi_rate'][1].tolist(),
            previous_bernoulli_psi_rate=old_psi.tolist(),
            maximum_h_h2_rate_difference=float(np.max(abs(rates['psi_rate'][0]-rates['psi_rate'][1]))),
            maximum_previous_rate_difference=float(np.max(abs(rates['psi_rate'][1]-old_psi))))
        report['port_work'] = dict(discrete_h2_newtons=rates['terms']['port_work'][1].tolist(),
            seam_order12_newtons=rates['port_seam_weak'].tolist(), seam_order6_newtons=rates['port_seam_weak_order6'].tolist(),
            maximum_discrete_seam_difference_newtons=float(np.max(abs(rates['terms']['port_work'][1]-rates['port_seam_weak']))),
            maximum_seam_order6_12_difference_newtons=float(np.max(abs(rates['port_seam_weak']-rates['port_seam_weak_order6']))))
        report['terms_h2_newtons'] = {name: value[1].tolist() for name, value in rates['terms'].items()}
        report['height_steps_m'] = rates['height_steps_m'].tolist()
        report['geometry_evaluations'] = rates['geometry_evaluations']
        for name, value in rates.items():
            if isinstance(value, np.ndarray): saved[name] = value
        saved.update(state=state, energy=assembly.energy, potential=field['potential'], right=field['right'],
                     gradient=gradient, gauge_gradient=gauge_gradient, gauge_psi_rate=gauge_rates['psi_rate'])
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
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))
