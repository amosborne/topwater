"""Analytic owned-chart action chain versus the saved nonlinear FD controls."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from topwater import InitialState, sea_states
from .boat_action_tangent import action_terms, refine_tangent
from .boat_chart_tangent import chart_tangent
from .boat_contact_reaction import load_snapshot
from .boat_world_field import current_field
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(reference, body_reference, output):
    reference, body_reference, output = Path(reference), Path(body_reference), fresh(output)
    parent_control, body_control = (json.loads(path.read_text()) for path in (reference, body_reference))
    if not parent_control['complete'] or not body_control['complete']:
        raise ValueError('complete parent/body FD controls required')
    for path, digest in parent_control['input_sha256'].items():
        if sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('control input changed: '+path)
    if body_control['input_sha256'][str(reference)] != sha256(reference.read_bytes()).hexdigest():
        raise ValueError('body control belongs to another parent control')
    for name in ('cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh', 'canoe.py', 'incident_baseline_benchmark.py'):
        if sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() != parent_control['source_sha256'][name]:
            raise ValueError('physical primal/load helper changed: '+name)
    for path, control in ((reference, parent_control), (body_reference, body_control)):
        for name, digest in control['source_sha256'].items():
            if sha256((path.with_name(path.stem+'-source')/name).read_bytes()).hexdigest() != digest:
                raise ValueError('original source archive changed: '+name)
    vertex = json.loads(Path(next(p for p in parent_control['input_sha256'] if p.endswith('.json') and '-geometry-' not in p)).read_text())
    snapshot = Path(next(p for p in vertex['input_sha256'] if p.endswith('.json') and '-geometry-' not in p))
    source, data, arrays, pilot, chart = load_snapshot(snapshot)
    with np.load(reference.with_suffix('.npz'), allow_pickle=False) as saved:
        state, energy, gradient, old_weak = (saved[key] for key in ('state', 'energy', 'gradient', 'corrected_action_weak'))
        old_rates, owned_phi, owned_right = saved['psi_rate'], saved['potential'], saved['right']
    np.testing.assert_array_equal(state, data['state'])
    metadata = json.loads(Path(source['chart_identity']['path']).with_suffix('.json').read_text())
    values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    time, twist = source['time_seconds'], state[7:13]
    field = current_field(arrays, energy, values, twist, time=time)
    np.testing.assert_array_equal(field['potential'], owned_phi)
    np.testing.assert_array_equal(field['right'], owned_right)
    report = dict(complete=False, qualified=False, order=parent_control['order'],
        scope='Exact current-chart geometry/load/action derivatives versus owned nonlinear h/h2 controls; no trajectory or continuum qualification.',
        input_sha256={str(path): sha256(path.read_bytes()).hexdigest() for path in (reference, reference.with_suffix('.npz'), body_reference)},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('boat_action_tangent_audit.py', 'boat_action_tangent.py', 'boat_chart_tangent.py', 'boat_parent_action.py',
             'boat_body_action.py', 'boat_world_field.py', 'boat_geometry_template.py', 'boat_contact_reaction.py',
             'boat_world_boundary.py', 'canoe.py', 'incident_baseline_benchmark.py')}, calls=[])
    output.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    try:
        for _ in range(4):
            before = perf_counter()
            coarse, derivative = chart_tangent(chart, pose=InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7])),
                parent_eta=arrays['parent_nodes'][:, 2], parent_psi=arrays['parent_psi'])
            fine, dx, ds = refine_tangent(coarse, derivative['point_jacobian'], derivative['parent_map_jacobian'],
                pilot['boundary_refinement']['levels'], local_body_refinement=source.get('local_body_refinement'))
            for key in ('points', 'faces', 'free_nodes', 'free_parent_indices', 'free_parent_weights'):
                np.testing.assert_array_equal(fine[key], arrays[key])
            chain = action_terms(fine, field, gradient, values, twist, time, dx, ds)
            seconds = perf_counter()-before
            body_terms = {key: np.array(value) for key, value in body_control['pose_terms'].items()}
            body_pose = sum(body_terms.values())
            parent_terms = {key: np.array(value) for key, value in parent_control['terms_h2_newtons'].items()}
            parent_terms['incident_port'] = parent_terms.pop('incident_load')+parent_terms.pop('port_work')
            body_terms['incident_port'] = body_terms.pop('port_work')
            # For body geometry the combined boundary term includes incident
            # cancellation. Compare the sum, not the separately defined terms.
            body_combined = body_terms.pop('boundary_load')+body_terms.pop('incident_port')
            report['calls'].append(dict(seconds=seconds,
                parent_rate_h_h2_errors=(chain['psi_rate'][None]-old_rates).tolist(),
                parent_weak_h_h2_errors=(chain['parent_weak_load'][None]-old_weak).tolist(),
                body_pose_h_h2_errors=(chain['body_pose_load'][None]-body_pose).tolist(),
                parent_term_h2_errors={key: (value-parent_terms[key]).tolist() for key, value in chain['parent_terms'].items()},
                body_term_h_h2_errors={key: (value-body_terms[key]).tolist() for key, value in chain['body_terms'].items()
                    if key not in ('boundary_load', 'incident_port')},
                body_boundary_plus_port_h_h2_errors=(chain['body_terms']['boundary_load'][None]+chain['body_terms']['incident_port'][None]-body_combined).tolist(),
                maximum_incident_derivative_residual=derivative['maximum_incident_derivative_residual']))
            np.testing.assert_allclose(chain['psi_rate'], old_rates[1], rtol=1e-6, atol=1e-5)
            np.testing.assert_allclose(chain['body_pose_load'], body_pose[1], rtol=1e-6, atol=.002)
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['seconds'] = perf_counter()-started
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--body-reference', required=True)
    parser.add_argument('--output', required=True)
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))
