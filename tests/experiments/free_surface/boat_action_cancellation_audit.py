"""Matched reference/combined incident-port derivatives at one owned boat state."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from topwater import sea_states
from .boat_body_action import body_derivatives
from .boat_contact_reaction import load_snapshot, refine_snapshot
from .boat_geometry_template import evaluate_template
from .boat_parent_action_reference import action_rates
from .boat_parent_action import port_transport_covector
from .boat_world_field import current_field
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(reference, body_reference, output):
    reference, body_reference, output = Path(reference), Path(body_reference), fresh(output)
    control, body_control = (json.loads(p.read_text()) for p in (reference, body_reference))
    if not control['complete'] or not body_control['complete']:
        raise ValueError('complete parent and body controls required')
    for path, digest in control['input_sha256'].items():
        if sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('reference input changed: '+path)
    if body_control['input_sha256'][str(reference)] != sha256(reference.read_bytes()).hexdigest():
        raise ValueError('body control belongs to a different parent audit')
    for name, digest in control['source_sha256'].items():
        if name in ('boat_world_field.py', 'boat_parent_action.py'):
            continue  # These are the explicitly compared changes below.
        if sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() != digest:
            raise ValueError('unmodified physical source changed: '+name)
    for path in (reference, body_reference):
        original = json.loads(path.read_text())
        for name, digest in original['source_sha256'].items():
            archived = path.with_name(path.stem+'-source')/name
            if sha256(archived.read_bytes()).hexdigest() != digest:
                raise ValueError('original producer archive changed: '+name)
    vertex = json.loads(Path(next(p for p in control['input_sha256'] if p.endswith('.json') and '-geometry-' not in p)).read_text())
    snapshot = Path(next(p for p in vertex['input_sha256'] if p.endswith('.json') and '-geometry-' not in p))
    source, data, arrays, pilot, chart = load_snapshot(snapshot)
    with np.load(reference.with_suffix('.npz'), allow_pickle=False) as saved:
        energy, gradient, state = (saved[key] for key in ('energy', 'gradient', 'state'))
        owned_phi, owned_right, owned_rates = saved['potential'], saved['right'], saved['psi_rate']
        owned_load = saved['corrected_action_weak']
    np.testing.assert_array_equal(state, data['state'])
    metadata = json.loads(Path(source['chart_identity']['path']).with_suffix('.json').read_text())
    values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    time, twist, psi = source['time_seconds'], state[7:13], arrays['parent_psi']
    def geometry(pose, eta):
        return refine_snapshot(evaluate_template(chart, pose=pose, parent_eta=eta, parent_psi=psi)[0], source, pilot)
    field = current_field(arrays, energy, values, twist, time=time)
    np.testing.assert_array_equal(field['potential'], owned_phi)
    np.testing.assert_array_equal(field['right'], owned_right)
    report = dict(complete=False, qualified=False,
        scope='Exact incident-load/port cancellation parity and matched component timing; no physical or trajectory qualification.',
        input_sha256={str(p): sha256(p.read_bytes()).hexdigest() for p in (reference, reference.with_suffix('.npz'), body_reference)},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('boat_action_cancellation_audit.py', 'boat_parent_action.py', 'boat_parent_action_reference.py', 'boat_body_action.py', 'boat_world_field.py',
             'boat_geometry_template.py', 'boat_contact_reaction.py', 'boat_world_boundary.py', 'canoe.py', 'incident_baseline_benchmark.py')},
        original_source_sha256=control['source_sha256'], original_body_source_sha256=body_control['source_sha256'],
        calls=[])
    output.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    reference_parent, reference_body = None, None
    try:
        for combined in (False, True, True, False):
            before = perf_counter()
            covector = port_transport_covector(arrays, field['potential'], values, time) if combined else None
            parent = action_rates(arrays, field, gradient, values, twist, time, geometry,
                combine_incident_ports=combined, port_covector=covector)
            body = body_derivatives(arrays, field, gradient, values, twist, time, parent['eta_rate'], geometry,
                combine_incident_ports=combined, port_covector=covector)
            seconds = perf_counter()-before
            np.testing.assert_allclose(parent['psi_rate'], owned_rates, rtol=1e-11, atol=1e-10)
            np.testing.assert_allclose(parent['corrected_action_weak'], owned_load, rtol=1e-11, atol=1e-7)
            expected_body = sum(np.array(term) for term in body_control['pose_terms'].values())
            np.testing.assert_allclose(body['pose_load'], expected_body, rtol=1e-10, atol=1e-6)
            if not combined:
                reference_parent, reference_body = parent, body
            parent_error = np.max(abs(parent['corrected_action_weak']-reference_parent['corrected_action_weak']))
            body_error = np.max(abs(body['pose_load']-reference_body['pose_load']))
            row = dict(combined=combined, seconds=seconds, maximum_parent_weak_load_difference=float(parent_error),
                maximum_parent_rate_difference=float(np.max(abs(parent['psi_rate']-reference_parent['psi_rate']))),
                maximum_body_pose_load_difference=float(body_error))
            if combined:
                row['maximum_parent_incident_port_term_difference'] = float(np.max(abs(parent['terms']['incident_port']-
                    reference_parent['terms']['incident_load']-reference_parent['terms']['port_work'])))
                row['maximum_body_boundary_plus_port_term_difference'] = float(np.max(abs(body['terms']['boundary_load']+
                    body['terms']['incident_port']-reference_body['terms']['boundary_load']-reference_body['terms']['port_work'])))
            report['calls'].append(row)
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
