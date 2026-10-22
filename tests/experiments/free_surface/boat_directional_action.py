"""Current energy-action derivative on the held-out nonlinear boat stencil."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.spatial.transform import Rotation

from topwater import InitialState, sea_states
from topwater.rankine import TrianglePanels
from .boat_geometry_template import evaluate_template
from .boat_world_field import current_field, surface_rates, refine_boundary
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(snapshot, output, order=12):
    from .cuda_operator import prepare_energy, cp
    from .cuda_directional import energy_direction_action
    if order not in (3, 6, 12):
        raise ValueError('bounded comparison supports orders 3, 6 and 12')
    snapshot, output = Path(snapshot), fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    previous = json.loads(snapshot.read_text())
    if not previous.get('complete') or previous['quadrature_order'] != 12:
        raise ValueError('requires the complete uniform-order-twelve prediction reference')
    folder = Path(__file__).parent
    for name, expected in previous['source_sha256'].items():
        if name != 'boat_operator_prediction.py' and sha256((folder/name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'stencil producer helper changed: {name}')
    with np.load(snapshot.with_suffix('.npz'), allow_pickle=False) as data:
        saved = dict(data)
    chart_path = Path(previous['chart']['path'])
    if sha256(chart_path.read_bytes()).hexdigest() != previous['chart']['sha256']:
        raise ValueError('saved geometry chart changed')
    with np.load(chart_path, allow_pickle=False) as data:
        chart = dict(data)
    metadata = json.loads(chart_path.with_suffix('.json').read_text())
    values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    y = saved['selected_states'][-1]
    size = (len(y)-13)//2
    pose, twist = InitialState(y[:3], y[3:7]/np.linalg.norm(y[3:7])), y[7:13]

    def geometry(position, orientation, eta, psi):
        arrays = evaluate_template(chart, pose=InitialState(position, orientation),
                                   parent_eta=eta, parent_psi=psi)[0]
        for _ in range(2):
            arrays = refine_boundary(arrays, region='body')
        return arrays

    target = geometry(pose.position, pose.orientation, y[13:13+size], y[13+size:])
    if not np.allclose(target['points'], saved['target_points'], rtol=0, atol=1e-12):
        raise ValueError('reconstructed target geometry differs')
    actual = current_field(target, saved['target_energy'], values, twist, time=.15)
    eta_rate, psi_rate, _ = surface_rates(target['points'], target['faces'][target['free_faces']],
        target['free_nodes'], actual['parent_map'], target['parent_psi'], actual['flux'])
    report = dict(complete=False, qualified=False, records=[], quadrature_order=order,
        scope='One actual current boat dE@phi validation/cost comparison; not differentiated pressure or a trajectory',
        input_sha256={str(p): sha256(p.read_bytes()).hexdigest()
                      for p in (snapshot, snapshot.with_suffix('.npz'))},
        source_sha256={name: sha256((folder/name).read_bytes()).hexdigest() for name in
            ('boat_directional_action.py', 'cuda_directional.py', 'cuda_directional.cu',
             'cuda_dual.cuh', 'cuda_pair.cuh', 'cuda_shape_gradient.py', 'cuda_shape.cu',
             'cuda_operator.py', 'boat_world_field.py', 'boat_geometry_template.py')})
    arrays = {}
    started = perf_counter()
    try:
        assembly = prepare_energy(target['faces'], TrianglePanels(target['points'][target['faces']], 1), order=order)
        if order == 12:
            report['primal_relative_replay_error'] = float(np.linalg.norm(assembly.energy-saved['target_energy'])
                                                          /np.linalg.norm(saved['target_energy']))
        for h in (.001, .0005):
            moved = []
            for dt in (h, -h):
                orientation = (Rotation.from_rotvec(dt*twist[3:])
                    *Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
                current = geometry(pose.position+dt*twist[:3], orientation,
                    y[13:13+size]+dt*eta_rate, y[13+size:]+dt*psi_rate)
                if not np.array_equal(current['faces'], target['faces']):
                    raise ValueError('direction requires one unchanged chart basis')
                moved.append(current)
            direction = (moved[0]['points']-moved[1]['points'])/(2*h)
            derivative, diagnostics = energy_direction_action(target['points'], target['faces'],
                direction, actual['potential'], assembly=assembly)
            seconds, finite_difference_seconds = [], []

            def action():
                before = perf_counter()
                result = energy_direction_action(target['points'], target['faces'],
                    direction, actual['potential'], assembly=assembly)
                seconds.append(perf_counter()-before)
                return result

            def finite_difference():
                before = perf_counter()
                matrices = [prepare_energy(a['faces'], TrianglePanels(a['points'][a['faces']], 1),
                                           order=order).energy for a in moved]
                result = ((matrices[0]-matrices[1])/(2*h))@actual['potential']
                finite_difference_seconds.append(perf_counter()-before)
                return result

            # Alternate method ordering to expose, rather than hide, changing
            # GPU clock/thermal conditions. Count complete current actions.
            for repetition in range(3):
                if repetition % 2:
                    reference = finite_difference()
                    derivative, diagnostics = action()
                else:
                    derivative, diagnostics = action()
                    reference = finite_difference()
            row = dict(step=h, action_norm=float(np.linalg.norm(derivative)),
                finite_difference_norm=float(np.linalg.norm(reference)),
                relative_action_error=float(np.linalg.norm(derivative-reference)/np.linalg.norm(reference)),
                maximum_absolute_action_error=float(np.max(abs(derivative-reference))),
                median_derivative_action_seconds=float(np.median(seconds)),
                derivative_action_seconds=seconds,
                median_two_assembly_difference_seconds=float(np.median(finite_difference_seconds)),
                two_assembly_difference_seconds=finite_difference_seconds, diagnostics=diagnostics)
            arrays[f'action_{h}'], arrays[f'reference_{h}'] = derivative, reference
            arrays[f'direction_{h}'] = direction
            report['records'].append(row)
            print(json.dumps(row), flush=True)
        # Current-primal reuse rejects geometric mismatches before any kernel.
        try:
            energy_direction_action(target['points']+1e-6, target['faces'], direction,
                                    actual['potential'], assembly=assembly)
        except ValueError:
            report['stale_geometry_rejected'] = True
        else:
            raise ValueError('stale geometry was accepted')
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['seconds'] = perf_counter()-started
        np.savez(output.with_suffix('.npz'), **arrays)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--order', type=int, choices=(3, 6, 12), default=12)
    run(**vars(parser.parse_args()))
