"""Matched current-boat gradient timing; width batching changes arithmetic only."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from topwater import sea_states
from topwater.rankine import TrianglePanels
from tests.reference_boat import WATER_DENSITY
from .boat_contact_reaction import load_snapshot
from .boat_world_field import current_field
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(reference, output):
    from .cuda_operator import prepare_energy
    from .cuda_shape_gradient import vertex_gradient, _kernels
    reference, output = Path(reference), fresh(output)
    original = json.loads(reference.read_text())
    if not original['complete'] or not original.get('vertex_gradient_check'):
        raise ValueError('complete scalar/directional/vertex audit required')
    for path, digest in original['input_sha256'].items():
        if sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('reference input changed: '+path)
    for name in ('boat_world_field.py', 'boat_world_boundary.py', 'boat_geometry_template.py',
                 'canoe.py', 'incident_baseline_benchmark.py', 'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh'):
        if sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() != original['source_sha256'][name]:
            raise ValueError('physical helper changed: '+name)
    snapshot = next(Path(path) for path in original['input_sha256']
                    if path.endswith('.json') and '-geometry-' not in path)
    source, data, arrays, pilot, chart = load_snapshot(snapshot)
    fixture = Path(__file__).resolve().parents[2]/'reference_boat.py'
    if sha256(fixture.read_bytes()).hexdigest() != pilot['fixture_sha256']:
        raise ValueError('physical fixture changed')
    with np.load(reference.with_suffix('.npz'), allow_pickle=False) as archive:
        expected = archive['energy_vertex_gradient']
        directions = {name: archive[name] for name in archive.files if name.endswith('_moving')}
        if not np.array_equal(archive['state'], data['state']):
            raise ValueError('reference state differs')
    metadata = json.loads(Path(source['chart_identity']['path']).with_suffix('.json').read_text())
    values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    order = original['order']
    assembly = prepare_energy(arrays['faces'], TrianglePanels(arrays['points'][arrays['faces']], 1), order=order)
    field = current_field(arrays, assembly.energy, values, data['state'][7:13],
                          time=source['time_seconds'], acceleration_columns=False)
    phi = field['potential']
    functional = WATER_DENSITY*(phi@assembly.energy@phi/2-phi@field['right'])
    expected_functional = original['operations'][0]['functional_joules']
    np.testing.assert_allclose(functional, expected_functional, rtol=1e-12, atol=1e-8)
    report = dict(complete=False, qualified=False, order=order, rows=[],
        scope='Two alternating matched pairs of complete vertex-gradient calls on one fixed actual state; no trajectory or full-RHS speed claim.',
        input_sha256={str(path): sha256(path.read_bytes()).hexdigest()
            for path in (reference, reference.with_suffix('.npz'), fixture)},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_surface_gradient_cost.py', 'boat_contact_reaction.py', 'boat_world_field.py',
                         'boat_geometry_template.py', 'boat_world_boundary.py', 'canoe.py',
                         'incident_baseline_benchmark.py', 'cuda_operator.py', 'cuda_energy.cu',
                         'cuda_shape_gradient.py', 'cuda_shape.cu', 'cuda_pair.cuh',
                         'cuda_dual.cuh', 'cuda_wide_dual.cuh')})
    started, saved = perf_counter(), {}
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        for index, width in enumerate((1, 3, 3, 1)):
            if perf_counter()-started > 240:
                raise RuntimeError('240s gradient comparison budget exhausted between calls')
            before = perf_counter()
            gradient = vertex_gradient(arrays['points'], arrays['faces'], phi, assembly=assembly,
                density=WATER_DENSITY, order=order, batch_width=width)
            seconds = perf_counter()-before
            delta = gradient-expected
            np.testing.assert_allclose(gradient, expected, rtol=1e-9, atol=1e-7)
            row = dict(index=index, batch_width=width, seconds=seconds,
                maximum_absolute_reference_error=float(np.max(abs(delta))),
                relative_reference_error=float(np.linalg.norm(delta)/np.linalg.norm(expected)),
                maximum_direction_contraction_error_newtons=float(max(abs(np.sum(delta*d)) for d in directions.values())),
                kernel_attributes=_kernels[width].attributes)
            report['rows'].append(row)
            saved[f'gradient_{index}_width{width}'] = gradient
            print(json.dumps(row), flush=True)
            output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        report['paired_width1_over_width3'] = [report['rows'][0]['seconds']/report['rows'][1]['seconds'],
                                               report['rows'][3]['seconds']/report['rows'][2]['seconds']]
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
    run(**vars(parser.parse_args()))
