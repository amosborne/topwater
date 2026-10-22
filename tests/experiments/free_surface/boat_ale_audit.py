"""Fresh h/h2 pressure comparison on a saved canonical-boat state and chart."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.spatial.transform import Rotation

from topwater import InitialState, sea_states
from topwater.rankine import TrianglePanels
from tests.reference_boat import reference_body
from .boat_geometry_template import evaluate_template, InvalidTemplate
from .boat_world_field import current_field, surface_rates, refine_boundary
from .boat_pressure_stencil import evaluate
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(pilot, output, state_case=None, time=.2, step=.001, chart_path=None):
    from .cuda_operator import prepare_energy
    pilot, output = Path(pilot), fresh(output)
    source = json.loads(pilot.read_text())
    if not np.isfinite(step) or step <= 0:
        raise ValueError('positive finite derivative step required')
    paths = [pilot, pilot.with_suffix('.npz')]
    if state_case:
        paths.append(Path(state_case))
        case = json.loads(Path(state_case).read_text())
        state, time = np.array(case['state']), case['time_seconds']
    else:
        with np.load(pilot.with_suffix('.npz'), allow_pickle=False) as data:
            matches = np.flatnonzero(np.abs(data['accepted_times']-time) < 1e-12)
            if len(matches) != 1:
                raise ValueError('requested accepted time missing or ambiguous')
            state = data['accepted_states'][matches[0]]
    fixture = Path(__file__).resolve().parents[2]/'reference_boat.py'
    fixture_hash = sha256(fixture.read_bytes()).hexdigest()
    if fixture_hash != source['fixture_sha256']:
        raise ValueError('pilot fixture differs from current mass/geometry')
    if state_case and case.get('fixture_sha256') != fixture_hash:
        raise ValueError('state fixture differs from current mass/geometry')
    if source['boundary_refinement']['region'] != 'body':
        raise ValueError('audit supports body-only refinement')
    report = dict(complete=False, qualified=False,
        scope='Prescribed nonlinear boat state, fresh same-chart h/h2 pressure; no trajectory or open-water qualification.',
        time_seconds=time, derivative_step_seconds=step, fixture_sha256=fixture_hash, phases=[], comparisons=[],
        input_sha256={str(p): sha256(p.read_bytes()).hexdigest() for p in paths},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_ale_audit.py', 'boat_ale_pressure.py', 'boat_geometry_template.py',
                         'boat_pressure_stencil.py', 'boat_world_field.py', 'boat_world_boundary.py',
                         'body_pressure.py', 'incident_baseline_benchmark.py', 'canoe.py',
                         'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh')})
    output.parent.mkdir(parents=True, exist_ok=True)
    started, saved = perf_counter(), dict(state=state)

    def energy(label, arrays):
        if perf_counter()-started > 90:
            raise RuntimeError('90 s audit budget exhausted')
        before = perf_counter()
        assembly = prepare_energy(arrays['faces'], TrianglePanels(arrays['points'][arrays['faces']], 1),
                                  order=source['operator_policy']['quadrature_order'])
        if not assembly.uniform_quadrature:
            raise ValueError('uniform quadrature required')
        saved[label] = assembly.energy
        report['phases'].append(dict(operation=label, seconds=perf_counter()-before))
        return assembly.energy

    def refine(arrays):
        for _ in range(source['boundary_refinement']['levels']):
            arrays = refine_boundary(arrays, region='body')
        return arrays

    try:
        size = (len(state)-13)//2
        pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
        for row in source['geometry']:
            path = Path(row['path']).with_suffix('.npz')
            if chart_path and path.resolve() != Path(chart_path).with_suffix('.npz').resolve():
                continue
            if sha256(path.read_bytes()).hexdigest() != row['sha256']:
                raise ValueError('pilot chart hash changed')
            with np.load(path, allow_pickle=False) as data:
                chart = dict(data)
            try:
                arrays, _ = evaluate_template(chart, pose=pose, parent_eta=state[13:13+size],
                                               parent_psi=state[13+size:])
            except InvalidTemplate:
                continue
            break
        else:
            raise ValueError('no saved chart covers state')
        body = reference_body()
        for name, expected in (('source_body_vertices', body.mesh.vertices),
                               ('source_body_faces', body.mesh.faces), ('source_body_cg', body.cg)):
            if not np.array_equal(chart[name], expected):
                raise ValueError('saved chart differs from canonical body')
        metadata_path = path.with_suffix('.json')
        metadata = json.loads(metadata_path.read_text())
        if metadata['fixture_sha256'] != fixture_hash:
            raise ValueError('chart fixture differs')
        report['input_sha256'].update({str(p): sha256(p.read_bytes()).hexdigest()
                                     for p in (path, metadata_path)})
        values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
        arrays = refine(arrays)
        center = energy('center_energy', arrays)
        field = current_field(arrays, center, values, state[7:13], time=time, acceleration_columns=False)
        er, pr, _ = surface_rates(arrays['points'], arrays['faces'][arrays['free_faces']],
            arrays['free_nodes'], field['parent_map'], arrays['parent_psi'], field['flux'])
        offsets, energies = {}, {}
        for dt in (step, -step, step/2, -step/2):
            orientation = (Rotation.from_rotvec(dt*state[10:13])
                *Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
            offsets[dt] = refine(evaluate_template(chart,
                pose=InitialState(pose.position+dt*state[7:10], orientation),
                parent_eta=state[13:13+size]+dt*er, parent_psi=state[13+size:]+dt*pr)[0])
            energies[dt] = energy(f'offset_{dt}_energy', offsets[dt])
        for charted in (False, True):
            row = dict(representation='charted nodal' if charted else 'material quadrature')
            report['comparisons'].append(row)
            try:
                row['pressures'] = evaluate(body, arrays, state, values, time, center,
                                            offsets, energies, charted=charted)[3]
            except ValueError as error:
                if charted:
                    raise
                row['failure'] = str(error)
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
    parser.add_argument('--pilot', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--state-case')
    parser.add_argument('--time', type=float, default=.2)
    parser.add_argument('--step', type=float, default=.001)
    parser.add_argument('--chart-path', help='Explicit chart owned by the saved pilot; no arbitrary replacement.')
    run(**vars(parser.parse_args()))
