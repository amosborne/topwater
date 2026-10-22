"""Fresh full RHS at a reused boat endpoint, separating policy from state error."""
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


def run(pilot, output):
    pilot, output = Path(pilot), fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    source = json.loads(pilot.read_text())
    if not source.get('complete') or source['rhs'][-1].get('operator_mode') != 'reused':
        raise ValueError('requires a complete pilot with a reused endpoint RHS')
    if source['controls']['derivative_step_seconds'] != .001:
        raise ValueError('endpoint audit requires the declared 1 ms pressure stencil')
    if source['boundary_refinement']['region'] != 'body':
        raise ValueError('endpoint audit requires body-only refinement')
    fixture_hash = sha256(Path(__file__).resolve().parents[2].joinpath('reference_boat.py').read_bytes()).hexdigest()
    if fixture_hash != source['fixture_sha256']:
        raise ValueError('boat fixture including mass/inertia differs from pilot')
    with np.load(pilot.with_suffix('.npz'), allow_pickle=False) as data:
        state, old_rate, time = data['accepted_states'][-1], data['last_rhs_rate'], float(data['accepted_times'][-1])
        if not np.array_equal(state, data['last_rhs_state']) or source['rhs'][-1]['time'] != time:
            raise ValueError('endpoint state/time do not match saved FSAL diagnostics')
    report = dict(complete=False, qualified=False,
        scope='One fresh full nonlinear RHS at the cached pilot exact endpoint; not trajectory, spatial or ocean qualification.',
        input_sha256={str(p): sha256(p.read_bytes()).hexdigest() for p in (pilot, pilot.with_suffix('.npz'))},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in (
            'boat_endpoint_audit.py', 'boat_geometry_template.py', 'boat_world_boundary.py',
            'boat_world_field.py', 'body_pressure.py', 'boat_pressure_stencil.py',
            'incident_baseline_benchmark.py', 'canoe.py', 'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh')},
        time_seconds=time, fixture_sha256=fixture_hash,
        quadrature_order=source['operator_policy']['quadrature_order'], phases=[])
    for name, digest in report['source_sha256'].items():
        if name in source['source_sha256'] and source['source_sha256'][name] != digest:
            raise ValueError('pilot physical helper changed: '+name)
    started, saved = perf_counter(), dict(state=state)

    def timed(name, function):
        if perf_counter()-started > 90:
            raise RuntimeError('90 s endpoint audit budget exhausted between operations')
        before = perf_counter()
        result = function()
        report['phases'].append(dict(operation=name, seconds=perf_counter()-before))
        return result

    try:
        size = (len(state)-13)//2
        pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
        valid = []
        for row in source['geometry']:
            path = Path(row['path']).with_suffix('.npz')
            if sha256(path.read_bytes()).hexdigest() != row['sha256']:
                raise ValueError('pilot chart changed')
            with np.load(path, allow_pickle=False) as data:
                chart = dict(data)
            try:
                arrays, _ = evaluate_template(chart, pose=pose, parent_eta=state[13:13+size],
                                               parent_psi=state[13+size:])
            except InvalidTemplate:
                continue
            valid.append((path, chart, arrays))
        if not valid:
            raise ValueError('no saved chart covers endpoint')
        path, chart, arrays = valid[0]
        for _, _, other in valid[1:]:
            if (not np.array_equal(arrays['faces'], other['faces'])
                    or not np.array_equal(arrays['body_source_group'], other['body_source_group'])
                    or not np.allclose(arrays['points'], other['points'], rtol=0, atol=1e-12)):
                raise ValueError('ambiguous endpoint chart basis')
        body = reference_body()
        if (not np.array_equal(chart['source_body_vertices'], body.mesh.vertices)
                or not np.array_equal(chart['source_body_faces'], body.mesh.faces)
                or not np.array_equal(chart['source_body_cg'], body.cg)):
            raise ValueError('canonical body differs from saved chart')
        metadata_path = path.with_suffix('.json')
        metadata = json.loads(metadata_path.read_text())
        if metadata['fixture_sha256'] != fixture_hash:
            raise ValueError('chart and current/pilot boat fixture hashes differ')
        report['chart'] = dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest(),
                              metadata_sha256=sha256(metadata_path.read_bytes()).hexdigest())
        values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))

        def refine(value):
            for _ in range(source['boundary_refinement']['levels']):
                value = refine_boundary(value, region='body')
            return value

        arrays = refine(arrays)
        if len(arrays['faces']) != source['rhs'][-1]['geometry'][0]['faces']:
            raise ValueError('endpoint refinement differs from pilot')

        def energy(label, value):
            from .cuda_operator import prepare_energy
            assembly = timed(label, lambda: prepare_energy(value['faces'],
                TrianglePanels(value['points'][value['faces']], 1), order=report['quadrature_order']))
            if not assembly.uniform_quadrature:
                raise ValueError('endpoint audit requires uniform quadrature')
            saved[label] = assembly.energy
            return assembly.energy

        center = energy('center_energy', arrays)
        field = current_field(arrays, center, values, state[7:13], time=time, acceleration_columns=False)
        er, pr, _ = surface_rates(arrays['points'], arrays['faces'][arrays['free_faces']],
            arrays['free_nodes'], field['parent_map'], arrays['parent_psi'], field['flux'])
        offsets, energies = {}, {}
        for dt in (.001, -.001):
            orientation = (Rotation.from_rotvec(dt*state[10:13])
                *Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
            moved = refine(evaluate_template(chart,
                pose=InitialState(pose.position+dt*state[7:10], orientation),
                parent_eta=state[13:13+size]+dt*er, parent_psi=state[13+size:]+dt*pr)[0])
            offsets[dt], energies[dt] = moved, energy(f'offset_{dt}_energy', moved)
        _, er, pr, pressures = evaluate(body, arrays, state, values, time, center, offsets, energies)
        raw = state[3:7]
        qdot = .5*np.r_[-state[10:13]@raw[1:], raw[0]*state[10:13]+np.cross(state[10:13], raw[1:])]
        rate = np.r_[state[7:10], qdot, pressures[0]['acceleration_world'], er, pr]
        saved.update(fresh_rate=rate, cached_rate=old_rate, points=arrays['points'], faces=arrays['faces'])
        report.update(complete=True, fresh_pressure=pressures[0],
            cached_minus_fresh_load_world=(np.array(source['rhs'][-1]['hydrodynamic_load_world'])
                                          -pressures[0]['load_world']).tolist(),
            cached_minus_fresh_rate=(old_rate-rate).tolist())
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
    run(**vars(parser.parse_args()))
