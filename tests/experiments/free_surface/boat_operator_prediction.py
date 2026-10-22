"""Falsify a geometry-secant operator on actual moving-boat pressure stencils.

Two accepted states train one rank-one geometry map. A later accepted state and
its actual pressure-direction offsets are held out. All current geometry, loads,
incident forcing and prescribed surface direction are shared between models.
Fresh operators/residuals are OFFLINE references, not an online refresh rule.
"""
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
from tests.reference_boat import WATER_DENSITY, GRAVITY, reference_body
from .boat_geometry_template import evaluate_template, InvalidTemplate
from .boat_world_boundary import body_material_samples, interpolate_body_trace
from .boat_world_field import current_field, surface_rates, refine_boundary
from .body_pressure import evaluate_body_pressure, quadrature
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(pilot, output, operator_from=None):
    pilot, output = Path(pilot), fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    recorded = json.loads(pilot.read_text())
    if not recorded.get('complete'):
        raise ValueError('requires a complete saved moving pilot')
    with np.load(pilot.with_suffix('.npz'), allow_pickle=False) as data:
        times, states = data['accepted_times'], data['accepted_states']
    selected = []
    for time in (.05, .10, .15):
        matches = np.flatnonzero(abs(times-time) <= 1e-14)
        if len(matches) != 1:
            raise ValueError('requires unique accepted .05/.10/.15 s states')
        selected.append(int(matches[0]))
    states = states[selected]
    size = (states.shape[1]-13)//2
    report = dict(complete=False, qualified=False,
        scope='One held-out moving-boat geometry/operator/pressure comparison; no ocean, temporal/spatial, cache-policy or throughput qualification',
        anchor_times_seconds=[.05, .10], target_time_seconds=.15, quadrature_order=12,
        input_sha256={str(p): sha256(p.read_bytes()).hexdigest()
                      for p in (pilot, pilot.with_suffix('.npz'))},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_operator_prediction.py', 'boat_geometry_template.py',
                         'boat_world_field.py', 'boat_world_boundary.py', 'body_pressure.py',
                         'canoe.py', 'incident_baseline_benchmark.py',
                         'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh')},
        phases=[], models=[])
    previous = None
    if operator_from is not None:
        operator_from = Path(operator_from)
        source = json.loads(operator_from.read_text())
        if not source.get('complete') or source['input_sha256'] != report['input_sha256']:
            raise ValueError('operator replay requires the same complete moving-state audit')
        for name, expected in source['source_sha256'].items():
            if name != Path(__file__).name and report['source_sha256'][name] != expected:
                raise ValueError(f'operator replay helper changed: {name}')
        with np.load(operator_from.with_suffix('.npz'), allow_pickle=False) as data:
            previous = dict(data)
        report['operator_replay'] = dict(path=str(operator_from),
            sha256=sha256(operator_from.read_bytes()).hexdigest(),
            arrays_sha256=sha256(operator_from.with_suffix('.npz').read_bytes()).hexdigest(),
            scope='Offline value/directional-derivative decomposition, no new operator assemblies')
    saved, started = {}, perf_counter()

    def timed(name, function):
        if perf_counter()-started > 120:
            raise RuntimeError('120 s prediction audit budget exhausted between operations')
        before = perf_counter()
        value = function()
        report['phases'].append(dict(operation=name, seconds=perf_counter()-before))
        return value

    def inputs(y):
        return dict(pose=InitialState(y[:3], y[3:7]/np.linalg.norm(y[3:7])),
                    parent_eta=y[13:13+size], parent_psi=y[13+size:])

    try:
        charts = []
        for row in recorded['geometry']:
            path = Path(row['path']).with_suffix('.npz')
            if sha256(path.read_bytes()).hexdigest() != row['sha256']:
                raise ValueError('saved geometry differs from pilot manifest')
            with np.load(path, allow_pickle=False) as data:
                charts.append((path, dict(data)))
        for chart_path, chart in charts:
            try:
                base = [evaluate_template(chart, **inputs(y))[0] for y in states]
            except InvalidTemplate:
                continue
            break
        else:
            raise ValueError('no common guarded geometry chart for the three states')
        report['chart'] = dict(path=str(chart_path), sha256=sha256(chart_path.read_bytes()).hexdigest())
        body = reference_body()
        if (not np.array_equal(chart['source_body_vertices'], body.mesh.vertices)
                or not np.array_equal(chart['source_body_faces'], body.mesh.faces)
                or not np.array_equal(chart['source_body_cg'], body.cg)):
            raise ValueError('saved chart differs from current canonical body')
        geometry_metadata = json.loads(chart_path.with_suffix('.json').read_text())
        values = coefficients(sea_states.loads_realization(json.dumps(geometry_metadata['actual_sea_realization'])))

        def refine(arrays):
            for _ in range(2):
                arrays = refine_boundary(arrays, region='body')
            return arrays

        base = [refine(a) for a in base]
        reference_ids = ('faces', 'body_faces', 'free_faces', 'free_nodes', 'body_source_group',
                         'free_source_parent', 'free_parent_indices')

        def same_basis(arrays):
            if any(not np.array_equal(arrays[key], base[0][key]) for key in reference_ids):
                raise ValueError('prediction requires identical chart basis and boundary ownership')

        for arrays in base:
            same_basis(arrays)

        def assemble(name, arrays, key):
            if previous is not None:
                return previous[key]
            from .cuda_operator import prepare_energy
            return timed(name, lambda: prepare_energy(arrays['faces'],
                TrianglePanels(arrays['points'][arrays['faces']], 1), order=12).energy)

        energies = [assemble(f'anchor{i}-energy', arrays, f'anchor{i}_energy')
                    for i, arrays in enumerate(base[:2])]
        target, y = base[2], states[2]
        if previous is not None:
            if (not np.array_equal(previous['selected_states'], states)
                    or any(not np.allclose(previous[key], arrays['points'], rtol=0, atol=1e-12)
                           for key, arrays in zip(('anchor0_points', 'anchor1_points', 'target_points'), base))
                    or source['chart'] != report['chart']):
                raise ValueError('operator replay differs from current chart/states/points')
        target_energy = assemble('target-energy', target, 'target_energy')
        twist, pose = y[7:13], inputs(y)['pose']
        actual = timed('target-field', lambda: current_field(target, target_energy, values, twist, time=.15))
        eta_rate, psi_rate, _ = surface_rates(target['points'], target['faces'][target['free_faces']],
            target['free_nodes'], actual['parent_map'], target['parent_psi'], actual['flux'])
        offsets = {}
        offset_energies = {}
        for h in (.001, .0005):
            for sign in (1, -1):
                dt = sign*h
                orientation = (Rotation.from_rotvec(dt*twist[3:])
                    *Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
                arrays = timed(f'offset-{dt}-geometry', lambda: refine(evaluate_template(chart,
                    pose=InitialState(pose.position+dt*twist[:3], orientation),
                    parent_eta=y[13:13+size]+dt*eta_rate,
                    parent_psi=y[13+size:]+dt*psi_rate)[0]))
                same_basis(arrays)
                offsets[dt] = arrays
                offset_energies[dt] = assemble(f'offset-{dt}-energy', arrays, f'offset_{dt}_energy')

        # Fixed lumped AREA weights define one explicit geometry projection.
        # Their derivative is zero: never silently change the metric at offsets.
        def nodal_area(arrays):
            panels = TrianglePanels(arrays['points'][arrays['faces']], 1)
            return np.bincount(arrays['faces'].ravel(),
                weights=np.repeat(panels.areas/3, 3), minlength=len(arrays['points']))

        weight = (nodal_area(base[0])+nodal_area(base[1]))/2
        weight /= weight.sum()
        origin, delta = base[0]['points'], base[1]['points']-base[0]['points']

        def inner(a, b):
            return float(np.sum(weight[:, None]*a*b))

        denominator = inner(delta, delta)
        if denominator <= 0:
            raise ValueError('zero anchor geometry displacement')

        def prediction(arrays):
            change = arrays['points']-origin
            alpha = inner(change, delta)/denominator
            return energies[0]+alpha*(energies[1]-energies[0]), dict(alpha=alpha,
                orthogonal_geometry_rms_m=float(np.sqrt(inner(change-alpha*delta, change-alpha*delta))))

        geometry_rate = (offsets[.0005]['points']-offsets[-.0005]['points'])/.001
        alignment = inner(delta, geometry_rate)/np.sqrt(denominator*inner(geometry_rate, geometry_rate))
        report.update(anchor_geometry_rms_m=float(np.sqrt(denominator)),
                      pressure_direction_geometry_cosine=float(alignment),
                      faces=len(target['faces']), nodes=len(target['points']))
        samples, groups = body_material_samples(target, quadrature()[0])
        hints = np.broadcast_to(np.arange(len(target['body_faces']))[:, None], groups.shape)
        rotation = rotation_matrix(pose.orientation)
        mesh_velocity = twist[:3]+np.cross(twist[3:], target['points']-pose.position)

        def surface_rms(delta_parent):
            faces = target['faces'][target['free_faces']]
            triangles = target['points'][faces]
            area = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])[:, 2]/2
            local = np.full(len(target['points']), -1)
            local[target['free_nodes']] = np.arange(len(target['free_nodes']))
            nodal = (actual['parent_map']@delta_parent)[local[faces]]
            square = (np.sum(nodal*nodal, axis=1)+np.sum(nodal, axis=1)**2)/12
            return float(np.sqrt(area@square/area.sum()))

        names = ['fresh', 'frozen', 'secant']
        if previous is not None:
            names += ['exact_center_secant_slope', 'secant_center_exact_slope', 'exact_center_exact_slope']
        predicted_target, predicted_description = prediction(target)
        for name in names:
            row = dict(name=name, complete=False, pressure=[])
            report['models'].append(row)
            before = perf_counter()
            try:
                if name in ('fresh', 'exact_center_secant_slope', 'exact_center_exact_slope'):
                    energy, description = target_energy, {}
                elif name == 'frozen':
                    energy, description = energies[1], {}
                else:
                    energy, description = predicted_target, predicted_description
                if (not np.isfinite(energy).all()
                        or np.linalg.norm(energy-energy.T) > 1e-10*np.linalg.norm(energy)):
                    raise ValueError('nonfinite or asymmetric predicted energy')
                field = current_field(target, energy, values, twist, time=.15)
                defect = (target_energy@field['potential']-field['right'])
                unknown = np.setdiff1d(np.arange(len(energy)), target['free_nodes'])
                eta, psi, _ = surface_rates(target['points'], target['faces'][target['free_faces']],
                    target['free_nodes'], field['parent_map'], target['parent_psi'], field['flux'])
                row.update(geometry_prediction=description,
                    true_current_interior_relative_defect=float(np.linalg.norm(defect[unknown])
                        /max(np.linalg.norm(field['right'][unknown]), 1e-30)),
                    parent_height_rate_error_surface_rms_m_per_s=surface_rms(eta-eta_rate),
                    parent_potential_rate_error_surface_rms_m2_per_s2=surface_rms(psi-psi_rate))
                for h in (.001, .0005):
                    traces = []
                    for dt in (h, -h):
                        arrays = offsets[dt]
                        if name == 'fresh':
                            offset_energy, offset_description = offset_energies[dt], {}
                        elif name == 'frozen':
                            offset_energy, offset_description = energies[1], {}
                        elif name == 'exact_center_secant_slope':
                            predicted_offset, offset_description = prediction(arrays)
                            offset_energy = target_energy+predicted_offset-predicted_target
                        elif name in ('secant_center_exact_slope', 'exact_center_exact_slope'):
                            exact_slope = (offset_energies[h]-offset_energies[-h])/(2*h)
                            offset_energy, offset_description = energy+dt*exact_slope, {}
                        else:
                            offset_energy, offset_description = prediction(arrays)
                        offset_field = current_field(arrays, offset_energy, values, twist,
                                                     acceleration_columns=False, time=.15+dt)
                        traces.append(interpolate_body_trace(arrays, offset_field['potential'], samples,
                                                             groups, face_hints=hints))
                        if offset_description:
                            row.setdefault('offset_predictions', {})[str(dt)] = offset_description
                    pressure = evaluate_body_pressure(target['points'], target['faces'], target['body_faces'],
                        pose.position, body.mass, rotation@body.inertia_cg@rotation.T, twist,
                        field['potential'], (traces[0]-traces[1])/(2*h), mesh_velocity,
                        field['acceleration_potentials'], density=WATER_DENSITY, gravity=GRAVITY)
                    row['pressure'].append(dict(step=h, load_world=pressure.hydrodynamic_load.tolist(),
                        acceleration_world=pressure.acceleration.tolist(), added_mass=pressure.added_mass.tolist()))
                row.update(complete=True, seconds=perf_counter()-before)
                saved[name+'_potential'] = field['potential']
            except Exception as error:
                row['failure'] = dict(type=type(error).__name__, message=str(error))
                if name == 'fresh':
                    raise
            print(json.dumps(row), flush=True)
        saved.update(anchor0_energy=energies[0], anchor1_energy=energies[1], target_energy=target_energy,
            anchor0_points=origin, anchor1_points=base[1]['points'], target_points=target['points'],
            selected_states=states, geometry_weight=weight)
        for dt, energy in offset_energies.items():
            saved[f'offset_{dt}_energy'] = energy
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
    parser.add_argument('--operator-from', help='Complete saved audit for offline value/slope decomposition')
    run(**vars(parser.parse_args()))
