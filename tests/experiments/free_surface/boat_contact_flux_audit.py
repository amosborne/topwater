"""One actual cut-boat snapshot: quadrature and compatible affine-flux controls."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from topwater import InitialState, sea_states
from topwater.rankine import TrianglePanels
from tests.reference_boat import reference_body
from .boat_geometry_template import evaluate_template, InvalidTemplate
from .boat_world_field import current_field, surface_rates, refine_boundary
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(pilot, output, chart_path, composite=False, orders=(6, 12, 24), edge_graded=False,
        free_refinement=0, local_body_refinement=0):
    from . import cuda_operator as gpu
    from .cuda_operator import prepare_energy
    pilot, output = Path(pilot), fresh(output)
    if not orders or orders[0] != 6 or any(order not in (6, 12, 24, 48, 96) for order in orders):
        raise ValueError('begin with replay order six; bounded orders are 6/12/24/48/96')
    if free_refinement not in (0, 1):
        raise ValueError('bounded spatial control allows zero or one free-surface refinement')
    if local_body_refinement not in (0, 1, 2) or (free_refinement and local_body_refinement):
        raise ValueError('separate free refinement from zero/one/two local body refinements')
    output.parent.mkdir(parents=True, exist_ok=True)
    source = json.loads(pilot.read_text())
    fixture = Path(__file__).resolve().parents[2]/'reference_boat.py'
    if sha256(fixture.read_bytes()).hexdigest() != source['fixture_sha256']:
        raise ValueError('current boat mass/geometry differs from pilot')
    with np.load(pilot.with_suffix('.npz'), allow_pickle=False) as data:
        state, time = data['accepted_states'][-1], float(data['accepted_times'][-1])
        if not np.array_equal(state, data['last_rhs_state']):
            raise ValueError('requires an accepted endpoint matching the last completed RHS')
    old = next(row for row in reversed(source['rhs']) if row['complete'])
    if old['time'] != time or source['boundary_refinement']['region'] != 'body':
        raise ValueError('requires matching endpoint and body-only refinement')
    report = dict(complete=False, qualified=False,
        scope='Actual accepted diagnostic-tank state; quadrature and affine harmonic controls, no trajectory or open-water qualification.',
        time_seconds=time, original_surface=old['surface'], rows=[],
        input_sha256={str(p): sha256(p.read_bytes()).hexdigest() for p in (pilot, pilot.with_suffix('.npz'))},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_contact_flux_audit.py', 'boat_contact_trace.py', 'boat_target_quadrature.py', 'boat_world_field.py', 'boat_world_boundary.py',
                         'boat_parent_clip_proof.py', 'boat_geometry_template.py', 'canoe.py',
                         'incident_baseline_benchmark.py', 'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh')})
    saved, started = dict(state=state), perf_counter()
    try:
        pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
        size = (len(state)-13)//2
        selected_path = Path(chart_path).with_suffix('.npz').resolve()
        owned = [item for item in source['geometry']
                 if Path(item['path']).with_suffix('.npz').resolve() == selected_path]
        if len(owned) != 1:
            raise ValueError('explicit chart must belong uniquely to the pilot geometry manifest')
        for item in owned:
            path = Path(item['path']).with_suffix('.npz')
            if sha256(path.read_bytes()).hexdigest() != item['sha256']:
                raise ValueError('saved geometry hash changed')
            with np.load(path, allow_pickle=False) as data:
                chart = dict(data)
            try:
                arrays, _ = evaluate_template(chart, pose=pose, parent_eta=state[13:13+size],
                                               parent_psi=state[13+size:])
            except InvalidTemplate:
                continue
            for _ in range(source['boundary_refinement']['levels']):
                arrays = refine_boundary(arrays, region='body')
            if len(arrays['faces']) == old['geometry'][0]['faces']:
                break
        else:
            raise ValueError('explicit saved chart does not cover the accepted endpoint')
        report['chart_identity'] = dict(path=str(path), node_chart_key=str(chart['node_chart_key']),
            selection='Explicit owned chart, hash verified; legacy pilot has no recorded base chart key. Order-six surface replay checked below.')
        body = reference_body()
        for name, expected in (('source_body_vertices', body.mesh.vertices),
                               ('source_body_faces', body.mesh.faces), ('source_body_cg', body.cg)):
            if not np.array_equal(chart[name], expected):
                raise ValueError('saved chart differs from canonical body')
        metadata_path = path.with_suffix('.json')
        metadata = json.loads(metadata_path.read_text())
        if metadata['fixture_sha256'] != source['fixture_sha256']:
            raise ValueError('chart fixture differs')
        report['input_sha256'].update({str(p): sha256(p.read_bytes()).hexdigest() for p in (path, metadata_path)})
        values = coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))

        def check_replay(surface):
            for key in ('maximum_fine_vertical_velocity_at_quadrature',
                        'maximum_parent_height_rate', 'maximum_parent_potential_rate'):
                if not np.isclose(surface[key], old['surface'][key], rtol=1e-10, atol=1e-10):
                    raise ValueError('explicit chart does not reproduce saved surface: '+key)

        if free_refinement or local_body_refinement:
            base_panels = TrianglePanels(arrays['points'][arrays['faces']], 1)
            base_energy = prepare_energy(arrays['faces'], base_panels, order=6).energy
            base_field = current_field(arrays, base_energy, values, state[7:13], time=time)
            _, _, base_surface = surface_rates(arrays['points'], arrays['faces'][arrays['free_faces']],
                arrays['free_nodes'], base_field['parent_map'], arrays['parent_psi'], base_field['flux'])
            check_replay(base_surface)
            if free_refinement:
                arrays = refine_boundary(arrays, region='free')
            else:
                contact_nodes = np.intersect1d(arrays['free_nodes'],
                    np.unique(arrays['faces'][arrays['body_faces']]))
                location = np.array(base_surface['peak_vertical_location_world'])
                node = int(contact_nodes[np.argmin(np.linalg.norm(arrays['points'][contact_nodes]-location, axis=1))])
                report['local_body_refinement'] = dict(levels=local_body_refinement, contact_node=node,
                    rule='Body faces incident to the nearest contact node to the ORIGINAL q6 vertical peak; same node at each level.')
                for _ in range(local_body_refinement):
                    selected = arrays['body_faces'][np.any(arrays['faces'][arrays['body_faces']] == node, axis=1)]
                    arrays = refine_boundary(arrays, region='body', marked_faces=selected)
        report['extra_free_refinement'] = free_refinement
        report['boundary_counts'] = dict(nodes=len(arrays['points']), faces=len(arrays['faces']),
            free_faces=len(arrays['free_faces']), body_faces=len(arrays['body_faces']))
        points, faces, fixed = arrays['points'], arrays['faces'], arrays['free_nodes']
        panels = TrianglePanels(points[faces], 1)
        free_faces = faces[arrays['free_faces']]
        local = np.full(len(points), -1, dtype=int)
        local[fixed] = np.arange(len(fixed))
        projected = points[free_faces, :2]
        area = np.linalg.det(projected[:, 1:]-projected[:, :1])/2
        mass = np.zeros((len(fixed), len(fixed)))
        for i in range(3):
            for j in range(3):
                np.add.at(mass, (local[free_faces[:, i]], local[free_faces[:, j]]), area*(2 if i == j else 1)/12)
        mass_factor = cho_factor(mass)
        unknown = np.setdiff1d(np.arange(len(points)), fixed)
        body_edges = {tuple(sorted((int(a), int(b)))) for face in faces[arrays['body_faces']]
                      for a, b in zip(face, np.roll(face, -1))}
        edges = Counter(tuple(sorted((int(a), int(b)))) for face in free_faces
                        for a, b in zip(face, np.roll(face, -1)))
        contact_edges = np.array([edge for edge, count in edges.items() if count == 1 and edge in body_edges])
        if not len(contact_edges):
            raise ValueError('no shared body/free contact edges')
        saved.update(points=points, faces=faces, free_nodes=fixed, free_mass=mass, contact_edges=contact_edges)
        a, constant = np.array([.3, -.2, .4]), .7
        exact_phi = points@a+constant
        boundary = np.zeros(len(points))
        free_boundary = np.zeros(len(points))
        face_flux = panels.areas*(panels.normals@a)/3
        for corner in range(3):
            np.add.at(boundary, faces[:, corner], face_flux)
            np.add.at(free_boundary, free_faces[:, corner], face_flux[arrays['free_faces']])
        nonfree_boundary = boundary-free_boundary
        expected_flux = free_boundary[fixed]
        expected_density = cho_solve(mass_factor, expected_flux)
        rules = [(order, 0, 0) for order in orders]
        if composite:
            rules += [(6, 1, 0), (6, 2, 0)]
        if edge_graded:
            rules += [(6, 0, 2), (6, 0, 4)]
        for order, levels, edge_power in rules:
            if perf_counter()-started > 90:
                raise RuntimeError('90s contact audit budget exhausted')
            before = perf_counter()
            captured = {}
            original_pairs = gpu._pairs

            def capture_pairs(grid, block, args):
                original_pairs(grid, block, args)
                captured['v'] = args[-2].get()
                captured['h'] = args[-1].get()

            gpu._pairs = capture_pairs
            try:
                if edge_power:
                    from .boat_target_quadrature import prepare_edge_energy
                    assembly = prepare_edge_energy(faces, panels, order, edge_power)
                elif levels:
                    from .boat_target_quadrature import prepare_composite_energy
                    assembly = prepare_composite_energy(faces, panels, order, levels)
                else:
                    assembly = prepare_energy(faces, panels, order=order)
            finally:
                gpu._pairs = original_pairs
            energy = assembly.energy
            if not np.isfinite(energy).all():
                raise ValueError('nonfinite current energy')
            field = current_field(arrays, energy, values, state[7:13], time=time)
            er, pr, surface = surface_rates(points, free_faces, fixed, field['parent_map'],
                                           arrays['parent_psi'], field['flux'])
            if order == 6 and levels == 0 and edge_power == 0 and not (free_refinement or local_body_refinement):
                check_replay(surface)
            peak = np.array(surface['peak_vertical_location_world'])
            start, direction = points[contact_edges[:, 0]], points[contact_edges[:, 1]]-points[contact_edges[:, 0]]
            parameter = np.clip(np.sum((peak-start)*direction, axis=1)/np.sum(direction**2, axis=1), 0, 1)
            surface['peak_distance_to_contact_m'] = float(np.linalg.norm(peak-start-parameter[:, None]*direction, axis=1).min())
            surface['peak_free_parent'] = int(arrays['free_source_parent'][surface['peak_vertical_free_face_local']])
            lifting = np.zeros(len(points)); lifting[fixed] = exact_phi[fixed]
            harmonic = lifting.copy()
            harmonic[unknown] = cho_solve(cho_factor(energy[np.ix_(unknown, unknown)]),
                                         (nonfree_boundary-energy@lifting)[unknown])
            harmonic_flux = (energy@harmonic-nonfree_boundary)[fixed]
            harmonic_density = cho_solve(mass_factor, harmonic_flux)
            defect = energy@exact_phi-boundary
            delta = harmonic_density-expected_density
            _, harmonic_rate, harmonic_surface = surface_rates(points, free_faces, fixed, field['parent_map'],
                arrays['parent_nodes']@a+constant, harmonic_flux, gravity=0.)
            _, expected_rate, expected_surface = surface_rates(points, free_faces, fixed, field['parent_map'],
                arrays['parent_nodes']@a+constant, expected_flux, gravity=0.)
            radiation = (energy@field['acceleration_potentials']-field['body_load'])[fixed]
            contributions = cho_solve(mass_factor, radiation*state[7:13])
            actual_density = cho_solve(mass_factor, field['flux'])
            stationary_density = actual_density-contributions.sum(axis=1)
            raw_v, raw_h = captured['v'], captured['h']
            from .boat_contact_trace import compare_traces
            normal_traces = cho_solve((assembly.factor.get(), True),
                                     raw_h@np.column_stack((field['potential'], harmonic, exact_phi)))
            actual_trace_comparison = compare_traces(points, free_faces, fixed,
                field['potential'], actual_density, normal_traces[arrays['free_faces'], 0], contact_edges)
            affine_trace_comparison = compare_traces(points, free_faces, fixed,
                harmonic, harmonic_density, normal_traces[arrays['free_faces'], 1], contact_edges, gravity=0.)
            exact_q = panels.normals@a
            green = raw_h@exact_phi-raw_v@exact_q
            reciprocal_change = (raw_v.T-raw_v)@exact_q/2
            green_scale = np.linalg.norm(raw_v@exact_q)
            opposite = np.roll(panels.triangles, -2, axis=1)-np.roll(panels.triangles, -1, axis=1)
            gradient = np.cross(panels.geometric_normals[:, None], opposite)/(2*panels.areas[:, None, None])
            curl = np.cross(panels.geometric_normals[:, None], gradient)
            curl_phi = np.einsum('fck,fc->fk', curl, exact_phi[faces])

            def weak_curl_action(matrix):
                face_vector = matrix@curl_phi
                result = np.zeros(len(points))
                for corner in range(3):
                    np.add.at(result, faces[:, corner], np.sum(curl[:, corner]*face_vector, axis=1))
                return result

            weak_raw = weak_curl_action(raw_v)+raw_h.T@exact_q-boundary
            weak_sym = weak_curl_action((raw_v+raw_v.T)/2)+raw_h.T@exact_q-boundary
            propagated = raw_h.T@cho_solve((assembly.factor.get(), True), green-reciprocal_change)
            decomposition_error = defect-weak_sym-propagated
            pair_control = dict(
                raw_green_relative=float(np.linalg.norm(green)/green_scale),
                symmetrized_green_relative=float(np.linalg.norm(green-reciprocal_change)/green_scale),
                reciprocity_action_relative=float(np.linalg.norm(reciprocal_change)/green_scale),
                raw_green_maximum=float(np.max(abs(green))),
                reciprocity_action_maximum=float(np.max(abs(reciprocal_change))),
                raw_weak_curl_relative=float(np.linalg.norm(weak_raw)/np.linalg.norm(boundary)),
                symmetric_weak_curl_relative=float(np.linalg.norm(weak_sym)/np.linalg.norm(boundary)),
                propagated_auxiliary_relative=float(np.linalg.norm(propagated)/np.linalg.norm(boundary)),
                full_defect_decomposition_maximum=float(np.max(abs(decomposition_error))))
            report['rows'].append(dict(order=order, target_subdivision_levels=levels,
                target_edge_power=edge_power,
                actual_trace_comparison=actual_trace_comparison,
                affine_trace_comparison=affine_trace_comparison,
                affine_auxiliary_exact_trace_error_maximum=float(np.max(abs(normal_traces[:, 2]-exact_q))),
                affine_pair_control=pair_control,
                integration_diagnostics=assembly.diagnostics,
                energy_constant_maximum=float(np.max(abs(energy@np.ones(len(points))))),
                seconds=perf_counter()-before,
                actual_surface=surface, relative_field_residual=field['relative_field_residual'],
                affine_control=dict(velocity=a.tolist(), gauge_constant=constant,
                    maximum_potential_error=float(np.max(abs(harmonic-exact_phi))),
                    weak_consistency_relative=float(np.linalg.norm(defect)/np.linalg.norm(boundary)),
                    maximum_graph_flux_error=float(np.max(abs(delta))),
                    graph_flux_mass_relative=float(np.sqrt(delta@mass@delta/(expected_density@mass@expected_density))),
                    exact_vertical_velocity_m_per_s=float(a[2]),
                    recovered_surface=harmonic_surface, independently_projected_surface=expected_surface,
                    maximum_parent_rate_error=float(np.max(abs(harmonic_rate-expected_rate)))),
                body_motion_decomposition=dict(
                    maximum_absolute_radiation_contributions_m_per_s=np.max(abs(contributions), axis=0).tolist(),
                    maximum_absolute_stationary_density_m_per_s=float(np.max(abs(stationary_density))))))
            saved.update({f'q{order}_sub{levels}_edge{edge_power}_{name}': value for name, value in dict(energy=energy,
                potential=field['potential'], flux=field['flux'], eta_rate=er, psi_rate=pr,
                harmonic_potential=harmonic, harmonic_flux=harmonic_flux, expected_flux=expected_flux,
                raw_green_defect=green, reciprocity_action=reciprocal_change,
                auxiliary_normal_traces=normal_traces,
                raw_weak_curl_defect=weak_raw, symmetric_weak_curl_defect=weak_sym,
                propagated_auxiliary_defect=propagated,
                radiation_density_contributions=contributions, stationary_density=stationary_density).items()})
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
    parser.add_argument('--chart-path', required=True)
    parser.add_argument('--composite', action='store_true')
    parser.add_argument('--orders', type=int, nargs='+', default=(6, 12, 24))
    parser.add_argument('--edge-graded', action='store_true')
    parser.add_argument('--free-refinement', type=int, choices=(0, 1), default=0)
    parser.add_argument('--local-body-refinement', type=int, choices=(0, 1, 2), default=0)
    run(**vars(parser.parse_args()))
