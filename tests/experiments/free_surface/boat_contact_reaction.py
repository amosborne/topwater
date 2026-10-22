"""Decompose a saved free reaction without another operator or fluid solve."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from topwater import InitialState
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from .boat_geometry_template import evaluate_template
from .boat_world_field import refine_boundary
from .canoe import body_load
from topwater import HullMesh
from .p4_fluid_ordering import fresh


DEPENDENCIES = ('boat_geometry_template.py', 'boat_world_field.py', 'boat_world_boundary.py', 'canoe.py')


def refine_snapshot(arrays, report, pilot_report):
    for _ in range(pilot_report['boundary_refinement']['levels']):
        arrays = refine_boundary(arrays, region='body')
    for _ in range(report.get('extra_free_refinement', 0)):
        arrays = refine_boundary(arrays, region='free')
    local_refinement = report.get('local_body_refinement', {})
    for _ in range(local_refinement.get('levels', 0)):
        node = local_refinement['contact_node']
        selected = arrays['body_faces'][np.any(arrays['faces'][arrays['body_faces']] == node, axis=1)]
        arrays = refine_boundary(arrays, region='body', marked_faces=selected)
    return arrays


def load_snapshot(audit):
    audit = Path(audit)
    report = json.loads(audit.read_text())
    if not report['complete']:
        raise ValueError('complete owned contact audit required')
    for path, digest in report['input_sha256'].items():
        if sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('contact audit input changed: '+path)
    for name in DEPENDENCIES:
        if sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() != report['source_sha256'][name]:
            raise ValueError('physical helper changed: '+name)
    data = dict(np.load(audit.with_suffix('.npz'), allow_pickle=False))
    state = data['state']
    count = (len(state)-13)//2
    chart_path = Path(report['chart_identity']['path'])
    chart = dict(np.load(chart_path, allow_pickle=False))
    if str(chart['node_chart_key']) != report['chart_identity']['node_chart_key']:
        raise ValueError('chart key changed')
    arrays, _ = evaluate_template(chart, pose=InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7])),
                                  parent_eta=state[13:13+count], parent_psi=state[13+count:])
    # Original audit requires exactly the pilot's body-only refinement.
    pilot = next(Path(path) for path in report['input_sha256']
                 if path.endswith('.json') and '-geometry-' not in path)
    pilot_report = json.loads(pilot.read_text())
    arrays = refine_snapshot(arrays, report, pilot_report)
    for name in ('points', 'faces', 'free_nodes'):
        if not np.array_equal(arrays[name], data[name]):
            raise ValueError('reconstructed saved geometry differs: '+name)
    return report, data, arrays, pilot_report, chart


def run(audit, output, order=48, probe_point=None):
    audit, output = Path(audit), fresh(output)
    report, data, arrays, pilot_report, chart = load_snapshot(audit)
    state = data['state']
    local_refinement = report.get('local_body_refinement', {})
    prefixes = [prefix for prefix in (f'q{order}_sub0_edge0_', f'q{order}_sub0_')
                if prefix+'energy' in data]
    if len(prefixes) != 1:
        raise ValueError('requires one unambiguous saved uniform rule')
    prefix = prefixes[0]
    energy, potential, reaction = (data[prefix+name] for name in ('energy', 'potential', 'flux'))
    auxiliary = data[prefix+'auxiliary_normal_traces'][:, 0]
    points, faces, fixed = data['points'], data['faces'], data['free_nodes']
    panels = TrianglePanels(points[faces], 1)
    free_ids = arrays['free_faces']
    auxiliary_all, auxiliary_free = np.zeros(len(points)), np.zeros(len(points))
    face_load = panels.areas*auxiliary/3
    for corner in range(3):
        np.add.at(auxiliary_all, faces[:, corner], face_load)
        np.add.at(auxiliary_free, faces[free_ids, corner], face_load[free_ids])
    energy_action = energy@potential
    imposed_nonfree = energy_action[fixed]-reaction
    body_ids = arrays['body_faces']
    auxiliary_body = np.zeros(len(points))
    for corner in range(3):
        np.add.at(auxiliary_body, faces[body_ids, corner], face_load[body_ids])
    imposed_body = (body_load(HullMesh(points, faces), panels, body_ids, arrays['pose_position'])@state[7:13])[fixed]
    components = dict(auxiliary_free=auxiliary_free[fixed],
        body_neumann_residual=auxiliary_body[fixed]-imposed_body,
        outer_neumann_residual=(auxiliary_all-auxiliary_free-auxiliary_body)[fixed]-(imposed_nonfree-imposed_body),
        calderon_consistency=(energy_action-auxiliary_all)[fixed])
    closure = sum(components.values())-reaction
    factor = cho_factor(data['free_mass'])
    densities = {name: cho_solve(factor, value) for name, value in components.items()}
    density = cho_solve(factor, reaction)
    local = np.full(len(points), -1, dtype=int)
    local[fixed] = np.arange(len(fixed))
    free_faces = faces[free_ids]
    shape = triangle_samples(np.eye(3)[None], 3)[0][0]
    # Match the actual vertical-velocity peak, not a maximum of another field.
    triangle = points[free_faces]
    delta = triangle[:, 1:, :2]-triangle[:, :1, :2]
    slope = np.linalg.solve(delta, (triangle[:, 1:, 2]-triangle[:, :1, 2])[..., None])[..., 0]
    tangent = np.linalg.solve(delta, (potential[free_faces[:, 1:]]-potential[free_faces[:, :1]])[..., None])[..., 0]
    jacobian_square = 1+np.sum(slope*slope, axis=1)
    normal = np.einsum('qc,fc->fq', shape, density[local[free_faces]])
    tangential = np.sum(slope*tangent, axis=1)
    vertical = (normal+tangential[:, None])/jacobian_square[:, None]
    face, sample = np.unravel_index(np.argmax(abs(vertical)), vertical.shape)
    nodal_weights = shape[sample]
    node_ids = local[free_faces[face]]
    free_support, body_support = np.zeros(len(points)), np.zeros(len(points))
    for corner in range(3):
        np.add.at(free_support, faces[free_ids, corner], panels.areas[free_ids]/3)
        np.add.at(body_support, faces[body_ids, corner], panels.areas[body_ids]/3)
    names = list(densities)
    density_columns = np.column_stack(list(densities.values()))
    def observation(face, nodal_weights):
        node_ids = local[free_faces[face]]
        return dict(location_world=(nodal_weights@triangle[face]).tolist(),
            reaction_vertical_m_per_s=float((nodal_weights@density[node_ids]+tangential[face])/jacobian_square[face]),
            tangential_vertical_contribution_m_per_s=float(tangential[face]/jacobian_square[face]),
            signed_vertical_contributions_m_per_s={name: float(nodal_weights@value[node_ids]/jacobian_square[face])
                                                  for name, value in densities.items()},
            sample_barycentric_weights=nodal_weights.tolist(),
            node_support=[dict(node=int(node), body_basis_support_m2=float(body_support[node]),
                free_basis_support_m2=float(free_support[node]),
                body_to_free_support_ratio=float(body_support[node]/free_support[node]))
                for node in free_faces[face]])

    result = dict(qualified=False,
        scope='Exact algebraic decomposition of saved discrete reaction; no component is declared physical truth.',
        input_sha256={str(p): sha256(p.read_bytes()).hexdigest() for p in (audit, audit.with_suffix('.npz'))},
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), order=order,
        dependency_sha256={name: report['source_sha256'][name] for name in DEPENDENCIES},
        extra_free_refinement=report.get('extra_free_refinement', 0),
        local_body_refinement=local_refinement,
        maximum_weak_closure_error=float(np.max(abs(closure))),
        recovered_graph_flux_mass_gram=dict(component_order=names,
            matrix=(density_columns.T@data['free_mass']@density_columns).tolist()),
        peak=observation(face, nodal_weights),
        components={name: dict(weak_euclidean_norm=float(np.linalg.norm(components[name])),
            recovered_mass_norm=float(np.sqrt(value@data['free_mass']@value)))
                    for name, value in densities.items()})
    if probe_point is not None:
        probe = np.asarray(probe_point, dtype=float)
        if probe.shape != (3,) or not np.isfinite(probe).all():
            raise ValueError('probe must be a finite world XYZ point')
        bary = np.linalg.solve(delta.swapaxes(1, 2), (probe[:2]-triangle[:, 0, :2])[..., None])[..., 0]
        bary = np.column_stack((1-bary.sum(axis=1), bary))
        matches = np.flatnonzero((bary.min(axis=1) >= -1e-9)
            & (abs(np.sum(bary*triangle[:, :, 2], axis=1)-probe[2]) < 1e-9))
        if not len(matches):
            raise ValueError('probe is not on the saved free surface')
        result['fixed_probe'] = dict(requested_world=probe.tolist(),
            observations=[observation(int(face), bary[face]) for face in matches])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--order', type=int, default=48)
    parser.add_argument('--probe-point', type=float, nargs=3)
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))
