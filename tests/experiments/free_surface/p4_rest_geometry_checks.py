"""Actual P4 canoe geometry-only prerequisites for an exact rest linear limit.

No CURRENT fluid stiffness, pressure factorization or PDE solve is performed.
The one whole-geometry build owns fixed reference ALE factors. Raw directional
errors are saved without claiming linear/nonlinear or Capytaine qualification.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from topwater.hydrostatics import hydrostatics
from .canonical_canoe import body_map
from .p4_whole_geometry import P4WholeCanoeGeometry
from .p4_canoe_field import surface_maps, body_load
from .p4_canoe import loads as affine_loads
from .p4_tetra import face_nodes, indices


def _maximum(value):
    return float(np.max(abs(value))) if np.size(value) else 0.


def _state(owner, q):
    orientation = (Rotation.from_rotvec(q[3:6])*Rotation.from_quat(
        owner.domain.rest.orientation, scalar_first=True)).as_quat(scalar_first=True)
    return owner.evaluate(q[:3], orientation, q[6:])


def _gravity_gradient(owner, q, state=None):
    state = _state(owner, q) if state is None else state
    moments = owner.boundary_moments(state, order=9)
    world_body, height = state.pullback(moments['gravity_gradient'])
    gradient = np.r_[body_map(q[3:6]).T@world_body, height]
    gradient[2] += owner.domain.body.mass*9.81
    return gradient, state, moments


def _selected_modes(owner):
    nodes = face_nodes(owner.mesh, owner.free_faces)
    near = owner.chart_to_global[owner.chart.free_nodes]
    candidates = np.flatnonzero(np.all(~np.isin(nodes, near), axis=1))
    if not len(candidates):
        raise ValueError('no independent far P4 surface triangle found')
    face = int(candidates[len(candidates)//2])
    far = np.zeros(len(owner.free_nodes))
    far[np.searchsorted(owner.free_nodes, nodes[face])] = 27*np.prod(indices(3)/4, axis=1)
    local = int(owner.chart.high_contact_nodes[len(owner.chart.high_contact_nodes)//2])
    global_node = int(owner.chart_to_global[local])
    contact = np.zeros(len(owner.free_nodes))
    contact[np.searchsorted(owner.free_nodes, global_node)] = 1.
    return [dict(name='far interior P4 bubble', vector=far, steps=(2e-4, 1e-4),
                 source=dict(free_face_index=face, global_trace_nodes=nodes[face].tolist())),
            dict(name='independent high contact height', vector=contact, steps=(2e-5, 1e-5),
                 source=dict(chart_node=local, global_node=global_node))]


def _revision():
    try:
        return subprocess.check_output(['git', 'rev-parse', 'HEAD'],
            cwd=Path(__file__).resolve().parents[3], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def run():
    owner = P4WholeCanoeGeometry.make()
    q = np.r_[owner.domain.rest.position, np.zeros(3+len(owner.free_nodes))]
    gradient, state, moments = _gravity_gradient(owner, q)
    eta_map, transport, mass = surface_maps(owner, state, order=9)
    eta_error = (eta_map-mass).tocoo()
    current_load = body_load(owner.mesh, owner.body_faces, state.points, q[:3], order=9)
    # Separate affine physical normals and degree5 analytic load quadrature,
    # rather than using the current curved-field routine as its own reference.
    reference_load, _ = affine_loads(owner.mesh, owner.body_faces, q[:3])
    load_error = current_load-reference_load
    expected_body = hydrostatics(owner.domain.body, owner.domain.rest).stiffness
    modes = _selected_modes(owner)
    vectors = np.array([mode['vector'] for mode in modes])
    arrays = dict(rest_gravity_gradient=gradient, hydrostatic_body_hessian=expected_body,
        eta_mass_error_row=eta_error.row, eta_mass_error_column=eta_error.col,
        eta_mass_error_value=eta_error.data, body_transport=transport,
        body_load_error=load_error, selected_height_modes=vectors,
        reference_mass_data=mass.data, reference_mass_indices=mass.indices,
        reference_mass_indptr=mass.indptr)
    # Physical closed-volume gradient = one-form tested with constant potential.
    world_body, volume_height = state.pullback(moments['volume_gradient'])
    volume_body = body_map(q[3:6]).T@world_body
    expected_volume_body = body_map(q[3:6]).T@(current_load.sum(axis=0)+transport.sum(axis=0))
    expected_volume_height = np.asarray(eta_map.T@np.ones(len(owner.free_nodes))).ravel()
    arrays.update(volume_body_gradient=volume_body,
        volume_body_one_form_error=volume_body-expected_volume_body,
        volume_height_gradient=volume_height,
        volume_height_one_form_error=volume_height-expected_volume_height)
    body_levels = []
    # Mixed units are handled by a declared length normalization, while raw
    # dimensional matrices and every selected directional action are retained.
    length = 3.5
    dimension = np.array([1., 1., 1., length, length, length])
    hessian_scale = 1000*9.81*length**2*dimension[:, None]*dimension[None, :]
    for level, step in enumerate((2e-4, 1e-4)):
        body_hessian = np.empty((6, 6))
        selected_cross = np.empty((len(modes), 6))
        height_cross_max = []
        for axis in range(6):
            direction = np.zeros_like(q); direction[axis] = step
            plus = _gravity_gradient(owner, q+direction)[0]
            minus = _gravity_gradient(owner, q-direction)[0]
            derivative = (plus-minus)/(2*step)
            body_hessian[:, axis] = derivative[:6]
            selected_cross[:, axis] = vectors@derivative[6:]
            height_cross_max.append(_maximum(derivative[6:]))
        error = body_hessian-expected_body
        arrays[f'body_hessian_level{level}'] = body_hessian
        arrays[f'body_hessian_error_level{level}'] = error
        arrays[f'height_body_selected_cross_level{level}'] = selected_cross
        record = dict(step=step, numeric_body_hessian=body_hessian.tolist(),
            body_hessian_error=error.tolist(), maximum_absolute_error=float(np.max(abs(error))),
            maximum_nondimensional_error=float(np.max(abs(error)/hessian_scale)),
            maximum_hessian_asymmetry=_maximum(body_hessian-body_hessian.T),
            height_body_cross_maximum_by_body_column=height_cross_max,
            selected_height_body_cross=selected_cross.tolist())
        body_levels.append(record)
        print(json.dumps(dict(progress=True, stage='geometry-only body Hessian complete', **record)), flush=True)
    selected_results = []
    for index, mode in enumerate(modes):
        vector = mode['vector']
        expected_height = 9810*(mass@vector)
        arrays[f'height_hessian_expected_mode{index}'] = expected_height
        levels = []
        for level, step in enumerate(mode['steps']):
            direction = np.r_[np.zeros(6), vector]
            plus = _gravity_gradient(owner, q+step*direction)[0]
            minus = _gravity_gradient(owner, q-step*direction)[0]
            derivative = (plus-minus)/(2*step)
            error = derivative[6:]-expected_height
            cross_reference = np.asarray(body_levels[level]['selected_height_body_cross'])[index]
            arrays[f'height_hessian_action_mode{index}_level{level}'] = derivative[6:]
            arrays[f'height_hessian_error_mode{index}_level{level}'] = error
            arrays[f'body_height_cross_mode{index}_level{level}'] = derivative[:6]
            levels.append(dict(step=step, maximum_absolute_height_action_error=_maximum(error),
                relative_height_action_l2_error=float(np.linalg.norm(error)/np.linalg.norm(expected_height)),
                body_height_cross=derivative[:6].tolist(),
                body_height_cross_maximum=_maximum(derivative[:6]),
                cross_symmetry_error=(derivative[:6]-cross_reference).tolist()))
        selected_results.append(dict(name=mode['name'], source=mode['source'], differences=levels,
            fine_over_coarse_height_error=levels[1]['maximum_absolute_height_action_error']
                /levels[0]['maximum_absolute_height_action_error'] if levels[0]['maximum_absolute_height_action_error'] else None))
        print(json.dumps(dict(progress=True, stage='geometry-only height Hessian complete',
                              result=selected_results[-1])), flush=True)
    norm = np.linalg.norm(reference_load, axis=0)
    here = Path(__file__).resolve().parent
    source_names = ('p4_rest_geometry_checks.py', 'p4_whole_geometry.py', 'p4_volume_map.py',
                    'p4_canoe_contact.py', 'canoe_contact.py', 'p4_canoe_field.py',
                    'p4_canoe.py', 'p4_surface.py')
    report = dict(qualified=False, case='actual canonical P4 canoe geometry-only rest identities',
        current_fluid_operator_assembled=False, fluid_PDE_solved=False,
        reference_geometry_factors_built_once=True, revision=_revision(),
        source_sha256={name: sha256((here/name).read_bytes()).hexdigest() for name in source_names},
        independent_hydrostatics_source_sha256=sha256(Path(hydrostatics.__code__.co_filename).read_bytes()).hexdigest(),
        reference_points_sha256=sha256(np.asarray(owner.mesh.points, dtype='<f8').tobytes()).hexdigest(),
        elements_sha256=sha256(np.asarray(owner.mesh.elements, dtype='<i8').tobytes()).hexdigest(),
        free_nodes_sha256=sha256(np.asarray(owner.free_nodes, dtype='<i8').tobytes()).hexdigest(),
        active_nodes=len(np.unique(owner.mesh.elements)), eta_modes=len(owner.free_nodes),
        rest_body_configuration=q[:6].tolist(), rest_reference_point_error=_maximum(state.points-owner.mesh.points),
        rest_gravity_gradient_body=gradient[:6].tolist(), rest_gravity_gradient_height_maximum=_maximum(gradient[6:]),
        eta_map_minus_mass_maximum=_maximum(eta_error.data), body_transport_maximum=_maximum(transport),
        body_load_maximum_error_by_mode=np.max(abs(load_error), axis=0).tolist(),
        body_load_relative_l2_error_by_mode=[float(np.linalg.norm(load_error[:, i])/value) if value else None for i, value in enumerate(norm)],
        free_surface_body_load_maximum_error=_maximum(load_error[owner.free_nodes]),
        volume_body_one_form_error=(volume_body-expected_volume_body).tolist(),
        volume_height_one_form_maximum_error=_maximum(volume_height-expected_volume_height),
        independent_expected_body_hessian=expected_body.tolist(), body_differences=body_levels,
        body_hessian_error_fine_over_coarse=body_levels[1]['maximum_nondimensional_error']
            /body_levels[0]['maximum_nondimensional_error'] if body_levels[0]['maximum_nondimensional_error'] else None,
        selected_height_differences=selected_results,
        body_hessian_normalization='rho*g*L^2*d_i*d_j, d=(1,1,1,L,L,L), L=3.5m; raw dimensional matrices also retained.',
        scope='Geometry/surface action/gravity prerequisites only. Current closed-boundary gradients '
            'are pulled back to canonical body coordinates with K(r).T; bodymgz added once. '
            'Expected hydrostatic stiffness comes from independently clipped flat-water physical hull. '
            'Selected height actions and cross blocks do not establish a complete Hessian or six-mode '
            'dynamic first variation. No current fluid PDE, trajectory, spatial, failure-event, '
            'Capytaine or performance qualification. Physical boundary guards apply; no new all-cell '
            'Jacobian scan is implied by these boundary-only moment checks.')
    return report, arrays


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-rest-geometry-checks.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        report, arrays = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output.with_suffix('.npz'), **arrays)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
