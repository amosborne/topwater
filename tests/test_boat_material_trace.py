"""Material samples survive independent wet-boundary triangulations."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

from topwater import InitialState
from topwater.frames import rotation_matrix
from .reference_boat import equilibrium, reference_body
from .experiments.free_surface.boat_world_boundary import (
    construct, body_material_samples, interpolate_body_trace,
)


def affine_material_potential(mesh, body, pose):
    coordinates = np.linalg.solve(rotation_matrix(pose.orientation),
                                  (mesh.vertices-pose.position).T).T+body.cg
    return .4+coordinates@np.array([.2, -.1, .3])


def test_offset_trace_uses_material_coverage_and_preserves_gauge(tmp_path):
    body, pose = reference_body(), equilibrium()[0]
    mesh, arrays, report = construct(tmp_path/'current.json')
    assert report['complete'] and report['geometry_roundoff_checks_passed']
    points, groups = body_material_samples(arrays, np.full((1, 3), 1/3))
    phi = affine_material_potential(mesh, body, pose)
    expected = .4+points@np.array([.2, -.1, .3])
    assert_allclose(interpolate_body_trace(arrays, phi, points, groups), expected, atol=2e-14, rtol=0)
    assert_allclose(interpolate_body_trace(arrays, phi+7, points, groups), expected+7, atol=4e-14, rtol=0)
    original = arrays['source_body_vertices'][arrays['source_body_faces'][arrays['body_vertex_source_triangle']]]
    reconstructed = np.einsum('fak,fakj->faj', arrays['body_vertex_source_barycentric'], original)
    assert_allclose(reconstructed, arrays['body_material_vertices'], atol=1e-15, rtol=0)

    offset = InitialState(position=pose.position+[0, 0, .001], orientation=pose.orientation)
    moved_mesh, moved, moved_report = construct(tmp_path/'offset.json', body=body, pose=offset,
        parent_xy=arrays['parent_nodes'][:, :2], parents=arrays['parents'],
        parent_eta=arrays['parent_nodes'][:, 2], parent_psi=arrays['parent_psi'])
    assert moved_report['supplied_inputs'] == dict(body=True, pose=True, parent_xy=True,
                                                parent_eta=True, parent_psi=True)
    # Deep material points remain wet in BOTH snapshots. Matching node labels
    # would be invalid even if the independent meshes happen to have equal size.
    deep = points[..., 2] < .1
    assert deep.any()
    moved_phi = affine_material_potential(moved_mesh, body, offset)
    assert_allclose(interpolate_body_trace(moved, moved_phi, points[deep], groups[deep]),
                    expected[deep], atol=2e-14, rtol=0)
    assert np.array_equal(moved['parent_psi'], arrays['parent_psi'])
    # This original transom point is physically dry, not an extrapolation site.
    with pytest.raises(ValueError, match='lack wet'):
        interpolate_body_trace(moved, moved_phi, body.mesh.vertices[2:3], np.array([8]))

