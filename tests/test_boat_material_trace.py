"""Material samples survive independent wet-boundary triangulations."""
from math import fsum
import numpy as np
import pytest
from numpy.testing import assert_allclose

from topwater import InitialState
from topwater.frames import rotation_matrix
from .reference_boat import equilibrium, reference_body
from .experiments.free_surface.boat_world_boundary import (
    body_material_samples, interpolate_body_trace,
)
from .reference_boat import boundary_snapshot
from .experiments.free_surface.boat_world_field import refine_boundary


def affine_material_potential(points, body, pose):
    coordinates = np.linalg.solve(rotation_matrix(pose.orientation),
                                  (points-pose.position).T).T+body.cg
    return .4+coordinates@np.array([.2, -.1, .3])


def test_two_marked_edges_use_conforming_oriented_three_triangle_split():
    # Two adjacent marked faces of a closed tetrahedron force BOTH neighbours
    # to use the two-edge pattern; the actual boat regions exercise 0/1/3 edges.
    points = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    arrays = dict(points=points, faces=faces, body_faces=np.array([0, 1]),
        free_faces=np.array([3]), wall_faces=np.array([2]), bottom_faces=np.array([], dtype=int),
        body_source_group=np.array([0, 1]), free_source_parent=np.array([0]),
        body_material_vertices=points[faces[:2]], parent_nodes=points[1:],
        free_nodes=np.array([1, 2, 3]),
        free_parent_indices=np.array([[0, -1, -1], [1, -1, -1], [2, -1, -1]]),
        free_parent_weights=np.array([[1., 0, 0]]*3))
    refined = refine_boundary(arrays, region='body')
    assert np.array_equal(np.bincount(refined['refinement_parent_face']), [4, 4, 3, 3])
    child = refined['points'][refined['faces']]
    parent = points[faces][refined['refinement_parent_face']]
    normal = np.cross(parent[:, 1]-parent[:, 0], parent[:, 2]-parent[:, 0])
    child_normal = np.cross(child[:, 1]-child[:, 0], child[:, 2]-child[:, 0])
    assert np.all(np.einsum('ij,ij->i', normal, child_normal) > 0)
    # The fixed combinatorial choice is repeatable, including symmetric faces.
    assert np.array_equal(refine_boundary(arrays, region='body')['faces'], refined['faces'])
    local = refine_boundary(arrays, region='body', marked_faces=np.array([0]))
    assert np.array_equal(np.bincount(local['refinement_parent_face']), [4, 2, 2, 2])
    assert_allclose(local['refinement_moment_errors'], 0., atol=1e-16)
    with pytest.raises(ValueError, match='subset'):
        refine_boundary(arrays, region='body', marked_faces=np.array([3]))
    with pytest.raises(ValueError, match='region'):
        refine_boundary(arrays, region='walls')


def test_offset_trace_uses_material_coverage_and_preserves_gauge(tmp_path):
    body, pose = reference_body(), equilibrium()[0]
    mesh, arrays, report = boundary_snapshot(tmp_path/'current.json')
    assert report['complete'] and report['geometry_roundoff_checks_passed']
    points, groups = body_material_samples(arrays, np.full((1, 3), 1/3))
    phi = affine_material_potential(mesh.vertices, body, pose)
    expected = .4+points@np.array([.2, -.1, .3])
    hints = np.broadcast_to(np.arange(len(arrays['body_faces']))[:, None], groups.shape)
    assert_allclose(interpolate_body_trace(arrays, phi, points, groups), expected, atol=2e-14, rtol=0)
    assert_allclose(interpolate_body_trace(arrays, phi, points, groups, face_hints=hints),
                    expected, atol=2e-14, rtol=0)
    # Incorrect indices are only hints: facet/containment checks must fall back.
    assert_allclose(interpolate_body_trace(arrays, phi, points, groups,
                    face_hints=np.full(groups.shape, len(arrays['body_faces'])+1)), expected, atol=2e-14, rtol=0)
    assert_allclose(interpolate_body_trace(arrays, phi+7, points, groups), expected+7, atol=4e-14, rtol=0)
    original = arrays['source_body_vertices'][arrays['source_body_faces'][arrays['body_vertex_source_triangle']]]
    reconstructed = np.einsum('fak,fakj->faj', arrays['body_vertex_source_barycentric'], original)
    assert_allclose(reconstructed, arrays['body_material_vertices'], atol=1e-15, rtol=0)

    offset = InitialState(position=pose.position+[0, 0, .001], orientation=pose.orientation)
    moved_mesh, moved, moved_report = boundary_snapshot(tmp_path/'offset.json', body=body, pose=offset,
        parent_xy=arrays['parent_nodes'][:, :2], parents=arrays['parents'],
        parent_eta=arrays['parent_nodes'][:, 2], parent_psi=arrays['parent_psi'])
    assert moved_report['input_sha256'] == report['input_sha256']
    np.testing.assert_array_equal(moved['pose_position'], offset.position)
    np.testing.assert_array_equal(moved['parent_nodes'], arrays['parent_nodes'])
    # Deep material points remain wet in BOTH snapshots. Matching node labels
    # would be invalid even if the independent meshes happen to have equal size.
    deep = points[..., 2] < .1
    assert deep.any()
    moved_phi = affine_material_potential(moved_mesh.vertices, body, offset)
    assert_allclose(interpolate_body_trace(moved, moved_phi, points[deep], groups[deep]),
                    expected[deep], atol=2e-14, rtol=0)
    assert_allclose(interpolate_body_trace(moved, moved_phi, points[deep], groups[deep],
                    face_hints=hints[deep]), expected[deep], atol=2e-14, rtol=0)
    assert np.array_equal(moved['parent_psi'], arrays['parent_psi'])
    # This original transom point is physically dry, not an extrapolation site.
    with pytest.raises(ValueError, match='lack wet'):
        interpolate_body_trace(moved, moved_phi, body.mesh.vertices[2:3], np.array([8]))
    with pytest.raises(ValueError, match='lack wet'):
        interpolate_body_trace(moved, moved_phi, body.mesh.vertices[2:3], np.array([8]),
                               face_hints=np.array([0]))


@pytest.mark.parametrize('region', ['all', 'body', 'free'])
def test_refinement_preserves_physical_moments_and_affine_material_world_traces(tmp_path, region):
    mesh, arrays, _ = boundary_snapshot(tmp_path/'coarse.json')
    refined = refine_boundary(arrays, region=region)
    selected = (mesh.faces if region == 'all' else mesh.faces[arrays[region+'_faces']])
    marked = np.unique(np.sort(np.concatenate([selected[:, [0, 1]],
        selected[:, [1, 2]], selected[:, [2, 0]]]), axis=1), axis=0)
    assert len(refined['faces']) == len(mesh.faces)+2*len(marked)
    assert refined['refinement_marked_edges'] == len(marked)
    counts = np.bincount(refined['refinement_parent_face'], minlength=len(mesh.faces))
    assert np.all(counts >= 1)
    if region == 'all':
        assert np.all(counts == 4)
        # Preserve the original red-refinement numbering, not just its geometry.
        edges, inverse = np.unique(np.sort(np.concatenate([mesh.faces[:, [0, 1]],
            mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]]]), axis=1),
            axis=0, return_inverse=True)
        ab, bc, ca = (inverse+len(mesh.vertices)).reshape(3, -1)
        a, b, c = mesh.faces.T
        expected = np.vstack([np.column_stack(triangle) for triangle in
            ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))])
        assert np.array_equal(refined['faces'], expected)
        assert np.array_equal(refine_boundary(arrays)['faces'], expected)
    else:
        assert np.all(counts[arrays[region+'_faces']] == 4)
        assert len(refined['faces']) < 4*len(mesh.faces)
    # Each child has its parent's orientation and the complete child areas sum
    # to its area. HullMesh inside refine_boundary independently checks closure.
    coarse = mesh.vertices[mesh.faces]
    normal = np.cross(coarse[:, 1]-coarse[:, 0], coarse[:, 2]-coarse[:, 0])
    children = refined['points'][refined['faces']]
    child_normal = np.cross(children[:, 1]-children[:, 0], children[:, 2]-children[:, 0])
    owner = refined['refinement_parent_face']
    assert np.all(np.einsum('ij,ij->i', child_normal, normal[owner]) > 0)
    area_sum = np.bincount(owner, weights=np.linalg.norm(child_normal, axis=1),
                           minlength=len(mesh.faces))
    assert_allclose(area_sum, np.linalg.norm(normal, axis=1), rtol=3e-14, atol=2e-15)
    np.savez(tmp_path/'refined.npz', **refined)
    with np.load(tmp_path/'refined.npz', allow_pickle=False) as archive:
        for key in archive.files:
            assert archive[key].dtype.kind != 'O'
    assert np.array_equal(refined['source_body_faces'], arrays['source_body_faces'])
    assert 'body_vertex_source_barycentric' not in refined

    def moments(data, ids):
        t = data['points'][data['faces'][ids]]
        normal = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
        area = np.linalg.norm(normal, axis=1)/2
        first = area[:, None]*t.mean(axis=1)
        return np.array([fsum(area), *(fsum(column) for column in first.T)])

    for group in np.unique(arrays['body_source_group']):
        before = arrays['body_faces'][arrays['body_source_group'] == group]
        after = refined['body_faces'][refined['body_source_group'] == group]
        assert_allclose(moments(refined, after), moments(arrays, before), rtol=2e-14, atol=3e-14)
    for kind in ('free', 'wall', 'bottom'):
        assert_allclose(moments(refined, refined[kind+'_faces']),
                        moments(arrays, arrays[kind+'_faces']), rtol=2e-14, atol=3e-14)
    # An arbitrary affine WORLD field must survive interpolation on subdivided
    # faces, including shared waterline nodes, independently of mesh numbering.
    parent = arrays['parent_nodes']
    values = parent[:, :2]@np.array([.2, -.4])+.7
    ids, weights = refined['free_parent_indices'], refined['free_parent_weights']
    trace = (weights*values[np.maximum(ids, 0)]).sum(axis=1)
    free = refined['free_nodes']
    assert_allclose(trace, refined['points'][free, :2]@np.array([.2, -.4])+.7, atol=2e-15)
    assert_allclose((weights*parent[np.maximum(ids, 0), 2]).sum(axis=1),
                    refined['points'][free, 2], atol=2e-15)

    body, pose = reference_body(), equilibrium()[0]
    phi = affine_material_potential(refined['points'], body, pose)
    samples, groups = body_material_samples(arrays, np.full((1, 3), 1/3))
    expected = .4+samples@np.array([.2, -.1, .3])
    assert_allclose(interpolate_body_trace(refined, phi, samples, groups), expected, atol=2e-14)

