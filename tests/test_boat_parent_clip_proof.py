"""Actual chine/vertical-transom incidence on two fixed world parents."""
import numpy as np
import pytest

from .experiments.free_surface.boat_parent_clip_proof import build


@pytest.mark.parametrize('heights,expected_groups', [
    ((0., 0., 0., 0.), [0, 8]),
    ((-.01, .005, .025, 0.), [0, 2, 8]),
])
def test_two_parents_share_cutpoints_and_close_across_actual_chine(heights, expected_groups):
    mesh, result = build(heights)
    assert result['finite_transom_alias']
    assert result['duplicate_coincident_body_faces_removed'] == 1
    assert result['shared_roof_edge_owned_points'] >= 2
    assert result['roof_incidence_exact']
    assert result['emitted_physical_groups'] == expected_groups
    # Independently count actual triangles: every edge twice, opposite winding.
    edges = np.concatenate([mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]]])
    _, inverse, counts = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)
    assert np.all(counts == 2)
    winding = np.bincount(inverse, weights=np.where(edges[:, 0] < edges[:, 1], 1, -1))
    assert np.all(winding == 0)
    assert abs(result['volume_error']) < 1e-12
    assert result['first_moment_max_error'] < 2e-12
    assert abs(result['affine_rigid_work_error']) < 1e-12
