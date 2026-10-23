"""Actual chine/vertical-transom incidence on two fixed world parents."""
import numpy as np
import pytest
from fractions import Fraction

from .experiments.free_surface.boat_parent_clip_proof import (
    build, hull_planes, PlanePoints, polygon_centroid, merge_body_polygons,
)
from .reference_boat import reference_body, GROUPS


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


def test_authored_coplanar_group_rejects_a_nonplanar_constituent():
    body = reference_body()
    points = np.array([[Fraction(float(x)) for x in p] for p in body.mesh.vertices], dtype=object)
    # Vertex5 belongs to the second triangle of the aft V-plane, but not the
    # first triangle used to define that plane. Convexity alone is insufficient.
    points[5, 2] += Fraction(1, 1000)
    with pytest.raises(ValueError, match='noncoplanar facet group'):
        hull_planes(PlanePoints(), points, body.mesh.faces, groups=GROUPS)


def test_fan_center_is_invariant_to_collinear_vertex_insertion_and_has_continuous_cut_limit():
    f = Fraction
    points = [(f(x), f(y), f(2*x+3*y)) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
    center = polygon_centroid(points, [0, 1, 2, 3])
    assert center == (f(1, 2), f(1, 2), f(5, 2))
    epsilon = f(1, 10**6)
    points.append((epsilon, f(0), 2*epsilon))
    assert polygon_centroid(points, [0, 4, 1, 2, 3]) == center
    points.append((f(0), epsilon, 3*epsilon))
    # Cutting off a tiny corner changes area and centroid by O(epsilon²),
    # whereas the old vertex mean tends to (.4,.4), an O(1) displacement.
    cut = polygon_centroid(points, [4, 1, 2, 3, 5])
    assert max(abs(a-b) for a, b in zip(cut, center)) < 3*epsilon**2


def test_body_union_removes_artificial_seam_but_preserves_shared_nonbody_cutpoint():
    registry = PlanePoints()
    registry.points = [(Fraction(x), Fraction(y), Fraction(0))
                       for x, y in [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (1, -1)]]
    registry.labels = [[('hull', (0,))]]
    registry.incidence = [frozenset({0}) for _ in registry.points]
    polygons = [([0, 1, 4, 3], 0, 1, 'body', 0), ([1, 2, 5, 4], 0, 1, 'body', 0),
                ([0, 6, 1], 1, 1, 'column', 0)]
    merged, report = merge_body_polygons(registry, polygons)
    body = [p for p, _, _, kind, _ in merged if kind == 'body']
    assert len(body) == 1
    assert set(body[0]) == {0, 1, 2, 3, 5}  # Required1 retained; artificial4 removed.
    assert report['body_polygons_before'] == 2
    assert report['body_polygons_after'] == 1
    assert report['group_area_and_first_moment_preserved_exactly']
    # A disjoint wet region must fail explicitly, never be bridged by convex hull.
    registry.points[2] = (Fraction(4), Fraction(0), Fraction(0))
    registry.points[5] = (Fraction(4), Fraction(1), Fraction(0))
    polygons[1] = ([2, 5, 4], 0, 1, 'body', 0)
    with pytest.raises(ValueError, match='multiple|branching'):
        merge_body_polygons(registry, polygons[:2])
