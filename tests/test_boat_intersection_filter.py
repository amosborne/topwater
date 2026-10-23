"""Rejection-only filtering must preserve the rational cold geometry owner."""
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest

from topwater import InitialState
from tests.reference_boat import equilibrium
from .experiments.free_surface.boat_intersection_filter import IntersectionRejectionFilter
from .experiments.free_surface.boat_parent_clip_proof import PlanePoints, cross, dot
from .experiments.free_surface.boat_world_boundary import construct
from .experiments.free_surface.boat_geometry_template import evaluate_template


def _stats():
    return dict(calls=0, candidates=0, interval_constructions=0, rejected=0, exact_fallback=0)


def _exact_outside(planes, triple, constraints):
    a, b, c = [planes[i] for i in triple]
    determinant = dot(a[:3], cross(b[:3], c[:3]))
    if determinant == 0:
        return False
    numerator = tuple(sum(d[3]*v[k] for d, v in zip((a, b, c),
        (cross(b[:3], c[:3]), cross(c[:3], a[:3])))) for k in range(3))
    return any(sign*(dot(planes[i][:3], numerator)-planes[i][3]*determinant)*determinant > 0
               for i, sign in constraints)


def test_filtered_rejections_match_exact_indirect_signs():
    rng = np.random.default_rng(23)
    # Include exact zeros, dependent triples, non-dyadic coefficients and large
    # cancelling terms; no random floating construction determines the oracle.
    raw = rng.integers(-7, 8, size=(14, 4))
    planes = [tuple(F(int(x), 3+row % 4) for x in entries)
              for row, entries in enumerate(raw)]
    planes.extend([(F(1), F(0), F(0), F(0)),
                   (F(1), F(0), F(0), F(1, 10**20)),
                   (F(10**20), F(10**20+1), F(1), F(1))])
    triples = list(combinations(range(len(planes)), 3))
    constraints = [(i, 1 if i % 2 else -1) for i in range(len(planes))]
    stats = _stats()
    filtered = IntersectionRejectionFilter(planes, stats)
    rejected = filtered.reject(triples, constraints)
    assert rejected.any() and (~rejected).any()
    for triple, skip in zip(triples, rejected):
        if skip:
            assert _exact_outside(planes, triple, constraints)
    assert stats['interval_constructions'] == len(triples)
    np.testing.assert_array_equal(filtered.reject(triples, constraints), rejected)
    assert stats['interval_constructions'] == len(triples)
    assert stats['rejected'] == 2*int(rejected.sum())


def test_uncertain_degenerate_and_nonfinite_cases_fall_through():
    planes = [(F(1), F(0), F(0), F(0)), (F(0), F(1), F(0), F(0)),
              (F(0), F(0), F(1), F(0)), (F(1), F(0), F(0), F(0))]
    filtered = IntersectionRejectionFilter(planes, _stats())
    assert not filtered.reject([(0, 1, 2), (0, 1, 3)], [(3, 1)]).any()
    planes[3] = (F(10**400), F(0), F(0), F(1))
    filtered = IntersectionRejectionFilter(planes, _stats())
    assert not filtered.reject([(0, 1, 2)], [(3, 1)]).any()


@pytest.mark.parametrize('normal_scale', [F(1), F(1, 2**540)])
def test_rational_underflow_never_causes_false_rejection(normal_scale):
    # Below the smallest float64 subnormal: coefficient conversion encloses it
    # by zero/minsubnormal rather than silently declaring an exact zero.
    tiny = F(1, 2**1100)
    planes = [(normal_scale, F(0), F(0), F(0)),
              (F(0), normal_scale, F(0), F(0)),
              (F(0), F(0), normal_scale, F(0)),
              (F(1), F(0), F(0), tiny),
              (F(1), F(0), F(0), -tiny)]
    # With the small normal scale, cross products and det(A) also underflow,
    # despite the three defining normals being exactly independent.
    triples = list(combinations(range(len(planes)), 3))
    for constraints in ([(3, 1)], [(3, -1)], [(4, 1)], [(4, -1)]):
        filtered = IntersectionRejectionFilter(planes, _stats())
        rejected = filtered.reject(triples, constraints)
        for triple, skip in zip(triples, rejected):
            if skip:
                assert _exact_outside(planes, triple, constraints)
        # The origin's classification cannot be certified at this scale.
        assert not rejected[triples.index((0, 1, 2))]
        expected_outside = constraints in ([(3, -1)], [(4, 1)])
        assert _exact_outside(planes, (0, 1, 2), constraints) == expected_outside
    assert tiny != 0 and float(tiny) == 0
    if normal_scale != 1:
        assert normal_scale**3 > 0 and float(normal_scale**2) == 0


@pytest.mark.parametrize('sign', [0, 2, -2, .5, np.nan])
def test_invalid_constraint_signs_reject(sign):
    planes = [(F(1), F(0), F(0), F(0))]*3
    filtered = IntersectionRejectionFilter(planes, _stats())
    with pytest.raises(ValueError, match='constraint signs'):
        filtered.reject([(0, 1, 2)], [(0, sign)])


def test_interval_cache_cap_clears_and_oversize_batch_falls_back(monkeypatch):
    planes = [(F(1), F(0), F(0), F(0)), (F(0), F(1), F(0), F(0)),
              (F(0), F(0), F(1), F(0)), (F(1), F(1), F(1), F(1))]
    stats = _stats()
    filtered = IntersectionRejectionFilter(planes, stats, max_entries=2)
    constraints = [(3, 1)]
    filtered.reject([(0, 1, 2), (0, 1, 3)], constraints)
    filtered.reject([(0, 2, 3), (1, 2, 3)], constraints)
    assert stats['cache_clears'] == 1
    assert stats['cache_entries'] == stats['cache_peak_entries'] == 2
    assert stats['cache_numeric_bytes'] == stats['cache_peak_numeric_bytes'] == 128
    assert len(filtered.homogeneous) <= 2
    def forbidden(*args):
        pytest.fail('oversize batch allocated homogeneous intervals')
    monkeypatch.setattr(filtered, '_construct_intervals', forbidden)
    triples = list(combinations(range(len(planes)), 3))
    rejected = filtered.reject(triples, constraints)
    assert not rejected.any()
    assert stats['oversize_batches'] == 1
    assert stats['cache_entries'] == 2


def test_exact_registry_points_incidence_and_relevant_constraints_are_unchanged():
    results = []
    definitions = [[(0, 0, 0), (0, 1, 0), (0, 0, 1)],
                   [(1, 0, 0), (1, 1, 0), (1, 0, 1)],
                   [(0, 0, 0), (1, 0, 0), (0, 0, 1)],
                   [(0, 1, 0), (1, 1, 0), (0, 1, 1)],
                   [(0, 0, 0), (1, 0, 0), (0, 1, 0)],
                   [(0, 0, 1), (1, 0, 1), (0, 1, 1)],
                   [(2, 0, 0), (0, 2, 0), (0, 0, 2)]]
    for enabled in (False, True):
        registry = PlanePoints(capture_constraints=True, cache_signs=True,
                               filter_intersections=enabled)
        planes = [registry.plane(points, ('test', i))[0] for i, points in enumerate(definitions)]
        constraints = list(zip(planes, (-1, 1, -1, 1, -1, 1, 1)))
        vertices = registry.polytope(constraints)
        results.append((registry, vertices))
    reference, actual = results
    assert reference[1] == actual[1]
    for attribute in ('planes', 'points', 'incidence', 'relevant', 'vertex_by_owner'):
        assert getattr(reference[0], attribute) == getattr(actual[0], attribute)
    assert actual[0].intersection_filter_stats['rejected'] > 0
    assert len(actual[0].cache) < len(reference[0].cache)


@pytest.mark.parametrize('generic', [False, True])
def test_canonical_cold_arrays_and_chart_replay_unchanged(tmp_path, generic):
    pose = equilibrium()[0]
    if generic:
        pose = InitialState(np.array([.123, .087, pose.position[2]]), pose.orientation)
    results = [construct(tmp_path/f'cold-{enabled}.json', pose=pose,
        export_template=True, filter_intersections=enabled) for enabled in (False, True)]
    _, reference, _ = results[0]
    _, actual, report = results[1]
    assert reference.keys() == actual.keys()
    for key in reference:
        np.testing.assert_array_equal(actual[key], reference[key], err_msg=key)
    for arrays in (reference, actual):
        replay, _ = evaluate_template(arrays, pose=pose,
            parent_eta=arrays['parent_nodes'][:, 2], parent_psi=arrays['parent_psi'])
        np.testing.assert_allclose(replay['points'], arrays['points'], atol=1e-12, rtol=0)
    stats = report['intersection_filter']
    assert stats['rejected'] > 0 and stats['exact_fallback'] > 0
    assert stats['rejected']+stats['exact_fallback'] == stats['candidates']
