"""Integer homogeneous predicates must preserve the rational geometry owner."""
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest

from topwater import InitialState
from .reference_boat import equilibrium
from .experiments.free_surface import boat_parent_clip_proof as proof
from .experiments.free_surface import boat_world_boundary as boundary
from .experiments.free_surface.boat_integer_intersections import IntegerIntersections


def _registry(planes, enabled, **limits):
    registry = proof.PlanePoints(integer_predicates=enabled, cache_signs=True,
                                 capture_constraints=True, **limits)
    registry.planes = [tuple(F(value) for value in row) for row in planes]
    registry.labels = [[] for _ in planes]
    return registry


def _ownership(registry):
    return {name: getattr(registry, name) for name in
        ('points', 'incidence', 'relevant', 'vertex_by_owner', 'registered_intersections', 'cache', 'sign_cache')}


def test_nondyadic_large_and_degenerate_intersections_and_signs_are_exact():
    planes = [(F(1, 3), F(2, 7), 0, F(1, 11)),
              (0, F(-3, 5), F(4, 9), F(2, 13)),
              (F(5, 17), 0, F(7, 19), F(-3, 23)),
              (F(10**50), F(10**50+1), 1, F(1, 10**30)),
              (1, 0, 0, 0), (1, 0, 0, F(1, 2**1100)),
              (0, 1, 0, 0), (0, 0, 1, 0), (2, 0, 0, 0)]
    old, new = [_registry(planes, mode) for mode in (False, True)]
    for triple in combinations(range(len(planes)), 3):
        a, b = old.intersection(triple), new.intersection(triple)
        assert a == b
        if a is not None:
            for index in range(len(planes)):
                assert new._sign(triple, b, index) == old._sign(triple, a, index)
    assert new.cache == old.cache
    assert new.sign_cache == old.sign_cache
    assert new.sign_requests == old.sign_requests
    assert new.sign_cache_hits == old.sign_cache_hits
    assert old.integer_geometry is None
    assert len(new.integer_geometry.homogeneous) == len(new.cache)
    # An external point must not silently inherit another triple's coordinates.
    external = (F(1, 3), F(-2, 7), F(5, 11))
    for index in range(len(planes)):
        value = old.value(index, external)
        assert new.integer_geometry.sign_point(external, index) == (value > 0)-(value < 0)


def test_polytope_registration_incidence_and_relevance_are_unchanged():
    planes = [(1, 0, 0, 0), (1, 0, 0, 1), (0, 1, 0, 0),
              (0, 1, 0, 1), (0, 0, 1, 0), (0, 0, 1, 1),
              (1, 1, 1, F(3, 2)), (1, 1, 1, 0)]
    old, new = [_registry(planes, mode) for mode in (False, True)]
    cube = [(0, -1), (1, 1), (2, -1), (3, 1), (4, -1), (5, 1)]
    for constraints in (cube, [*cube, (6, 1)], [*cube, (7, -1)], [*cube, (7, 1)]):
        assert old.polytope(constraints) == new.polytope(constraints)
        assert _ownership(new) == _ownership(old)
        assert list(new.registered_intersections) == list(old.registered_intersections)


@pytest.mark.parametrize('enabled', [False, True])
def test_existing_intersection_cap_and_plane_freeze_are_preserved(enabled):
    registry = _registry([(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (1, 1, 1, 1)],
                         enabled, max_intersections=1)
    assert registry.intersection((0, 1, 2)) == (F(0), F(0), F(0))
    assert registry.intersection((2, 1, 0)) == (F(0), F(0), F(0))
    with pytest.raises(RuntimeError, match='intersection cap'):
        registry.intersection((0, 1, 3))
    with pytest.raises(RuntimeError, match='freeze every plane'):
        registry.plane([(0, 0, 0), (1, 0, 0), (0, 1, 0)], ('late', 0))


def test_integer_owner_rejects_changed_frozen_coefficients():
    planes = [(F(1), F(0), F(0), F(0)), (F(0), F(1), F(0), F(0)),
              (F(0), F(0), F(1), F(0))]
    owner = IntegerIntersections(planes)
    owner.check_planes(planes)
    planes[0] = (F(1), F(0), F(0), F(1))
    with pytest.raises(ValueError, match='unchanged frozen'):
        owner.check_planes(planes)


def test_dense_witness_signs_include_exact_fan_points_without_sign_cache_changes():
    planes = [(1, 0, 0, F(1, 3)), (0, 1, 0, F(-2, 7)), (0, 0, 1, 0),
              (F(10**50), F(10**50+1), 1, F(1, 2**1100))]
    points = [(F(0), F(0), F(0)), (F(1, 3), F(-2, 7), F(0)),
              (F(2, 11), F(4, 13), F(-3, 17))]
    points.append(tuple(sum(point[index] for point in points)/3 for index in range(3)))
    old, new = [_registry(planes, mode) for mode in (False, True)]
    expected = np.array([[int(np.sign(old.value(index, point))) for index in range(len(planes))]
                         for point in points], dtype=np.int8)
    for registry in (old, new):
        np.testing.assert_array_equal(registry.exact_signs(points), expected)
        assert registry.exact_signs(points).dtype == np.int8
        assert registry.exact_signs([]).shape == (0, len(planes))
        assert registry.sign_cache == {}
        assert registry.sign_requests == registry.sign_cache_hits == 0


@pytest.mark.parametrize('heights', [(0.,)*4, (-.01, .005, .025, 0.)])
def test_actual_adjacent_transom_chine_geometry_is_bit_identical(monkeypatch, heights):
    old_mesh, old_report = proof.build(heights)
    original = proof.PlanePoints
    with monkeypatch.context() as context:
        context.setattr(proof, 'PlanePoints', lambda **options: original(**options, integer_predicates=True))
        new_mesh, new_report = proof.build(heights)
    np.testing.assert_array_equal(new_mesh.vertices, old_mesh.vertices)
    np.testing.assert_array_equal(new_mesh.faces, old_mesh.faces)
    for key in old_report.keys()-{'seconds'}:
        assert new_report[key] == old_report[key], key


@pytest.mark.parametrize('generic', [False, True])
def test_complete_cold_arrays_and_template_are_bit_identical(monkeypatch, tmp_path, generic):
    pose = equilibrium()[0]
    if generic:
        pose = InitialState(pose.position+np.array([.123, .087, 0.]), pose.orientation)
    _, old, _ = boundary.construct(tmp_path/'old.json', pose=pose, export_template=True)
    original = boundary.PlanePoints
    captured = []
    def factory(**options):
        registry = original(**options, integer_predicates=True)
        captured.append(registry)
        return registry
    with monkeypatch.context() as context:
        context.setattr(boundary, 'PlanePoints', factory)
        _, new, _ = boundary.construct(tmp_path/'new.json', pose=pose, export_template=True)
    assert old.keys() == new.keys()
    for key in old:
        np.testing.assert_array_equal(new[key], old[key], err_msg=key)
    registry = captured[0]
    assert len(registry.integer_geometry.homogeneous) == len(registry.cache)
