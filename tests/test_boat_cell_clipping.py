"""Exact cell replacement controls: geometry AND registration ownership."""
from fractions import Fraction as F

import numpy as np
import pytest

from topwater import InitialState
from .reference_boat import equilibrium
from .experiments.free_surface import boat_parent_clip_proof as reference
from .experiments.free_surface import boat_world_boundary as boundary
from .experiments.free_surface import boat_cell_clipping as candidate


def _ownership(registry):
    return {key: getattr(registry, key) for key in
            ('planes', 'points', 'incidence', 'relevant', 'vertex_by_owner', 'registered_intersections')}


@pytest.mark.parametrize('cut,lower_dimension', [
    ((1, 0, 0, 2), False),       # No cut.
    ((1, 0, 0, F(1, 3)), False),
    ((1, 0, 0, 0), True),        # Planar tangent set.
    ((1, 1, 0, 0), True),        # Tangent vertical edge.
    ((1, 1, 1, 0), True),        # Tangent point.
    ((1, 0, 0, -1), False),      # Strictly empty.
])
def test_exact_clip_preserves_tangent_sets_and_registration(cut, lower_dimension):
    planes = [(0, 0, 1, 1), (0, 0, 1, 0), (1, 0, 0, 0),
              (0, 1, 0, 0), (1, 1, 0, 1), cut]
    constraints = [(0, 1), (1, -1), (2, -1), (3, -1), (4, 1)]
    registries = [reference.PlanePoints(capture_constraints=True) for _ in range(2)]
    for registry in registries:
        registry.planes = [tuple(F(value) for value in plane) for plane in planes]
        registry.labels = [[] for _ in planes]
    old, new = registries
    column = candidate.known_column(new, constraints)
    assert column is not None
    assert old.polytope(constraints) == candidate.register_column(new, column, constraints)
    clipped = candidate.clip(new, column, (5, 1))
    assert (clipped is None) == lower_dimension
    full_constraints = [*constraints, (5, 1)]
    assert old.polytope(full_constraints) == candidate.register_column(new, clipped, full_constraints)
    assert _ownership(old) == _ownership(new)


def test_successive_lower_dimensional_cuts_keep_reference_ownership():
    planes = [(0, 0, 1, 1), (0, 0, 1, 0), (1, 0, 0, 0),
              (0, 1, 0, 0), (1, 1, 0, 1), (1, 0, 0, 0),
              (0, 1, 0, F(1, 2)), (0, 0, 1, 0), (0, 1, 0, 0), (1, 0, 0, -1)]
    constraints = [(0, 1), (1, -1), (2, -1), (3, -1), (4, 1)]
    old, new = [reference.PlanePoints(capture_constraints=True) for _ in range(2)]
    for registry in (old, new):
        registry.planes = [tuple(F(value) for value in plane) for plane in planes]
        registry.labels = [[] for _ in planes]
    column = candidate.known_column(new, constraints)
    for index in range(5, len(planes)):
        column = candidate.clip(new, column, (index, 1))
        assert column is None  # Planar -> smaller planar -> edge -> point -> empty.
        constraints.append((index, 1))
        expected = old.polytope(constraints)
        assert candidate.register_column(new, column, constraints) == expected
        assert _ownership(new) == _ownership(old)
    assert not expected


def test_nonsimple_vertex_registers_all_active_triples_in_legacy_order():
    planes = [(0, 0, 1, 1), (0, 0, 1, 0), (1, 0, 0, 0),
              (0, 1, 0, 0), (1, 1, 0, 1), (1, 1, 1, 0)]
    constraints = [(0, 1), (1, -1), (2, -1), (3, -1), (4, 1)]
    old, new = [reference.PlanePoints(capture_constraints=True) for _ in range(2)]
    for registry in (old, new):
        registry.planes = [tuple(F(value) for value in plane) for plane in planes]
        registry.labels = [[] for _ in planes]
    column = candidate.known_column(new, constraints)
    assert old.polytope(constraints) == candidate.register_column(new, column, constraints)
    clipped = candidate.clip(new, column, (5, -1))  # Redundant, touches only origin.
    assert clipped is not None
    constraints.append((5, -1))
    assert old.polytope(constraints) == candidate.register_column(new, clipped, constraints)
    assert _ownership(old) == _ownership(new)
    origin = new.points.index((F(0), F(0), F(0)))
    triples = [triple for triple, vertex in new.registered_intersections.items() if vertex == origin]
    assert len(triples) == 4


@pytest.mark.parametrize('heights', [(0.,)*4, (-.01, .005, .025, 0.)])
def test_adjacent_transom_chine_cells_and_registration_are_identical(monkeypatch, heights):
    snapshots, emissions, original = [], [], reference.cell_polygons
    def record(owner):
        def wrapped(registry, *args, **kwargs):
            result = owner(registry, *args, **kwargs)
            # Copy mutable ownership before the next cell modifies it.
            from copy import deepcopy
            snapshots.append(deepcopy(_ownership(registry)))
            emissions.append(deepcopy(result))
            return result
        return wrapped
    with monkeypatch.context() as context:
        context.setattr(reference, 'cell_polygons', record(original))
        old_mesh, old_report = reference.build(heights)
    old_snapshots = snapshots[:]
    old_emissions = emissions[:]
    snapshots.clear()
    emissions.clear()
    with monkeypatch.context() as context:
        context.setattr(reference, 'cell_polygons', record(candidate.cell_polygons))
        new_mesh, new_report = reference.build(heights)
    assert snapshots == old_snapshots
    for actual, expected in zip(emissions, old_emissions, strict=True):
        assert actual[0] == expected[0]  # Ordered polygons and physical owners.
        assert actual[1] == expected[1]
        np.testing.assert_array_equal(actual[2], expected[2])
        assert actual[3] == expected[3]
    np.testing.assert_array_equal(new_mesh.vertices, old_mesh.vertices)
    np.testing.assert_array_equal(new_mesh.faces, old_mesh.faces)
    for key in old_report.keys()-{'seconds'}:
        assert new_report[key] == old_report[key], key


def test_canonical_generic_full_arrays_and_template_are_identical(monkeypatch, tmp_path):
    pose = equilibrium()[0]
    pose = InitialState(pose.position+np.array([.123, .087, 0.]), pose.orientation)
    _, old, _ = boundary.construct(tmp_path/'old.json', pose=pose, export_template=True)
    with monkeypatch.context() as context:
        context.setattr(boundary, 'cell_polygons', candidate.cell_polygons)
        _, new, _ = boundary.construct(tmp_path/'new.json', pose=pose, export_template=True)
    assert old.keys() == new.keys()
    for key in old:
        np.testing.assert_array_equal(new[key], old[key], err_msg=key)
