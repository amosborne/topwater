"""Actual canonical facet subdivisions and equal-area overlap counterexample."""
import json

import numpy as np
import pytest

from topwater import SimulationInputs
from topwater._clipping import clip_flat_water
from topwater.frames import rotation_matrix
from topwater.hydrostatics import equilibrate_heave
from .reference_hull import input_data
from .experiments.capytaine_mesh_cost import subdivide
from .experiments.free_surface.canoe_wet_reference import verify_coverage


def test_actual_canoe_subdivisions_cover_identical_clipped_physical_facets():
    data = SimulationInputs.from_json(json.dumps(input_data(2)))
    rest, _ = equilibrate_heave(data.body, data.initial_state)
    world = (data.body.mesh.vertices-data.body.cg)@rotation_matrix(rest.orientation).T+rest.position
    wet, parent, *_ = clip_flat_water(world, data.body.mesh.faces, data.body.mesh.physical_faces, 0., True)
    candidate, parents = wet.copy(), parent.copy()
    for _ in range(2):
        candidate, parents = subdivide(candidate), np.tile(parents, 4)
        report = verify_coverage(candidate, parents, wet, parent)
        assert report['candidate_triangles'] == len(candidate)
        assert report['physical_parents'] == len(np.unique(parent))
        assert report['field_accuracy_proven'] is False
    # Missing panels and reversed orientation must not pass just by lying on
    # the original parent plane.
    with pytest.raises(ValueError, match='gap or overlap|loses'):
        verify_coverage(candidate[1:], parents[1:], wet, parent)
    bad = candidate.copy();bad[0] = bad[0, ::-1]
    with pytest.raises(ValueError, match='reverses'):
        verify_coverage(bad, parents, wet, parent)


def test_equal_area_hole_and_duplicate_fail_explicit_overlap_check():
    target = np.array([[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]])
    small = np.array([[0., 0., 0.], [.5, 0., 0.], [0., 1., 0.]])
    # Twice half the area equals the target area, but half the facet is absent.
    with pytest.raises(ValueError, match='gap or overlap'):
        verify_coverage(np.array([small, small]), np.array([7, 7]), target, np.array([7]))
