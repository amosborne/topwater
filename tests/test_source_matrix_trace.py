"""Matrix parity on two real convex hulls and changed clipping poses."""
import numpy as np
import pytest
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.reference_boat import simulation_inputs
from tests.reference_convex_boat import box_inputs
from tests.experiments.free_surface.boat_graph_geometry import constant_twist_pose
from tests.experiments.free_surface.boat_source_matrix_audit import old_source_matrix, matrix_accuracy_request
from tests.experiments.free_surface.boat_spatial_work import source_matrix, subdivide


@pytest.mark.parametrize('factory', [simulation_inputs, box_inputs])
@pytest.mark.parametrize('time', [0., .15])
def test_complete_source_matrix_matches_old_actual_hull(factory, time):
    inputs = factory()
    pose = constant_twist_pose(inputs.initial_state, time)
    wet = hydrostatics(inputs.body, pose)
    panels = TrianglePanels(subdivide(wet.wetted_triangles), -1)
    old, new = old_source_matrix(panels), source_matrix(panels)
    # Old unused-gK rejection can trigger Decimal recovery of gS as well.
    # Use existing componentwise gS requests, contracted with target normals;
    # do not pick a tolerance by padding the observed matrix difference.
    budget = matrix_accuracy_request(panels)
    rounding = 64*np.finfo(float).eps*np.maximum(1., abs(old))
    assert np.all(abs(new-old) <= budget+rounding)
    np.testing.assert_array_equal(np.diag(new), np.diag(old))
