"""Motion-frame kinematics and valid calm-rest zero-potential startup."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

pytest.importorskip('cupy')

from topwater import InitialState
from topwater.frames import rotation_matrix
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.reference_boat import simulation_inputs
from tests.experiments.free_surface.boat_graph_bvp import GraphBody
from tests.experiments.free_surface.boat_graph_time_derivative import prescribed_pose
from tests.experiments.free_surface.cuda_graph import cp


def test_prescribed_world_velocity_matches_moving_material_points():
    quaternion = np.array([.96, .12, -.16, .2])
    initial = InitialState([.3, -.2, .1], quaternion/np.linalg.norm(quaternion),
                           [.03, -.02, .01], [.01, -.015, .02])
    material = np.array([[1.2, -.7, .2], [-.5, .4, -.3]])
    time = .13
    pose = prescribed_pose(initial, time)
    points = pose.position+material@rotation_matrix(pose.orientation).T
    omega = rotation_matrix(pose.orientation)@pose.angular_velocity
    expected = pose.linear_velocity+np.cross(omega, points-pose.position)
    for step in [1e-3, 5e-4]:
        offsets = [prescribed_pose(initial, time+s) for s in [step, -step]]
        world = [p.position+material@rotation_matrix(p.orientation).T for p in offsets]
        assert_allclose((world[0]-world[1])/(2*step), expected, atol=2e-10, rtol=0)


def test_calm_rest_body_bvp_handles_exact_zero_rhs_without_nan():
    inputs = simulation_inputs()
    wet = hydrostatics(inputs.body, inputs.initial_state)
    panels = TrianglePanels(wet.wetted_triangles, -1)
    context = GraphBody(panels, np.zeros((17, 17)), 64.,
                        inputs.initial_state.position[:2]-32., 0., 0., 4)
    solution = context.solve(np.zeros(17**2), np.zeros(len(panels)))
    assert_allclose(solution.potential, 0., atol=0, rtol=0)
    assert_allclose(solution.density, 0., atol=0, rtol=0)
    assert_allclose(cp.asnumpy(solution.spectrum), 0., atol=0, rtol=0)
    assert solution.diagnostics['true_relative_residual'] == 0.
    assert solution.diagnostics['gmres_iterations'] == 0
