"""Static clipping-frame covariance and the reported one-ulp grid regression."""
import numpy as np
from numpy.testing import assert_allclose

from tests.reference_convex_boat import box_inputs
from tests.experiments.free_surface.boat_nonplanar_geometry import local_geometry_inputs, local_axis
from tests.experiments.free_surface.incident_baseline_benchmark import coefficients, vectorized


def test_local_geometry_frame_preserves_the_physical_wave_realization():
    inputs = box_inputs()
    origin, pose, waves = local_geometry_inputs(inputs)
    assert_allclose(pose.position+origin, inputs.initial_state.position, atol=0, rtol=0)
    for field in ['orientation', 'linear_velocity', 'angular_velocity']:
        assert_allclose(getattr(pose, field), getattr(inputs.initial_state, field), atol=0, rtol=0)
    points = np.array([[1.2, -.4, -.2], [-2.1, .7, 0.], [0., 0., .1]])
    expected = vectorized(points+origin, coefficients(inputs.waves), .37)
    actual = vectorized(points, coefficients(waves), .37)
    assert_allclose(actual, expected, atol=2e-13, rtol=2e-13)
    assert waves.seed == inputs.waves.seed


def test_local_parent_grid_keeps_intended_zero_and_binary_half_steps():
    axis = local_axis(2.4, 4)
    assert axis[2] == 0.
    assert axis[1] == -1.2 and axis[3] == 1.2
    assert_allclose(axis, -axis[::-1], atol=0, rtol=0)
