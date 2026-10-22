"""Fixed grid/ownership controls; no mesher or fluid operator."""
import numpy as np
import pytest

from tests.experiments.free_surface.quarter_canoe import POCKET
from tests.experiments.free_surface.quarter_canoe_preflight import (
    WIDTH, BUFFER, envelope_grid, reflected_body)


@pytest.mark.parametrize('envelope,intervals', [(16, (15, 26)), (18, (16, 28))])
def test_actual_grid_retains_phase_envelope_and_all_physical_planes(envelope, intervals):
    x, y, z = envelope_grid(envelope)
    assert (len(x)-1, len(y)-1) == intervals
    np.testing.assert_array_equal(z, [-4., -2.25, -1., -.25, 0.])
    for axis, values in enumerate((x, y)):
        for plane in (0., POCKET[axis], POCKET[axis]+BUFFER, POCKET[axis]+BUFFER+WIDTH):
            assert plane in values
        assert np.max(np.diff(values)) <= (4 if axis == 0 else 2)*np.pi/envelope*(1+8*np.finfo(float).eps)


def test_body_reflections_preserve_orientation_and_reflect_parent_ids():
    triangle = np.array([[[.2, .1, -.3], [.7, .1, -.2], [.2, .4, -.1]]])
    parent_maps = np.array([[3], [8], [9], [12]])
    copies, parents = reflected_body(triangle, np.array([0]), parent_maps)
    np.testing.assert_array_equal(parents, [3, 8, 9, 12])
    normal = np.cross(triangle[0, 1]-triangle[0, 0], triangle[0, 2]-triangle[0, 0])
    for reflected, sign in zip(copies, ((1, 1), (1, -1), (-1, 1), (-1, -1))):
        np.testing.assert_allclose(np.cross(reflected[1]-reflected[0], reflected[2]-reflected[0]),
                                   normal*[*sign, 1], atol=1e-16)
