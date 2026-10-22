"""Composite integration must preserve original panels and polynomial moments."""
import numpy as np
import pytest

from tests.experiments.free_surface.boat_target_quadrature import composite_samples
from tests.laplace_reference import triangle_samples


@pytest.mark.parametrize('levels', [1, 2])
def test_composite_preserves_quadratic_moments(levels):
    triangles = np.array([[[0., 0., 0.], [3., 0., 1.], [0., 2., -.4]],
                          [[2., 1., 3.], [2.001, 1., 3.], [2., 4., 2.]]])
    original = triangles.copy()
    points, weights = composite_samples(triangles, 3, levels)
    reference, reference_weights = triangle_samples(triangles, 5)
    assert points.shape == (2, 9*4**levels, 3)
    assert np.all(weights > 0)
    np.testing.assert_allclose(weights.sum(), 1., atol=1e-15)
    # Includes tilted and very thin panels; normalized integrals use the same
    # original physical area, once, outside this rule.
    for i in range(3):
        for j in range(3):
            np.testing.assert_allclose((points[:, :, i]*points[:, :, j])@weights,
                (reference[:, :, i]*reference[:, :, j])@reference_weights,
                rtol=2e-15, atol=2e-15)
    np.testing.assert_array_equal(triangles, original)
