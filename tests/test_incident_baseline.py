"""Independent Airy identities for the actual canonical-sea incident evaluator."""
import numpy as np
from numpy.testing import assert_allclose

from topwater.sea_states import realize
from tests.experiments.free_surface.incident_baseline_benchmark import (
    B, G, coefficients, fused_cpu, vectorized,
)


def test_actual_sea_airy_identities_and_potential_gradients():
    values = coefficients(realize(beam_m=B, seed=42, gravity=G))
    points = np.array([[0., 0., 0.], [2.7, -.8, -.3], [-4., .6, -1.9]])
    time, step = .37, 1e-5
    field = fused_cpu(points, values, time)
    assert_allclose(field, vectorized(points, values, time), rtol=2e-13, atol=3e-14)
    for axis in range(3):
        direction = np.zeros(3)
        direction[axis] = step
        derivative = (fused_cpu(points+direction, values, time)
                      -fused_cpu(points-direction, values, time))/(2*step)
        # Spatial velocity must be the gradient of the independently differenced
        # potential, including vertical decay rather than a surface-only formula.
        assert_allclose(derivative[:, 4], field[:, 6+axis], rtol=2e-8, atol=2e-10)
        if axis < 2:
            assert_allclose(derivative[:, 0], field[:, 1+axis], rtol=2e-8, atol=2e-10)
    derivative = (fused_cpu(points, values, time+step)-fused_cpu(points, values, time-step))/(2*step)
    assert_allclose(derivative[:, 4], field[:, 5], rtol=2e-8, atol=2e-10)
    assert_allclose(derivative[:, 0], field[:, 3], rtol=2e-8, atol=2e-10)
    surface = points.copy()
    surface[:, 2] = 0
    trace = fused_cpu(surface, values, time)
    # Linear free-surface kinematic and dynamic conditions, each pointwise.
    assert_allclose(trace[:, 3], trace[:, 8], rtol=2e-13, atol=3e-14)
    assert_allclose(trace[:, 5]+G*trace[:, 0], 0, atol=3e-14)
