"""Optional-GPU kinetic work and action-consistent geometric Taylor controls."""
import numpy as np
import pytest

pytest.importorskip('cupy')
from .experiments.free_surface import cuda_operator  # Existing Windows CUDA setup.
from .experiments.free_surface.p4_metric_action import MetricAction, TangentMetricError
from .experiments.free_surface.p4_geometry import integrals, InvalidGeometryError
from .test_p4_geometry import curved_cell


def test_exact_cached_kinetic_gradient_matches_independent_physical_work():
    mesh, points = curved_cell(.7);nodes = mesh.elements[0]
    phi = np.sin(points[:, 0]+.3*points[:, 1]-.2*points[:, 2])
    cache = MetricAction(points[mesh.elements], mesh.elements, len(points))
    current = cache.kinetic(phi);oracle = integrals(points[nodes], phi[nodes])
    assert current['kinetic_energy'] == pytest.approx(oracle['kinetic_energy'], rel=3e-13)
    np.testing.assert_allclose(current['kinetic_gradient'][nodes], oracle['kinetic_gradient'], rtol=3e-12, atol=3e-10)
    direction = np.column_stack((points[:, 0]*points[:, 1], points[:, 1]*points[:, 2], points[:, 0]*points[:, 2]))
    expected = np.sum(current['kinetic_gradient']*direction)
    errors = []
    for h in (1e-3, 5e-4):
        fd = (integrals((points+h*direction)[nodes], phi[nodes])['kinetic_energy']
              -integrals((points-h*direction)[nodes], phi[nodes])['kinetic_energy'])/(2*h)
        errors.append(abs(fd-expected))
    assert errors[1] < .4*errors[0]+1e-9
    action, _ = cache.action(phi)
    assert 500*phi@action == pytest.approx(current['kinetic_energy'], rel=3e-13)


def test_tangent_energy_has_matching_base_force_and_distinct_guards():
    mesh, points = curved_cell(.7);nodes = mesh.elements[0]
    phi = np.sin(points[:, 0]+.3*points[:, 1]-.2*points[:, 2])
    cache = MetricAction(points[mesh.elements], mesh.elements, len(points))
    direction = np.column_stack((points[:, 0]*points[:, 1], points[:, 1]*points[:, 2], points[:, 0]*points[:, 2]))
    errors = []
    for delta in (.01, .005):
        moved = points+delta*direction
        tangent = cache.tangent(moved[mesh.elements]);candidate = tangent.kinetic(phi)
        oracle = integrals(moved[nodes], phi[nodes])
        exact_force = np.zeros_like(points);exact_force[nodes] = oracle['kinetic_gradient']
        errors.append((abs(candidate['kinetic_energy']-oracle['kinetic_energy']),
                       np.linalg.norm(candidate['kinetic_gradient']-exact_force)))
        action, _ = tangent.action(phi)
        assert 500*phi@action == pytest.approx(candidate['kinetic_energy'], rel=3e-13)
    assert errors[1][0] < .3*errors[0][0]  # O(delta²) energy, not a field-only gate.
    assert errors[1][1] < .55*errors[0][1]  # O(delta) shape-force limitation retained.
    h = 1e-4
    plus = cache.tangent((moved+h*direction)[mesh.elements]).kinetic(phi)['kinetic_energy']
    minus = cache.tangent((moved-h*direction)[mesh.elements]).kinetic(phi)['kinetic_energy']
    expected = np.sum(candidate['kinetic_gradient']*direction)
    assert (plus-minus)/(2*h) == pytest.approx(expected, abs=3e-8)
    # A large physical one-axis stretch remains valid, while its first-order
    # tensor loses positivity. Refresh, rather than a geometry verdict.
    stretched = points.copy();stretched[:, 0] *= 3
    with pytest.raises(TangentMetricError):
        cache.tangent(stretched[mesh.elements])
    inverted = points.copy();inverted[:, 0] *= -1
    with pytest.raises(InvalidGeometryError):
        cache.tangent(inverted[mesh.elements])
