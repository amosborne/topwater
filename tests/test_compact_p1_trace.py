"""Independent Poisson-integral and differential checks, including shallow targets."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.integrate import quad
from tests.experiments.free_surface.compact_p1_trace import CompactP1Trace


def patch():
    x, y = np.linspace(-1., 1., 5), np.linspace(-.8, .8, 5)
    xx, yy = np.meshgrid(x, y)
    r = (xx/.8)**2+(yy/.65)**2
    values = np.zeros_like(r)
    active = r < 1
    values[active] = np.exp(1-1/(1-r[active]))
    return CompactP1Trace(x, y, values)


def independent_poisson(patch, target):
    """Adaptive edge parameter after exact radial integration of Poisson's kernel.

    Signed fans from target XY cover each triangle. This directly integrates
    f*h/(2pi*(r²+h²)^(3/2)); it calls no Rankine/solid-angle routine. The
    adaptive rule also applies when depth is much smaller than a source cell.
    """
    projection, height = target[:2], -target[2]
    total = 0.
    for ids in patch.faces:
        xy = np.array([[patch.x[i % len(patch.x)], patch.y[i//len(patch.x)]] for i in ids])
        values = patch.values.ravel()[ids]
        plane = np.linalg.solve(np.column_stack((np.ones(3), xy)), values)
        f0 = plane[0]+plane[1:]@projection
        for a, b in zip(xy-projection, np.roll(xy-projection, -1, axis=0)):
            determinant = a[0]*b[1]-a[1]*b[0]
            if determinant == 0:
                continue
            def integrand(t):
                radial = a+t*(b-a)
                length = np.linalg.norm(radial)
                root = np.hypot(length, height)
                first = 1/(root*(root+height))
                second = height*(np.arcsinh(length/height)-length/root)/length**3
                return determinant/(2*np.pi)*(f0*first+(plane[1:]@radial)*second)
            value, error = quad(integrand, 0., 1., epsabs=2e-13, epsrel=2e-12, limit=100)
            assert error < 2e-11
            total += value
    return total


def test_bulk_matches_independent_integral_and_gradient_difference():
    p = patch()
    points = np.array([[.17, -.09, -.02], [-.27, .13, -.35], [1.4, .2, -.07]])
    potential, velocity = p.bulk(points)
    expected = np.array([independent_poisson(p, point) for point in points])
    assert_allclose(potential, expected, rtol=2e-11, atol=3e-12)
    for axis in range(3):
        errors = []
        for step in [2e-4, 1e-4]:
            direction = np.eye(3)[axis]*step
            derivative = (p.bulk(points+direction)[0]-p.bulk(points-direction)[0])/(2*step)
            errors.append(np.max(abs(derivative-velocity[:, axis])))
        assert errors[1] < .35*errors[0]+3e-10


def test_owned_surface_limits_and_global_return_tail():
    p = patch()
    indices = np.array([1, 8, 19])
    surface = p.surface_velocity(indices)
    errors = []
    for depth in [2e-4, 1e-4]:
        points = p.panels.centroids[indices]-[0., 0., depth]
        errors.append(np.max(abs(p.bulk(points)[1]-surface)))
    assert errors[1] < .6*errors[0]
    xy = np.array([[1.7, .31], [20., 0.], [40., 0.]])
    returned = p.outside_velocity(xy)
    assert_allclose(returned[:, :2], 0., atol=0)
    assert np.all(returned[:, 2] < 0)
    mass = np.sum(p.panels.areas*p.values.ravel()[p.faces].mean(axis=1))
    assert_allclose(returned[1:, 2], -mass/(2*np.pi*xy[1:, 0]**3), rtol=.01)
    assert_allclose(p.bulk(np.column_stack((xy, np.full(len(xy), -1e-5))))[1], returned, atol=4e-6)


def test_owner_and_boundary_query_contracts():
    p = patch()
    with pytest.raises(ValueError, match='perimeter'):
        CompactP1Trace(p.x, p.y, np.ones_like(p.values))
    with pytest.raises(ValueError, match='below'):
        p.bulk(p.panels.centroids[:1])
    with pytest.raises(ValueError, match='outside'):
        p.outside_velocity([[p.x[0], p.y[0]], [0., 0.]])
    with pytest.raises(ValueError, match='indices'):
        p.surface_velocity([-1])
    x, y, values = p.x.copy(), p.y.copy(), p.values.copy()
    owned = CompactP1Trace(x, y, values)
    x[:] = 100.; y[:] = 100.; values[:] = 100.
    assert_allclose(owned.bulk([[.1, .2, -.1]])[0], p.bulk([[.1, .2, -.1]])[0], atol=0, rtol=0)
