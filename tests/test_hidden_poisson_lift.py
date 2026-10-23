"""Independent adaptive-integral and derivative controls for the hidden lift."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.integrate import quad

from tests.experiments.free_surface.boat_hidden_extension import poisson_lift


def test_axis_value_and_gradient_against_adaptive_radial_integral():
    center, radius = np.array([1.2, -.4]), .65
    depths = [.1, .7, 3.]
    points = np.array([[*center, -h] for h in depths])
    value, gradient = poisson_lift(points, center, radius, 64)
    expected_value, expected_vertical = [], []
    for h in depths:
        def measure(r):
            return r*np.exp(1-1/(1-(r/radius)**2)) if r < radius else 0.
        expected_value.append(quad(lambda r: h*measure(r)/(r*r+h*h)**1.5,
                                   0, radius, epsabs=1e-12, epsrel=1e-12)[0])
        expected_vertical.append(quad(lambda r: (2*h*h-r*r)*measure(r)/(r*r+h*h)**2.5,
                                      0, radius, epsabs=1e-12, epsrel=1e-12)[0])
    assert_allclose(value, expected_value, rtol=2e-11, atol=1e-12)
    assert_allclose(gradient[:, 2], expected_vertical, rtol=2e-10, atol=1e-12)
    assert_allclose(gradient[:, :2], 0, atol=1e-14)


def test_off_axis_gradient_against_potential_differences_and_plane_contract():
    center, radius = np.array([.2, -.3]), .65
    point, step = np.array([.7, .1, -.2]), 1e-5
    _, gradient = poisson_lift(point[None], center, radius, 64)
    difference = []
    for direction in np.eye(3):
        values, _ = poisson_lift(point+np.array([-2., -1., 1., 2.])[:, None]*step*direction,
                                 center, radius, 64)
        difference.append(np.array([1., -8., 8., -1.])@values/(12*step))
    assert_allclose(gradient[0], difference, rtol=1e-8, atol=1e-10)
    with pytest.raises(ValueError, match='boundary limit'):
        poisson_lift(np.array([[*center, 0.]]), center, radius, 24)
    value, gradient = poisson_lift(np.array([[2., 1., 0.]]), center, radius, 64)
    assert value[0] == 0
    assert_allclose(gradient[0, :2], 0, atol=0)
    assert gradient[0, 2] < 0
