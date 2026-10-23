"""Independent normalization and derivative checks for the work audit."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad

from tests.experiments.free_surface.boat_surface_work import disk_quadrature, fourier_lift


def test_disk_weights_include_full_angular_area():
    radius = 1.3
    _, weights = disk_quadrature(np.array([.2, -.4]), radius, 96)
    integral = 2*np.pi*quad(lambda r: r*np.exp(1-1/(1-(r/radius)**2)),
                            0, radius, epsabs=1e-13, epsrel=1e-13)[0]
    assert_allclose(weights.sum(), integral, rtol=1e-12, atol=1e-13)


def test_fourier_normal_matches_independent_potential_difference():
    points = np.array([[.3, .1, -.4], [1.1, -.7, -.9]])
    normals = np.array([[.3, .4, .5], [-.2, .5, .3]])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    center, origin = np.array([.2, -.3]), np.array([-8., -8.])
    _, actual, _ = fourier_lift(points, normals, center, 1.3, 33, 16., origin)
    for step in [2e-5, 1e-5]:
        plus = fourier_lift(points+step*normals, normals, center, 1.3, 33, 16., origin)[0]
        minus = fourier_lift(points-step*normals, normals, center, 1.3, 33, 16., origin)[0]
        assert_allclose(actual, (plus-minus)/(2*step), rtol=2e-9, atol=1e-11)
