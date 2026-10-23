"""Map-free harmonic actions against an analytic oblique Fourier mode."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

pytest.importorskip('cupy')
from tests.experiments.free_surface.boat_frozen_field import harmonic_action
from tests.experiments.free_surface.cuda_graph import cp


def test_harmonic_action_matches_oblique_mode_and_its_normal_derivative():
    size, length = 35, 64.
    x, y = cp.meshgrid(cp.arange(size), cp.arange(size))
    spectrum = cp.fft.fft2(cp.cos(2*np.pi*(2*x+y)/size)).ravel()
    origin = np.array([-32., -32.])
    points = np.array([[1., -.2, -.3], [-2., .5, -.1], [.3, 1.2, -2.]])
    normals = np.array([[0., 0., 1.], [.6, 0., .8], [0., -.8, .6]])
    k = 2*np.pi*np.array([2., 1.])/length
    wave = np.linalg.norm(k)
    phase = (points[:, :2]-origin)@k
    amplitude = np.exp(wave*points[:, 2])
    potential = amplitude*np.cos(phase)
    normal = amplitude*(normals[:, 2]*wave*np.cos(phase)
                         -(normals[:, :2]@k)*np.sin(phase))
    assert_allclose(harmonic_action(points, normals, spectrum, size, length, origin),
                    np.column_stack((potential, normal)), atol=2e-15, rtol=2e-14)
