"""Independent reflected-source integration and chain-rule controls."""
import numpy as np
from numpy.testing import assert_allclose

from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.experiments.free_surface.boat_image_plane import elevated_trace


def test_elevated_pair_matches_independent_source_quadrature():
    panels = TrianglePanels(np.array([[[-.4, -.3, -.5], [.4, -.3, -.5], [0, .4, -.6]]]), -1)
    targets = np.array([[.2, .1, 0], [.8, -.2, -.1]])
    height = .6
    value, vertical = elevated_trace(targets, panels, height)
    sources, weights = triangle_samples(panels.triangles, 24)
    image = sources.copy()
    image[..., 2] = 2*height-image[..., 2]
    real_delta = targets[:, None, None]-sources
    image_delta = targets[:, None, None]-image
    real_distance, image_distance = np.linalg.norm(real_delta, axis=-1), np.linalg.norm(image_delta, axis=-1)
    measure = panels.areas[:, None]*weights/(4*np.pi)
    reference = np.sum(measure*(1/real_distance-1/image_distance), axis=-1)
    derivative = np.sum(measure*(-real_delta[..., 2]/real_distance**3
                                 +image_delta[..., 2]/image_distance**3), axis=-1)
    assert_allclose(value, reference, rtol=2e-12, atol=2e-14)
    assert_allclose(vertical, derivative, rtol=2e-12, atol=2e-14)


def test_reflected_target_vertical_chain_rule():
    panels = TrianglePanels(np.array([[[-.4, -.3, -.5], [.4, -.3, -.5], [0, .4, -.6]]]), -1)
    targets = np.array([[.2, .1, 0], [.8, -.2, -.1]])
    _, derivative = elevated_trace(targets, panels, .6)
    for step in [2e-5, 1e-5]:
        delta = np.array([0, 0, step])
        difference = (elevated_trace(targets+delta, panels, .6)[0]
                      -elevated_trace(targets-delta, panels, .6)[0])/(2*step)
        assert_allclose(derivative, difference, rtol=2e-9, atol=2e-11)
