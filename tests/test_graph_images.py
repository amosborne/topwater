"""Interior-point checks of the local Fourier-graph image bound."""
import numpy as np
from numpy.testing import assert_allclose

from tests.laplace_reference import triangle_samples
from tests.experiments.free_surface.boat_graph_images import image_height


def test_flat_contact_recovers_the_flat_image_plane_to_roundoff():
    triangles = np.array([[[-1., -.5, 1.2], [1., -.5, 1.2], [0., .5, .9]]])
    height, bound = image_height(triangles, np.full((33, 33), 1.2), [16., 16.], [-8., -8.])
    assert_allclose(height, 1.2, rtol=0, atol=3e-11)
    assert bound['minimum_image_graph_gap_bound_m'] >= 0


def test_faraway_wave_crest_does_not_set_image_height():
    size, length = 129, 64.
    x = -length/2+np.arange(size)*length/size
    k = 2*np.pi/length
    eta = np.tile(.8*np.cos(k*(x-20)), (size, 1))
    triangles = np.array([[[-1., -.5, -.8], [1., -.5, -.8], [0., .5, -.8]]])
    height, bound = image_height(triangles, eta, [length, length], [-length/2, -length/2])
    points, _ = triangle_samples(triangles, 24)
    image_z = 2*height-points[..., 2]
    graph = .8*np.cos(k*(points[..., 0]-20))
    assert np.all(image_z > graph)
    assert bound['minimum_image_graph_gap_bound_m'] >= 0
    assert bound['global_graph_node_maximum_m'] > .79
    assert height < 0
