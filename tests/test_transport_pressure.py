"""Rigid-transport derivative against an independent analytic harmonic mode."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

pytest.importorskip('cupy')

from topwater.frames import rotation_matrix
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.reference_boat import simulation_inputs
from tests.reference_convex_boat import box_inputs
from tests.experiments.free_surface.boat_graph_bvp import GraphBody
from tests.experiments.free_surface.boat_transport_pressure import normal_transport_lift
from tests.experiments.free_surface.cuda_graph import cp


@pytest.mark.parametrize('make_inputs', [simulation_inputs, box_inputs], ids=['center-console', 'box'])
def test_normal_rigid_transport_lift_matches_spatial_difference_on_actual_hull(make_inputs):
    inputs = make_inputs()
    state = inputs.initial_state
    panels = TrianglePanels(hydrostatics(inputs.body, state).wetted_triangles, -1)
    origin = state.position[:2]-32.
    context = GraphBody(panels, np.zeros((17, 17)), 64., origin, 0., 0., 4)
    size = context.graph.padded_shape[0]
    i, j = cp.meshgrid(cp.arange(size), cp.arange(size))
    spectrum = cp.fft.fft2(cp.cos(2*np.pi*(2*i+j)/size)).ravel()
    wavevector = 2*np.pi*np.array([2., 1.])/64.
    wave = np.linalg.norm(wavevector)
    omega = rotation_matrix(state.orientation)@state.angular_velocity

    def analytic(points):
        phase = (points[:, :2]-origin)@wavevector
        gradient = np.exp(wave*points[:, 2, None])*np.column_stack((
            -wavevector[0]*np.sin(phase), -wavevector[1]*np.sin(phase), wave*np.cos(phase)))
        speed = state.linear_velocity+np.cross(omega, points-state.position)
        return gradient, np.sum(speed*gradient, axis=1)

    velocity, derivative = normal_transport_lift(context, spectrum, state)
    assert_allclose(velocity, analytic(panels.centroids)[0], atol=2e-14, rtol=3e-13)
    for step in [1e-3, 5e-4]:
        difference = (analytic(panels.centroids+step*panels.normals)[1]
                      -analytic(panels.centroids-step*panels.normals)[1])/(2*step)
        assert_allclose(derivative, difference, atol=3e-10, rtol=3e-8)
