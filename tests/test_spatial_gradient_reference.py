"""Independent integration/differencing checks on the spatial audit reference."""
import numpy as np
from numpy.testing import assert_allclose

from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels, panel_influence
from tests.reference_boat import simulation_inputs
from tests.laplace_reference import triangle_samples
from tests.experiments.free_surface.panel_gradient_reference import checked_influence


def panel():
    inputs = simulation_inputs()
    wet = hydrostatics(inputs.body, inputs.initial_state, density=1025.)
    # Actual long panel abutting the reference boat's waterline.
    triangles = wet.wetted_triangles
    wet_edge = np.sum(triangles[:, :, 2] == 0, axis=1) == 2
    eligible = triangles[wet_edge]
    diameter = np.linalg.norm(eligible[:, 1]-eligible[:, 0], axis=1)
    return TrianglePanels(eligible[np.argmax(diameter)][None], -1)


def test_bulk_gradients_against_independent_area_quadrature():
    panels = panel()
    target = panels.centroids[0]+np.array([12., 5., -2.])
    actual = checked_influence(target[None], panels)
    samples, weights = triangle_samples(panels.triangles, 32)
    delta = target-samples[0]
    radius = np.linalg.norm(delta, axis=1)
    normal = panels.normals[0]
    height = delta@normal
    measure = panels.areas[0]*weights/(4*np.pi)
    gs = np.sum(-measure[:, None]*delta/radius[:, None]**3, axis=0)
    gk = np.sum(measure[:, None]*(normal/radius[:, None]**3
        -3*height[:, None]*delta/radius[:, None]**5), axis=0)
    assert_allclose(actual.gradient_single_layer[0, 0], gs, rtol=1e-9, atol=1e-12)
    assert_allclose(actual.gradient_double_layer[0, 0], gk, rtol=1e-8, atol=1e-12)


def test_near_contact_gradients_against_independent_five_point_derivative():
    panels = panel()
    triangle = panels.triangles[0]
    edge = triangle[triangle[:, 2] == 0]
    point = .63*edge[0]+.37*edge[1]
    tangent = edge[1, :2]-edge[0, :2]
    outward = np.array([tangent[1], -tangent[0]])/np.linalg.norm(tangent)
    if outward@(point[:2]-panels.centroids[0, :2]) < 0:
        outward = -outward
    point[:2] += .001*outward
    actual = checked_influence(point[None], panels)
    expected = []
    for step in [2e-6, 1e-6]:
        values = []
        for axis in np.eye(3):
            targets = point+np.array([-2., -1., 1., 2.])[:, None]*step*axis
            result = panel_influence(targets, panels)
            entries = np.column_stack((result.single_layer[:, 0], result.double_layer[:, 0]))
            values.append(np.array([1., -8., 8., -1.])@entries/(12*step))
        expected.append(np.array(values).T)
    assert_allclose(expected[0], expected[1], rtol=2e-6, atol=2e-7)
    assert_allclose(np.stack((actual.gradient_single_layer[0, 0],
                             actual.gradient_double_layer[0, 0])), expected[1], rtol=2e-6, atol=2e-7)
