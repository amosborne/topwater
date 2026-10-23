"""Independent integrable fan reference at scalar edges/vertices."""
from math import asinh, log, pi, sqrt
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
import pytest
from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels, panel_influence, QuadratureError
from tests.experiments.free_surface.plane_single_layer import plane_single_layer


def fan_reference(point, triangle):
    result = 0.
    for a, b in zip(triangle[:, :2]-point, np.roll(triangle[:, :2], -1, axis=0)-point):
        determinant = a[0]*b[1]-a[1]*b[0]
        if determinant == 0:
            continue
        integral = quad(lambda t: 1/np.linalg.norm((1-t)*a+t*b), 0, 1,
                        epsabs=1e-11, epsrel=1e-11, limit=128)[0]
        result += determinant*integral
    return abs(result)/(4*pi)


def test_vertex_edge_and_near_edge_against_independent_positive_fans():
    triangle = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    panels = TrianglePanels(triangle[None], 1)
    points = np.array([[0, 0], [.31, 0], [.31, 1e-14], [.31, -1e-14], [.2, .3]])
    actual, metadata = plane_single_layer(points, panels)
    assert_allclose(actual[0, 0], sqrt(2)*asinh(1)/(4*pi), rtol=3e-15)
    # Ordinary edge/interior fans are well resolved. A microscopic signed fan
    # near the edge asks Quadpack to resolve an unnecessary narrow peak.
    expected = [fan_reference(points[i], triangle) for i in (0, 1, 4)]
    assert_allclose(actual[[0, 1, 4], 0], expected, atol=3e-13, rtol=3e-12)
    # Independent kernel-continuity envelope: split the plane into r<2*delta
    # (two disks contribute <=2.5*delta) and the far annulus, where the kernel
    # difference integrates to <=2*delta*log(D/(2*delta)). No near-edge adaptive
    # quadrature result or purported certified Quadpack bound is required.
    for i in (2, 3):
        delta = np.linalg.norm(points[i]-points[1])
        diameter = np.linalg.norm(triangle[:, :2]-points[1], axis=1).max()
        bound = delta*(2.5+2*log(diameter/(2*delta)))
        assert abs(actual[i, 0]-expected[1]) < bound+3e-15
    assert not metadata['certified_error_bound']
    # No small-distance projection: the opposite sides retain their distinct
    # input values, though continuity makes their potentials extremely close.
    assert actual[2, 0] != actual[3, 0]


def test_rotation_reversal_scale_and_existing_scalar_parity():
    triangle = np.array([[.1, -.2, 0.], [1.2, -.1, 0.], [.3, .9, 0.]])
    points = np.array([[1.7, .3], [-.4, -.7], [.8, -.6]])
    panels = TrianglePanels(triangle[None], 1)
    actual, _ = plane_single_layer(points, panels)
    expected = panel_influence(np.column_stack((points, np.zeros(len(points)))), panels).single_layer
    assert_allclose(actual, expected, atol=2e-15, rtol=2e-13)
    angle = .73
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    moved = triangle.copy()
    moved[:, :2] = 2.3*triangle[:, :2] @ rotation.T+np.array([.137, -.081])
    moved_points = 2.3*points @ rotation.T+np.array([.137, -.081])
    scaled, _ = plane_single_layer(moved_points, TrianglePanels(moved[::-1][None], -1))
    assert_allclose(scaled, 2.3*actual, atol=3e-15, rtol=3e-13)


def test_storage_accuracy_and_fallback_budgets_fail_closed():
    panels = TrianglePanels(np.array([[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]]), 1)
    with pytest.raises(MemoryError):
        plane_single_layer([[.2, .2]], panels, LaplaceControls(max_operator_bytes=1))
    tight = LaplaceControls(quadrature_atol=1e-30, quadrature_rtol=1e-30)
    with pytest.raises(QuadratureError, match='fallback budget'):
        plane_single_layer([[.2, .2]], panels, tight, maximum_fallbacks=0)
    with pytest.raises(QuadratureError, match='unresolved'):
        plane_single_layer([[.2, .2]], panels, tight)
    with pytest.raises(ValueError):
        plane_single_layer([[0, 0]], TrianglePanels(panels.triangles+np.array([0, 0, .1]), 1))
