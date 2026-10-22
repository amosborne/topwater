"""Moving polynomial graph charts, not fabricated rigid-body canoe cases."""
import numpy as np
import pytest
from scipy.sparse import csr_matrix

from .experiments.free_surface.p4_tetra import indices
from .experiments.free_surface.p4_surface import integrals


def chart(eta, transport):
    """Every height moves z; height-dependent XY bubbles fix the outer edge.

    Two declared interior XY transport parameters are chart motions, NOT
    physical six-DOF body motions. Four unused twist slots remain zero.
    """
    barycentric = indices(3)/4
    r, s = barycentric[:, 1], barycentric[:, 2]
    bubble = 27*np.prod(barycentric, axis=1)
    second = bubble*(r-s)
    a = .025*np.sin(np.arange(15)+.4)
    b = .018*np.cos(np.arange(15)*.7)
    coordinates = np.column_stack((r+bubble*(a@eta+transport[0]),
                                   s+second*(b@eta+transport[1]), eta))
    height = (np.outer(bubble, a), np.outer(second, b), np.eye(15))
    body = np.zeros((15, 3, 6))
    body[:, 0, 0], body[:, 1, 1] = bubble, second
    return coordinates, height, body


def test_planar_graph_exact_moments_and_transport_signs():
    barycentric = indices(3)/4
    r, s = barycentric[:, 1], barycentric[:, 2]
    c, hx, hy = .3, .17, -.11
    transport = np.array([.02, -.01])
    # Solve the actual chart for a plane z=c+hx*x+hy*y; no stale geometry.
    coordinates, height, body = chart(np.zeros(15), transport)
    eta = np.linalg.solve(np.eye(15)-hx*height[0]-hy*height[1],
                          c+hx*coordinates[:, 0]+hy*coordinates[:, 1])
    coordinates, height, body = chart(eta, transport)
    result = integrals(coordinates, tuple(csr_matrix(value) for value in height), body)
    expected_eta = result['mass']@(height[2]-hx*height[0]-hy*height[1])
    expected_body = result['mass']@(body[:, 2]-hx*body[:, 0]-hy*body[:, 1])
    np.testing.assert_allclose(result['eta_map'], expected_eta, atol=2e-14, rtol=1e-12)
    np.testing.assert_allclose(result['body_transport'], expected_body, atol=2e-14, rtol=1e-12)
    # The physical projected triangle stays fixed despite moving interior XY.
    assert result['volume_boundary'] == pytest.approx(c/2+(hx+hy)/6, abs=2e-14)
    first_moment = .5*(c*c/2+c*(hx+hy)/3+(hx*hx+hy*hy+hx*hy)/12)
    assert result['gravity_boundary'] == pytest.approx(9810*first_moment, abs=2e-10)
    assert np.linalg.eigvalsh(result['mass']).min() > 0
    assert np.linalg.matrix_rank(result['eta_map']) == 15
    assert np.max(abs(result['eta_map']-result['mass'])) > 1e-5
    assert np.linalg.norm(result['body_transport'][:, :2]) > 1e-3


def test_every_high_height_mode_has_actual_geometry_and_physical_volume_work():
    barycentric = indices(3)/4
    r, s = barycentric[:, 1], barycentric[:, 2]
    eta = .2+.12*r-.09*s+.07*r*s+.04*r**3
    transport = np.array([.025, -.02])
    coordinates, height, body = chart(eta, transport)
    result = integrals(coordinates, height, body)
    # A pure high mode vanishes on ALL triangle edges and all original corners.
    high_mode = 27*np.prod(barycentric, axis=1)
    direction = high_mode+.2*np.sin(np.arange(15)*.8)
    twist = np.array([.3, -.2, 0., 0., 0., 0.])
    flux = result['eta_map']@direction+result['body_transport']@twist
    volume_rate = np.ones(15)@flux
    expected_velocity = np.column_stack([matrix@direction for matrix in height])+body@twist
    errors = []
    for step in (2e-3, 1e-3):
        plus = chart(eta+step*direction, transport+step*twist[:2])
        minus = chart(eta-step*direction, transport-step*twist[:2])
        actual_velocity = (plus[0]-minus[0])/(2*step)
        np.testing.assert_allclose(actual_velocity, expected_velocity, atol=8e-14, rtol=2e-11)
        actual_rate = (integrals(*plus)['volume_boundary']-integrals(*minus)['volume_boundary'])/(2*step)
        errors.append(abs(actual_rate-volume_rate))
    assert errors[-1] < 2e-8
    assert errors[-1] < .4*errors[0]+1e-12
    # This mode would be silently lost by a corner-only geometry lift.
    plus = chart(eta+1e-3*high_mode, transport)[0]
    assert np.max(abs(plus[:, 2]-coordinates[:, 2])) > 5e-4
    assert np.linalg.norm(result['eta_map']@high_mode) > .02
    assert np.ones(15)@(result['eta_map']@high_mode) > .1


def test_curved_surface_polynomial_quadrature_and_sampled_graph_guard():
    barycentric = indices(3)/4
    eta = .2+.03*np.sin(np.arange(15))
    data = chart(eta, [.02, -.015])
    coarse, fine = integrals(*data, order=8), integrals(*data, order=12)
    for key in ('eta_map', 'body_transport', 'mass', 'volume_boundary', 'gravity_boundary'):
        np.testing.assert_allclose(coarse[key], fine[key], atol=2e-11, rtol=2e-12)
    reversed_coordinates = data[0].copy()
    reversed_coordinates[:, 0] *= -1
    with pytest.raises(ValueError, match='graph Jacobian'):
        integrals(reversed_coordinates, data[1], data[2])
