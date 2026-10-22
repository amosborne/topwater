"""Known 3D waves test signed PML corners, conormals and surface measure.

An independent analytic conormal closes integration by parts for every total
degree <=4 Cartesian test polynomial. No global matrix or PDE solve is used.
"""
from itertools import product

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from .experiments.free_surface.pml_corner import coordinates, coefficients


def _rule(bounds, order, fixed=None):
    t, w = leggauss(order)
    axes = [i for i in range(3) if fixed is None or i != fixed[0]]
    grids = np.meshgrid(*[.5*((bounds[i, 1]-bounds[i, 0])*t+bounds[i].sum()) for i in axes], indexing='ij')
    weights = np.meshgrid(*[.5*(bounds[i, 1]-bounds[i, 0])*w for i in axes], indexing='ij')
    points = np.empty((order**len(axes), 3))
    for i, grid in zip(axes, grids):
        points[:, i] = grid.ravel()
    if fixed is not None:
        points[:, fixed[0]] = fixed[1]
    return points, np.prod(weights, axis=0).ravel()


def _polynomials(points):
    powers = np.array([p for p in product(range(5), repeat=3) if sum(p) <= 4])
    value = np.prod(points[:, None, :]**powers[None, :, :], axis=-1)
    gradient = np.zeros(value.shape+(3,))
    for axis in range(3):
        reduced = powers.copy();reduced[:, axis] = np.maximum(0, reduced[:, axis]-1)
        gradient[..., axis] = powers[:, axis]*np.prod(points[:, None, :]**reduced[None, :, :], axis=-1)
    return value, gradient


def _wave(points, omega, angle, onset, width):
    k = brentq(lambda k: 9.81*k*np.tanh(4*k)-omega**2, 1e-6, 10.)
    kx, ky = k*np.cos(angle), k*np.sin(angle)
    mapped, stretch = coordinates(points, onset, width, 8., omega)
    phase = np.exp(1j*(kx*mapped[:, 0]+ky*mapped[:, 1]))
    vertical = np.cosh(k*(points[:, 2]+4))/np.cosh(4*k)
    vertical_z = k*np.sinh(k*(points[:, 2]+4))/np.cosh(4*k)
    phi = phase*vertical
    sx, sy = stretch.T
    gradient = np.column_stack((1j*kx*sx*phi, 1j*ky*sy*phi, phase*vertical_z))
    # Written directly from transformed physical normal flux, independent of
    # multiplying the volume tensor by computational gradients.
    conormal = np.column_stack((1j*kx*sy*phi, 1j*ky*sx*phi, sx*sy*phase*vertical_z))
    return phi, gradient, conormal, stretch


def _balance(omega, angle, bounds, order, *, wrong_corner=False):
    onset, width = np.array([.6, .4]), np.array([.8, .7])
    points, weights = _rule(bounds, order)
    _, gradient = _polynomials(points)
    _, wave_gradient, _, stretch = _wave(points, omega, angle, onset, width)
    tensor, _ = coefficients(stretch)
    if wrong_corner:
        tensor[:, 2] = 1.  # Deliberate missing Jacobian in the vertical term.
    volume = np.einsum('qai,qi,qi,q->a', gradient, wave_gradient, tensor, weights)
    boundary = np.zeros(35, complex)
    surface = np.zeros(35, complex)
    for axis, side in product(range(3), (0, 1)):
        points, weights = _rule(bounds, order, (axis, bounds[axis, side]))
        value, _ = _polynomials(points)
        phi, _, conormal, stretch = _wave(points, omega, angle, onset, width)
        if axis == 2 and side == 1:
            assert bounds[axis, side] == 0.
            _, measure = coefficients(stretch)
            surface += omega**2/9.81*np.einsum('qa,q,q,q->a', value, phi, measure, weights)
        else:
            boundary += (2*side-1)*np.einsum('qa,q,q->a', value, conormal[:, axis], weights)
    scale = np.maximum(1., abs(volume)+abs(surface)+abs(boundary))
    return float(np.max(abs(volume-surface-boundary)/scale))


def test_signed_coordinates_match_derivatives_and_preserve_physical_pocket():
    points = np.array([[-1.1, -.8, -1.], [-.6, .4, -.2], [.1, -.2, -.3], [1.1, .8, -.5]])
    onset, width = [.6, .4], [.8, .7]
    mapped, slope = coordinates(points, onset, width, 8., 4.)
    np.testing.assert_array_equal(mapped[1:3], points[1:3])
    np.testing.assert_array_equal(slope[1:3], np.ones((2, 2)))
    assert np.all(mapped[0, :2].imag < 0) and np.all(mapped[-1, :2].imag > 0)
    for axis in (0, 1):
        step = 1e-5
        plus = points.copy();plus[:, axis] += step
        minus = points.copy();minus[:, axis] -= step
        derivative = (coordinates(plus, onset, width, 8., 4.)[0][:, axis]
                      -coordinates(minus, onset, width, 8., 4.)[0][:, axis])/(2*step)
        np.testing.assert_allclose(derivative, slope[:, axis], atol=3e-9, rtol=3e-9)
    tensor, surface = coefficients(np.ones((2, 2), complex))
    np.testing.assert_array_equal(tensor, np.ones((2, 3)))
    np.testing.assert_array_equal(surface, np.ones(2))


def test_3d_corner_weak_balance_and_independent_conormal():
    # Separate quiet, one-axis layer, and all four corners. Bounds stop at
    # coefficient breakpoints; quadrature does not integrate across a ramp kink.
    regions = [np.array([[-.6, .6], [-.4, .4], [-4., 0.]]),
               np.array([[.6, 1.4], [-.4, .4], [-4., 0.]])]
    for sx, sy in product((-1, 1), repeat=2):
        regions.append(np.array([sorted(sx*np.array([.6, 1.4])),
                                 sorted(sy*np.array([.4, 1.1])), [-4., 0.]]))
    for omega in (2., 4., 6.):
        for bounds in regions:
            angle = np.arctan2(np.mean(bounds[1]), np.mean(bounds[0]))+.17
            coarse = _balance(omega, angle, bounds, 12)
            fine = _balance(omega, angle, bounds, 20)
            assert fine < 2e-11, (omega, bounds, coarse, fine)
            assert fine <= coarse+2e-13
    # A plausible but incorrect corner tensor must fail this physical control.
    assert _balance(4., np.pi/4, regions[-1], 20, wrong_corner=True) > 1e-3
