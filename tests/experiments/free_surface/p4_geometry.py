"""Contained nonlinear P4 cell calculus; no adopted canoe solver.

The primal integrals and their fixed-potential geometry derivatives use the
same quadrature. Positive sampled Jacobians are necessary, not a certificate
against inversion between samples. Curved kinetic integrals require quadrature
refinement because the inverse geometry Jacobian makes them rational.
"""
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss

from .p4_tetra import basis


@lru_cache(None)
def quadrature(order):
    x, w = leggauss(order)
    x, w = (x+1)/2, w/2
    u, v, t = np.meshgrid(x, x, x, indexing='ij')
    wu, wv, wt = np.meshgrid(w, w, w, indexing='ij')
    xyz = np.column_stack((u.ravel(), ((1-u)*v).ravel(), ((1-u)*(1-v)*t).ravel()))
    values, derivative = basis(np.column_stack((1-xyz.sum(axis=1), xyz)))
    return values, derivative[:, :, 1:]-derivative[:, :, :1], (wu*wv*wt*(1-u)**2*(1-v)).ravel()


def integrals(coordinates, potential, *, order=8):
    """One isoparametric P4 tetrahedron; gradients hold nodal phi fixed.

    Coordinates follow p4_tetra.indices(4), not corner-first ordering.
    Density is 1000 kg/m3 and gravity 9.81 m/s2, as in the canoe experiments.
    """
    coordinates, potential = np.asarray(coordinates), np.asarray(potential)
    if (coordinates.shape != (35, 3) or potential.shape != (35,)
            or not np.all(np.isfinite(coordinates)) or not np.all(np.isfinite(potential))):
        raise ValueError('finite P4 geometry(35,3) and potential(35) required')
    values, derivative, weights = quadrature(order)
    # J[reference_axis, physical_axis]; grad_x = J^-1 grad_reference.
    jacobian = np.einsum('qaj,ai->qji', derivative, coordinates)
    determinant = np.linalg.det(jacobian)
    if np.any(determinant <= 0) or not np.all(np.isfinite(determinant)):
        raise ValueError('nonpositive or invalid sampled P4 geometry Jacobian')
    gradient = np.einsum('qaj,qij->qai', derivative, np.linalg.inv(jacobian))
    velocity = np.einsum('a,qai->qi', potential, gradient)
    squared = np.einsum('qi,qi->q', velocity, velocity)
    measure = weights*determinant
    height = values@coordinates[:, 2]
    stress = .5*squared[:, None, None]*np.eye(3)-velocity[:, :, None]*velocity[:, None, :]
    kinetic_gradient = 1000*np.einsum('q,qij,qaj->ai', measure, stress, gradient)
    volume_gradient = np.einsum('q,qai->ai', measure, gradient)
    gravity_gradient = np.einsum('q,qai->ai', measure*height, gradient)
    gravity_gradient[:, 2] += values.T@measure
    return dict(kinetic_energy=float(500*measure@squared),
        volume=float(measure.sum()), gravity=float(9810*measure@height),
        kinetic_gradient=kinetic_gradient, volume_gradient=volume_gradient,
        gravity_gradient=9810*gravity_gradient,
        minimum_sampled_jacobian=float(determinant.min()))
