"""Contained nonlinear P4 cell calculus; no adopted canoe solver.

The primal integrals and their fixed-potential geometry derivatives use the
same quadrature. Positive sampled Jacobians are necessary, not a certificate
against inversion between samples. Curved kinetic integrals require quadrature
refinement because the inverse geometry Jacobian makes them rational.
"""
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix

from .p4_tetra import basis


class InvalidGeometryError(ValueError):
    """A sampled current P4 cell Jacobian is invalid; a trial may be rejected."""


@lru_cache(None)
def quadrature(order):
    x, w = leggauss(order)
    x, w = (x+1)/2, w/2
    u, v, t = np.meshgrid(x, x, x, indexing='ij')
    wu, wv, wt = np.meshgrid(w, w, w, indexing='ij')
    xyz = np.column_stack((u.ravel(), ((1-u)*v).ravel(), ((1-u)*(1-v)*t).ravel()))
    values, derivative = basis(np.column_stack((1-xyz.sum(axis=1), xyz)))
    return values, derivative[:, :, 1:]-derivative[:, :, :1], (wu*wv*wt*(1-u)**2*(1-v)).ravel()


def metrics(coordinates, order):
    """One cell's physical basis gradients and signed integration measure."""
    coordinates = np.asarray(coordinates)
    if coordinates.shape != (35, 3) or not np.all(np.isfinite(coordinates)):
        raise ValueError('finite P4 geometry(35,3) required')
    values, derivative, weights = quadrature(order)
    # J[reference_axis, physical_axis]; grad_x = J^-1 grad_reference.
    jacobian = np.einsum('qaj,ai->qji', derivative, coordinates)
    determinant = np.linalg.det(jacobian)
    if np.any(determinant <= 0) or not np.all(np.isfinite(determinant)):
        raise InvalidGeometryError('nonpositive or invalid sampled P4 geometry Jacobian')
    gradient = np.einsum('qaj,qij->qai', derivative, np.linalg.inv(jacobian))
    return values, gradient, weights*determinant, float(determinant.min())


def element_stiffness(coordinates, *, order=8):
    """Curved-cell primal from the same quadrature as its energy derivative."""
    _, gradient, measure, _ = metrics(coordinates, order)
    return np.einsum('qai,qbi,q->ab', gradient, gradient, measure, optimize=True)


def stiffness(mesh, coordinates, *, order=8):
    """Current curved geometry, with allocated orphan rows retained as zero.

    This contained reference assembles every current cell. It does not freeze
    geometry, reuse a fluid operator or claim qualified runtime.
    """
    coordinates = np.asarray(coordinates)
    if coordinates.shape != mesh.points.shape or not np.all(np.isfinite(coordinates)):
        raise ValueError('finite current coordinates matching the P4 mesh required')
    local = np.empty((len(mesh.elements), 35, 35))
    for i, nodes in enumerate(mesh.elements):
        local[i] = element_stiffness(coordinates[nodes], order=order)
    nodes = mesh.elements
    return coo_matrix((local.ravel(), (np.repeat(nodes, 35, axis=1).ravel(),
                                      np.tile(nodes, (1, 35)).ravel())),
                       shape=(len(coordinates), len(coordinates))).tocsr()


def integrals(coordinates, potential, *, order=8):
    """One isoparametric P4 tetrahedron; gradients hold nodal phi fixed.

    Coordinates follow p4_tetra.indices(4), not corner-first ordering.
    Density is 1000 kg/m3 and gravity 9.81 m/s2, as in the canoe experiments.
    """
    coordinates, potential = np.asarray(coordinates), np.asarray(potential)
    if potential.shape != (35,) or not np.all(np.isfinite(potential)):
        raise ValueError('finite P4 potential(35) required')
    values, gradient, measure, minimum = metrics(coordinates, order)
    velocity = np.einsum('a,qai->qi', potential, gradient)
    squared = np.einsum('qi,qi->q', velocity, velocity)
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
        minimum_sampled_jacobian=minimum)
