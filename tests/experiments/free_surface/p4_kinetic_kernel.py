"""Same-quadrature kinetic-only cell candidate; reference remains p4_geometry.

No geometry, quadrature, Jacobian guard or physical derivative is changed.
This API computes only the quantities used by the nonlinear kinetic action.
"""
import numpy as np

from .p4_geometry import quadrature, InvalidGeometryError


def metrics(coordinates, order):
    """Same physical metrics, using batched small matrix products."""
    coordinates = np.asarray(coordinates)
    if coordinates.shape != (35, 3) or not np.all(np.isfinite(coordinates)):
        raise ValueError('finite P4 geometry(35,3) required')
    values, derivative, weights = quadrature(order)
    jacobian = derivative.transpose(0, 2, 1)@coordinates
    determinant = np.linalg.det(jacobian)
    if np.any(determinant <= 0) or not np.all(np.isfinite(determinant)):
        raise InvalidGeometryError('nonpositive or invalid sampled P4 geometry Jacobian')
    gradient = derivative@np.linalg.inv(jacobian).transpose(0, 2, 1)
    return values, gradient, weights*determinant, float(determinant.min())


def integrals(coordinates, potential, *, order=8):
    """Physical kinetic energy and full fixed-nodal-potential shape derivative."""
    potential = np.asarray(potential)
    if potential.shape != (35,) or not np.all(np.isfinite(potential)):
        raise ValueError('finite P4 potential(35) required')
    _, gradient, measure, minimum = metrics(coordinates, order)
    velocity = np.einsum('a,qai->qi', potential, gradient)
    squared = np.einsum('qi,qi->q', velocity, velocity)
    # sigma_ij grad_j(Na) = .5|u|² grad_i(Na) - u_i(u.grad(Na)).
    # Contract directly instead of materializing/contracting q x3x3 stresses.
    directional = np.einsum('qai,qi->qa', gradient, velocity)
    force = (gradient.reshape(len(measure), -1).T@(.5*measure*squared)).reshape(35, 3)
    force -= (directional*measure[:, None]).T@velocity
    return dict(kinetic_energy=float(500*measure@squared),
                kinetic_gradient=1000*force, minimum_sampled_jacobian=minimum)


def element_stiffness(coordinates, *, order=8):
    """Unchanged primal contraction, with the same candidate physical metrics."""
    _, gradient, measure, _ = metrics(coordinates, order)
    return np.einsum('qai,qbi,q->ab', gradient, gradient, measure, optimize=True)
