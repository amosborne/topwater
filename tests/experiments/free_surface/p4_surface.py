"""One physical P4 graph triangle's surface action, not a canoe solver.

Geometry and potential use p4_tetra.indices(3) ordering. Geometry Jacobians
are supplied by the actual configuration chart; height modes must move the
geometry rather than being added only to a nominal mass matrix. Positive
sampled projected Jacobians do not certify orientation between samples.
"""
from functools import lru_cache

import numpy as np
from tests.laplace_reference import triangle_samples
from .p4_tetra import basis


@lru_cache(None)
def quadrature(order):
    barycentric, weights = triangle_samples(np.eye(3)[None], order)
    values, derivative = basis(barycentric[0])
    # triangle_samples weights sum to one; reference triangle area is 1/2.
    return values, derivative[:, :, 1:]-derivative[:, :, :1], weights/2


def integrals(coordinates, height_jacobians, body_jacobian, *, order=9):
    """Return C_eta, C_v, projected mass and physical boundary moments.

    coordinates: (15,3); height_jacobians: three dense/CSR (15,n_eta)
    matrices; body_jacobian: (15,3,6), with world linear/angular columns
    when the caller's chart represents an actual body. The one-form is
    psi.T@(eta_map@eta_dot + body_transport@world_twist).

    volume_boundary=int z*n_z dS and gravity_boundary=rho*g/2*int z²*n_z
    are this surface's contributions to CLOSED-boundary physical moments.
    A general isolated open patch does not enclose a volume. Density1000,
    gravity9.81. Degree14 general P4 surface products require adequate
    quadrature; order9 Gauss-Duffy integrates these polynomial products.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    height = tuple(np.asarray(matrix.toarray() if hasattr(matrix, 'toarray') else matrix,
                              dtype=float) for matrix in height_jacobians)
    body = np.asarray(body_jacobian, dtype=float)
    if (coordinates.shape != (15, 3) or len(height) != 3
            or any(matrix.ndim != 2 or matrix.shape[0] != 15 for matrix in height)
            or len({matrix.shape[1] for matrix in height}) != 1
            or height[0].shape[1] == 0 or body.shape != (15, 3, 6)
            or not all(np.all(np.isfinite(value)) for value in (coordinates, *height, body))):
        raise ValueError('finite P4 coordinates(15,3) and matching height/body Jacobians required')
    values, derivative, weights = quadrature(order)
    tangent = np.einsum('qaj,ai->qji', derivative, coordinates)
    normal_area = np.cross(tangent[:, 0], tangent[:, 1])
    if np.any(normal_area[:, 2] <= 0) or not np.all(np.isfinite(normal_area)):
        raise ValueError('nonpositive or invalid sampled P4 graph Jacobian')
    # At a quadrature point, each chart velocity is itself interpolated in
    # the geometry basis. Integrate against every potential test function.
    eta_velocity = np.stack([values@matrix for matrix in height], axis=1)
    body_velocity = np.einsum('qa,aik->qik', values, body)
    eta_map = np.einsum('q,qa,qi,qik->ak', weights, values, normal_area, eta_velocity)
    transport = np.einsum('q,qa,qi,qik->ak', weights, values, normal_area, body_velocity)
    mass = np.einsum('q,q,qa,qb->ab', weights, normal_area[:, 2], values, values)
    z = values@coordinates[:, 2]
    return dict(eta_map=eta_map, body_transport=transport, mass=mass,
        volume_boundary=float(weights@(z*normal_area[:, 2])),
        gravity_boundary=float(4905*weights@(z*z*normal_area[:, 2])),
        minimum_sampled_projected_jacobian=float(normal_area[:, 2].min()))
