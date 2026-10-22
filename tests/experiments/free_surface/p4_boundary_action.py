"""One oriented curved P4 triangle's full-boundary velocity one-form.

F=int phi*U dot (X_u cross X_v), with fixed nodal potential. This applies to
the entire physical boundary, including walls and downward-facing hulls; no
projected-normal sign condition is imposed. Sampled nondegeneracy cannot
certify orientation/injectivity between samples; the geometry chart owns that.
"""
import numpy as np

from .p4_surface import quadrature


def integrals(coordinates, potential, velocity, *, order=9):
    """Return F, partial_X F and partial_U F, in indices(3) face ordering.

    Coordinates and velocity are (15,3); potential is (15,). No density or
    sign of a Lagrangian is included. General P4 products reach degree14;
    order9 Gauss-Duffy covers these polynomial products exactly. All nodal
    quantities are held independent in the returned partial derivatives.

    For X(q), U=J(q)*qdot, the configuration derivative at fixed qdot is
    J.T*partial_X F + (D_qdot J).T*partial_U F, by Hessian symmetry. The
    caller assembles this with its actual geometry map, including the local
    rotation-coordinate/world-angular-velocity conversion.
    """
    coordinates, potential, velocity = [np.asarray(value, dtype=float)
        for value in (coordinates, potential, velocity)]
    if (coordinates.shape != (15, 3) or potential.shape != (15,)
            or velocity.shape != (15, 3)
            or not all(np.all(np.isfinite(value)) for value in (coordinates, potential, velocity))):
        raise ValueError('finite P4 coordinates/velocity(15,3) and potential(15) required')
    values, derivative, weights = quadrature(order)
    tangent = np.einsum('qaj,ai->qji', derivative, coordinates)
    normal_area = np.cross(tangent[:, 0], tangent[:, 1])
    norm = np.linalg.norm(normal_area, axis=1)
    if np.any(norm <= 0) or not np.all(np.isfinite(norm)):
        raise ValueError('degenerate or invalid sampled P4 boundary Jacobian')
    phi = values@potential
    current_velocity = values@velocity
    measure = weights*phi
    velocity_gradient = np.einsum('q,qa,qi->ai', measure, values, normal_area)
    # delta X_u dot (X_v cross U) + delta X_v dot (U cross X_u).
    position_gradient = np.einsum('q,qa,qi->ai', measure, derivative[:, :, 0],
                                  np.cross(tangent[:, 1], current_velocity))
    position_gradient += np.einsum('q,qa,qi->ai', measure, derivative[:, :, 1],
                                   np.cross(current_velocity, tangent[:, 0]))
    return dict(action=float(np.einsum('q,qi,qi->', measure, current_velocity, normal_area)),
        position_gradient=position_gradient, velocity_gradient=velocity_gradient,
        minimum_sampled_jacobian_norm=float(norm.min()))
