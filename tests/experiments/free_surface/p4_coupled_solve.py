"""Current fluid/body solve without six auxiliary body-response fields.

Contained algebra experiment, not an adopted time integrator. The current
Dirichlet fluid matrix K is positive definite. Eliminating the six rigid
velocities gives a positive rank-six correction to K; it is applied as an
action, never assembled densely. An existing fluid factor can precondition
the solve while the residual always uses the current operator.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse.linalg import LinearOperator, cg


def solve(matrix, load, rigid_mass, effective_momentum, surface_rhs,
          preconditioner, *, initial=None, rtol=1e-11, maxiter=100):
    """Solve K*phi-B*V=f and M*V+rho*B.T*phi=p, with rho=1000.

    ``p`` has already removed the prescribed free-surface momentum term;
    ``f`` is -K_IF*psi. ``preconditioner`` applies an existing SPD approximate
    inverse of K. Its construction and refresh policy belong to the caller.
    One solve replaces construction of all six K^-1*B columns. Both physical
    block residuals are returned so CG convergence alone cannot hide a sign
    or scaling error. This helper does not change fluid geometry or quadrature.
    """
    load, mass, momentum, forcing = [np.asarray(value, dtype=float) for value in
        (load, rigid_mass, effective_momentum, surface_rhs)]
    size = matrix.shape[0]
    if (matrix.shape != (size, size) or load.shape != (size, 6)
            or mass.shape != (6, 6) or momentum.shape != (6,)
            or forcing.shape != (size,) or not all(np.all(np.isfinite(value))
            for value in (load, mass, momentum, forcing))):
        raise ValueError('matching finite current fluid/body blocks required')
    if not np.isfinite(rtol) or not 0 < rtol < 1 or maxiter < 1:
        raise ValueError('positive relative tolerance below one and iteration limit required')
    # Body mass is only 6x6. The fluid correction remains a matrix action.
    mass_factor = cho_factor(mass)
    inverse_mass = lambda value: cho_solve(mass_factor, value)

    def action(potential):
        return matrix@potential+1000*load@inverse_mass(load.T@potential)

    operator = LinearOperator(matrix.shape, matvec=action, dtype=float)
    inverse = LinearOperator(matrix.shape, matvec=preconditioner, dtype=float)
    rhs = forcing+load@inverse_mass(momentum)
    iterations = 0

    def count(_):
        nonlocal iterations
        iterations += 1

    potential, info = cg(operator, rhs, M=inverse, x0=initial, rtol=rtol,
                         atol=0., maxiter=maxiter, callback=count)
    if info != 0:
        raise RuntimeError(f'current coupled fluid/body CG did not converge: info={info}')
    velocity = inverse_mass(momentum-1000*load.T@potential)
    fluid_residual = matrix@potential-load@velocity-forcing
    body_residual = mass@velocity+1000*load.T@potential-momentum
    if not all(np.all(np.isfinite(value)) for value in
               (potential, velocity, fluid_residual, body_residual)):
        raise FloatingPointError('nonfinite current coupled fluid/body solution or residual')
    return dict(potential=potential, velocity=velocity, iterations=iterations,
                coupled_relative_residual=float(np.linalg.norm(action(potential)-rhs)
                    /max(np.linalg.norm(rhs), np.finfo(float).tiny)),
                fluid_residual=fluid_residual, body_residual=body_residual)
