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


def body_preconditioner(fluid_inverse, response, effective_mass):
    """Woodbury inverse of a cached fluid/body system in unchanged master IDs.

    ``response`` is cached K^-1*B and ``effective_mass`` is its matching
    M+rho*B.T*K^-1*B. These are setup data, not current-field approximations;
    ``solve`` still applies and checks the current coupled operator.
    """
    mass_factor = cho_factor(effective_mass)

    def inverse(value):
        return fluid_inverse(value)-1000*response@cho_solve(mass_factor, response.T@value)

    return inverse


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


def benchmark(path, output):
    """Seven current-matrix cases with one reusable factor; no new geometry."""
    from hashlib import sha256
    import json
    from pathlib import Path
    from time import perf_counter
    from scipy.sparse.linalg import splu
    from .p4_fluid_ordering import load_input, fresh

    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    matrix, data, manifest = load_input(path)
    load, mass = data['body_load'], data['rigid_mass']
    report = dict(complete=False, qualified=False,
        scope='Current saved-canoe coupled algebra only; no changing operator or time integration',
        input_npz_sha256=manifest['npz_sha256'],
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        rtol=1e-11, ordering='MMD_AT_PLUS_A', cases=[])
    fields = {}
    failure = None
    try:
        start = perf_counter()
        factor = splu(matrix, permc_spec='MMD_AT_PLUS_A')
        report['preconditioner_setup_wall_seconds'] = perf_counter()-start
        cases = [('saved', data['p_effective'], data['fluid_surface_rhs'],
                  data['reference_velocity'], data['reference_potential'])]
        for column in range(6):
            velocity = np.eye(6)[column]
            cases.append((f'body_axis_{column}', (mass+data['added_mass'])@velocity,
                          np.zeros(len(load)), velocity, data['response'][:, column]))
        for name, momentum, forcing, reference_v, reference_phi in cases:
            start = perf_counter()
            current = solve(matrix, load, mass, momentum, forcing, factor.solve,
                            rtol=report['rtol'])
            elapsed = perf_counter()-start
            potential, velocity = current['potential'], current['velocity']
            fluid_scale = np.linalg.norm(matrix@potential)+np.linalg.norm(load@velocity)+np.linalg.norm(forcing)
            body_scale = np.linalg.norm(mass@velocity)+1000*np.linalg.norm(load.T@potential)+np.linalg.norm(momentum)
            report['cases'].append(dict(name=name, wall_seconds=elapsed, iterations=current['iterations'],
                condensed_relative_residual=current['coupled_relative_residual'],
                fluid_scaled_residual=float(np.linalg.norm(current['fluid_residual'])/max(fluid_scale, 1e-300)),
                body_scaled_residual=float(np.linalg.norm(current['body_residual'])/max(body_scale, 1e-300)),
                velocity_maximum_absolute_difference=float(np.max(abs(velocity-reference_v))),
                potential_relative_difference=float(np.linalg.norm(potential-reference_phi)
                    /max(np.linalg.norm(reference_phi), 1e-300))))
            fields[name+'_potential'] = potential
            fields[name+'_velocity'] = velocity
            print(json.dumps(report['cases'][-1]), flush=True)
        report['complete'] = True
    except BaseException as error:
        failure = error
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        try:
            np.savez(output.with_suffix('.npz'), **fields)
            output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        except BaseException as save_error:
            if failure is not None:
                failure.add_note('Coupled-solve evidence could not be saved: '+repr(save_error))
            else:
                raise
    return report


if __name__ == '__main__':
    import argparse
    from threadpoolctl import threadpool_limits
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        benchmark(args.input, args.output)
