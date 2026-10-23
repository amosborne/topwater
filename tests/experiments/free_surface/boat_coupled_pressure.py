"""Implicit six-acceleration coupling through one current fluid solve."""
from time import perf_counter
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres

from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.frames import rotation_matrix
from .cuda_graph import cp


def body_dynamics(body, state):
    """World-axis CG mass matrix and gravity/gyroscopic known load."""
    rotation = rotation_matrix(state.orientation)
    inertia = rotation@body.inertia_cg@rotation.T
    omega = rotation@state.angular_velocity
    mass = np.zeros((6, 6))
    mass[:3, :3] = body.mass*np.eye(3)
    mass[3:, 3:] = inertia
    load = np.r_[0., 0., -body.mass*GRAVITY, -np.cross(omega, inertia@omega)]
    return mass, load


def coupled_pressure(context, body, state, normal_rhs, known_phi_t, velocity):
    """Solve A sigma_dot=b_known+W a and M a+C sigma_dot=f_known.

    known_phi_t is the frozen-density Eulerian derivative, excluding the
    current density response. normal_rhs excludes all prescribed acceleration.
    This function owns neither offset geometry construction nor its guards.
    """
    started = perf_counter()
    panels = context.panels
    relative = panels.centroids-state.position
    modes = np.column_stack((panels.normals, np.cross(relative, panels.normals)))
    weighted = panels.areas[:, None]*modes
    mass, external = body_dynamics(body, state)
    known_pressure = -WATER_DENSITY*(known_phi_t+.5*np.sum(velocity**2, axis=1)
                                    +GRAVITY*panels.centroids[:, 2])
    known_load = weighted.T@known_pressure+external
    applications = 0

    def response(value):
        spectrum = -context.spectrum(context.trace@cp.asarray(value))
        potential = context.body_potential(value, spectrum)
        normal = cp.asnumpy(context.matrix@cp.asarray(value)+(context.normal@spectrum).real)
        force = WATER_DENSITY*(weighted.T@potential)
        return potential, normal, force

    def apply(value):
        nonlocal applications
        applications += 1
        _, normal, force = response(value)
        return normal+modes@np.linalg.solve(mass, force)

    operator = LinearOperator((len(panels), len(panels)), matvec=apply, dtype=np.float64)
    right = normal_rhs+modes@np.linalg.solve(mass, known_load)
    history = []
    density_dot, info = gmres(operator, right, M=context.precondition,
        rtol=1e-9, atol=0., restart=min(40, len(panels)), maxiter=5,
        callback=lambda residual: history.append(float(residual)), callback_type='pr_norm')
    potential, normal, force = response(density_dot)
    acceleration = np.linalg.solve(mass, known_load-force)
    field_defect = normal-normal_rhs-modes@acceleration
    relative_residual = float(np.linalg.norm(apply(density_dot)-right)/max(np.linalg.norm(right), 1e-30))
    if info != 0 or not np.isfinite(relative_residual) or relative_residual > 2e-9:
        raise RuntimeError(f'coupled current solve failed: info={info}, residual={relative_residual}')
    pressure = known_pressure-WATER_DENSITY*potential
    load = weighted.T@pressure
    return dict(acceleration=acceleration, density_dot=density_dot, pressure=pressure,
        pressure_functionals=load, mass_matrix=mass,
        diagnostics=dict(total_seconds=perf_counter()-started, gmres_iterations=len(history),
            operator_applications=applications, true_relative_residual=relative_residual,
            body_balance_residual=(mass@acceleration-load-external).tolist(),
            fluid_normal_residual_max=float(np.max(abs(field_defect))),
            additional_operator_builds=0, acceleration_response_basis_builds=0,
            scope='Current-domain implicit algebra only; geometry/state/pressure qualification is separate'))
