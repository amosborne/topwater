"""Unqualified current-geometry P4 canonical stage for the reference canoe.

This composes the same discrete fluid energy, actual geometry surface action
and rigid-body energy. Fixed-chart validity must be established by the geometry
owner before use. No trajectory, flooding/remeshing or Capytaine qualification
is implied by the existence of this stage.
"""
import numpy as np
from scipy.spatial.transform import Rotation

from topwater.frames import rotation_matrix
from .canonical_canoe import body_map
from .p4_canoe_field import prepare, solve
from .p4_geometry import kinetic_integrals as cell_integrals
from .p4_boundary_action import integrals as boundary_integrals
from .p4_tetra import face_nodes


def encode(data, psi, velocity):
    """Physical psi/world twist to momenta using this exact prepared geometry."""
    psi, velocity = np.asarray(psi), np.asarray(velocity)
    if (psi.shape != (len(data['free']),) or velocity.shape != (6,)
            or not np.all(np.isfinite(psi)) or not np.all(np.isfinite(velocity))):
        raise ValueError('finite P4 surface potential and world six-twist required')
    pi = 1000*data['eta_map'].T@psi
    body = data['body_map'].T@(data['effective_mass']@velocity+1000*data['coupling'].T@psi)
    return np.r_[data['q'], body, pi]


def decode(data, momentum):
    """Current physical state and configuration velocity without dense Schur."""
    momentum = np.asarray(momentum)
    if momentum.shape != data['q'].shape or not np.all(np.isfinite(momentum)):
        raise ValueError('finite body and P4 surface canonical momenta required')
    psi = data['eta_factor'].solve(momentum[6:]/1000, trans='T')
    velocity = np.linalg.solve(data['effective_mass'],
        np.linalg.solve(data['body_map'].T, momentum[:6])-1000*data['coupling'].T@psi)
    field = solve(data, psi, velocity)
    qdot = np.r_[np.linalg.solve(data['body_map'], velocity), field['eta_rate']]
    physical = np.r_[data['q'][:3], data['orientation'], velocity, data['q'][6:], psi]
    return dict(field, physical_state=physical, qdot=qdot)


def _configuration_state(geometry, q):
    orientation = (Rotation.from_rotvec(q[3:6])*Rotation.from_quat(
        geometry.domain.rest.orientation, scalar_first=True)).as_quat(scalar_first=True)
    return geometry.evaluate(q[:3], orientation, q[6:])


def _pullback(state, q, force):
    body, height = state.pullback(force)
    return np.r_[body_map(q[3:6]).T@body, height]


def _rigid_lagrangian(geometry, q, qdot):
    domain = geometry.domain
    orientation = (Rotation.from_rotvec(q[3:6])*Rotation.from_quat(
        domain.rest.orientation, scalar_first=True)).as_quat(scalar_first=True)
    rotation = rotation_matrix(orientation)
    velocity = body_map(q[3:6])@qdot[:6]
    return (.5*domain.body.mass*np.dot(velocity[:3], velocity[:3])
            +.5*velocity[3:]@(rotation@domain.body.inertia_cg@rotation.T)@velocity[3:])


def evaluate(geometry, data, momentum, *, derivative_step=1e-4):
    """One full nonlinear stage using a prepared CURRENT configuration.

    The fixed-potential shape derivative uses the stationary unreduced action.
    The only directional geometry difference differentiates J.T at fixed nodal
    force along qdot; Hessian symmetry includes every configuration coordinate.
    Three small rigid rotation differences supply the explicit body-chart term.
    """
    return action_rate(geometry, data, decode(data, momentum), derivative_step=derivative_step)


def action_rate(geometry, data, current, *, derivative_step=1e-4, kinetic=None):
    """Full action gradient from a decoded field, independent of its solver.

    ``kinetic`` optionally supplies (physical kinetic energy, full nodal shape
    gradient) evaluated from the same current field/geometry/quadrature. This
    keeps boundary action, chart pullbacks and rigid dynamics in one place.
    """
    if not np.isfinite(derivative_step) or derivative_step <= 0:
        raise ValueError('finite positive geometry derivative step required')
    volume_order, boundary_order = data['volume_order'], data['boundary_order']
    q, qdot = data['q'], current['qdot']
    points, potential = data['state'].points, current['potential']
    velocity = data['state'].velocity(qdot[6:], current['physical_state'][7:13])
    position_force = np.zeros_like(points)
    velocity_force = np.zeros_like(points)
    boundary_flux = np.zeros(len(points))
    action = 0.
    for nodes in face_nodes(geometry.mesh, geometry.faces):
        term = boundary_integrals(points[nodes], potential[nodes], velocity[nodes], order=boundary_order)
        action += term['action']
        np.add.at(position_force, nodes, term['position_gradient'])
        np.add.at(velocity_force, nodes, term['velocity_gradient'])
        np.add.at(boundary_flux, nodes, term['potential_gradient'])
    if kinetic is None:
        kinetic_gradient = np.zeros_like(points)
        physical_kinetic = 0.
        for nodes in geometry.mesh.elements:
            term = cell_integrals(points[nodes], potential[nodes], order=volume_order)
            physical_kinetic += term['kinetic_energy']
            np.add.at(kinetic_gradient, nodes, term['kinetic_gradient'])
    else:
        physical_kinetic, kinetic_gradient = kinetic
        kinetic_gradient = np.asarray(kinetic_gradient)
        if (not np.isfinite(physical_kinetic) or kinetic_gradient.shape != points.shape
                or not np.all(np.isfinite(kinetic_gradient))):
            raise ValueError('finite physical kinetic energy and full nodal shape force required')
    gravity = data['physical_boundary_moments']
    # The complete body/surface velocity term is int_boundary phi*Xdot.n:
    # body-chart tangential sliding does no normal work, and walls have zero
    # normal velocity. Those are physical chart invariants, not corrections.
    momentum_rate = _pullback(data['state'], q,
        1000*position_force-kinetic_gradient-gravity['gravity_gradient'])
    step = derivative_step
    plus, minus = q+step*qdot, q-step*qdot
    jdot_force = (_pullback(_configuration_state(geometry, plus), plus, velocity_force)
                  -_pullback(_configuration_state(geometry, minus), minus, velocity_force))/(2*step)
    momentum_rate += 1000*jdot_force
    # Explicit q dependence of rigid kinetic energy at fixed coordinate rates.
    rotation_step = 1e-5
    for axis in range(3, 6):
        direction = np.zeros_like(q); direction[axis] = rotation_step
        momentum_rate[axis] += (_rigid_lagrangian(geometry, q+direction, qdot)
                                -_rigid_lagrangian(geometry, q-direction, qdot))/(2*rotation_step)
    momentum_rate[2] -= geometry.domain.body.mass*9.81
    rate = np.r_[qdot, momentum_rate]
    if not np.all(np.isfinite(rate)):
        raise FloatingPointError('nonfinite P4 canonical stage')
    normal_residual = float(np.linalg.norm(data['energy']@potential-boundary_flux))
    normal_scale = float(np.linalg.norm(boundary_flux))
    return dict(current, rate=rate, fluid_volume=data['fluid_volume'],
        full_boundary_velocity_action=action, physical_fluid_kinetic_energy=physical_kinetic,
        full_weak_normal_velocity_residual=normal_residual,
        full_weak_normal_velocity_relative_residual=normal_residual/normal_scale if normal_scale else None,
        fluid_work_identity_absolute_residual=float(abs(1000*action-2*physical_kinetic)))


def rhs(geometry, canonical, *, derivative_step=1e-4, diagnostics=None):
    """Fresh current operator for each call; no lagging or energy projection."""
    canonical = np.asarray(canonical)
    size = 6+len(geometry.free_nodes)
    if canonical.shape != (2*size,) or not np.all(np.isfinite(canonical)):
        raise ValueError('finite P4 canoe canonical state required')
    data = prepare(geometry, canonical[:size])
    result = evaluate(geometry, data, canonical[size:], derivative_step=derivative_step)
    if diagnostics is not None:
        diagnostics.update(result)
    return result['rate']
