"""Prescribed-motion reaction from the same nonlinear variational canoe action.

The controller supplies only body generalized force. Fluid eta/Pi evolve with
the current variational equations; the reaction includes current added mass,
geometry transport and the prescribed acceleration. No frozen coefficients or
separate hydrostatic load is inserted into the fluid/body calculation.
"""
import numpy as np
from scipy.sparse.linalg import splu
from scipy.spatial.transform import Rotation

from topwater.frames import rotation_matrix
from .canonical_canoe import body_map, configuration, encode, _surface_maps
from .moving_canoe import geometry
from .surface_damping import matrices
from .variational_canoe import rhs

RHO = 1000.


def _surface(domain, q):
    orientation = (Rotation.from_rotvec(q[3:6])*Rotation.from_quat(
        domain.rest.orientation, scalar_first=True)).as_quat(scalar_first=True)
    points, height, body = geometry(domain, q[:3], orientation, q[6:])
    eta_map, _ = _surface_maps(domain, points, height, body)
    return points, eta_map


def reaction(domain, body_state, eta, pi, acceleration, *, damping=None, step=1e-4):
    """Return imposed-motion fluid rates and six world-axis hydrodynamic loads.

    ``body_state`` is (CG position, unit wxyz quaternion, WORLD linear/angular
    velocity), length 13. Eta and canonical Pi use domain.free_nodes order;
    acceleration is world linear/angular acceleration. The optional existing
    weak sponge must vanish on the local pocket. It fixes potential gauge.
    ``step`` controls geometry and inverse-momentum directional differences.
    """
    if not np.isfinite(step) or step <= 0:
        raise ValueError('reaction difference step must be finite and positive')
    body_state, eta, pi, acceleration = [np.asarray(value, dtype=float)
        for value in (body_state, eta, pi, acceleration)]
    size = len(domain.free_nodes)
    if (body_state.shape != (13,) or eta.shape != (size,) or pi.shape != (size,)
            or acceleration.shape != (6,) or not all(np.all(np.isfinite(value))
            for value in (body_state, eta, pi, acceleration))):
        raise ValueError('finite body_state(13), eta/Pi(surface nodes), acceleration(6) required')
    physical = np.r_[body_state, eta, np.zeros(size)]
    q = configuration(domain, physical)
    points, eta_map = _surface(domain, q)
    factor = splu(eta_map.tocsc())
    psi = factor.solve(pi/RHO, trans='T')
    physical[13+size:] = psi
    canonical = encode(domain, physical)
    canonical[-size:] = pi
    details = {}
    free = rhs(domain, canonical, diagnostics=details)
    configuration_size = len(q)
    qdot = free[:configuration_size]
    eta_rate = qdot[6:].copy()
    pi_rate = free[configuration_size+6:].copy()

    def map_rate(direction):
        return (_surface(domain, q+step*direction)[1]
            -_surface(domain, q-step*direction)[1])/(2*step)

    # Pi=rho*Ceta.T*psi has no dependence on body velocity/momentum.
    # Recover fluid psi_dot without calculating the free acceleration which
    # will be replaced by the prescribed acceleration.
    psi_rate = factor.solve(pi_rate/RHO-map_rate(qdot).T@psi, trans='T')
    if damping is not None:
        damping = np.asarray(damping)
        mass, weighted = matrices(points, domain.free_faces, domain.free_nodes, damping)
        if np.any(damping[domain.local_free_indices] != 0):
            raise ValueError('prescribed sponge must vanish on the complete local pocket')
        eta_damping = -factor.solve(weighted@eta)
        psi_damping = -splu(mass.tocsc()).solve(weighted@psi)
        direction = np.r_[np.zeros(6), eta_damping]
        # Consistent mass projection may leak into the local patch even when
        # sigma vanishes there. Retain this geometric Pi pushforward term.
        pi_rate += RHO*(eta_map.T@psi_damping+map_rate(direction).T@psi)
        eta_rate += eta_damping
        psi_rate += psi_damping

    velocity = body_state[7:13]
    quat = body_state[3:7]
    quaternion_rate = .5*np.r_[-velocity[3:]@quat[1:],
        quat[0]*velocity[3:]+np.cross(velocity[3:], quat[1:])]
    physical_rate = np.r_[velocity[:3], quaternion_rate, acceleration, eta_rate, psi_rate]
    offsets = []
    for sign in (1, -1):
        moved = physical+sign*step*physical_rate
        moved[3:7] /= np.linalg.norm(moved[3:7])
        offsets.append(encode(domain, moved))
    encoded_rate = (offsets[0]-offsets[1])/(2*step)
    controller = np.linalg.solve(body_map(q[3:6]).T,
        encoded_rate[configuration_size:configuration_size+6]
        -free[configuration_size:configuration_size+6])
    rotation = rotation_matrix(quat)
    inertia = rotation@domain.body.inertia_cg@rotation.T
    rigid_load = np.r_[domain.body.mass*acceleration[:3], inertia@acceleration[3:]
        +np.cross(velocity[3:], inertia@velocity[3:])]
    weight = np.array([0., 0., -domain.body.mass*9.81, 0., 0., 0.])
    hydrodynamic = rigid_load-weight-controller
    if not np.all(np.isfinite(np.r_[hydrodynamic, eta_rate, pi_rate, psi_rate])):
        raise FloatingPointError('nonfinite prescribed canoe reaction')
    return dict(eta_rate=eta_rate, pi_rate=pi_rate, psi_rate=psi_rate,
        hydrodynamic_load=hydrodynamic, controller_load=controller,
        physical_state=physical, physical_rate=physical_rate,
        fluid_momentum_pushforward_error=encoded_rate[-size:]-pi_rate,
        hamiltonian=details['hamiltonian'], fluid_volume=details['fluid_volume'],
        minimum_rim_clearance=details['minimum_rim_clearance'])
