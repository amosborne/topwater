"""Closed-tank prescribed motion from the full nonlinear P4 canoe action.

Correctness baseline only: THREE fresh current fluid factorizations per call,
base plus two sequential inverse-momentum differences. The base fluid LU is
released before either offset is prepared. No sponge, frozen operator, new
time integrator or open-water qualification is supplied by this helper.
"""
import numpy as np
from scipy.spatial.transform import Rotation

from .canonical_canoe import body_map
from .p4_canoe_field import prepare, surface_maps
from .p4_variational_canoe import encode, evaluate, _configuration_state
from .p4_canoe_release import _rim_bound
from .p4_tetra import face_nodes


def reaction(geometry, body_state, eta, pi, acceleration, *, step=1e-4):
    """Return fluid eta/Pi rates and six WORLD hydrodynamic forces/moments.

    body_state=(CG position, unit wxyz quaternion, WORLD linear/angular
    velocity), length13. Acceleration is WORLD linear/angular acceleration.
    eta and canonical Pi follow geometry.free_nodes, including all high modes.
    The caller owns the fixed reference geometry; no current fluid factor is
    returned or retained. ``step`` controls directional geometry/inverse-map
    differences, not a simulation timestep. P4 chart/cell guards and the same
    conservative positive rim-separation bound as the release remain active.
    General flooding/loss-contact/remeshing coverage remains unqualified.
    """
    body_state, eta, pi, acceleration = [np.asarray(value, dtype=float)
        for value in (body_state, eta, pi, acceleration)]
    count = len(geometry.free_nodes)
    if (body_state.shape != (13,) or eta.shape != (count,) or pi.shape != (count,)
            or acceleration.shape != (6,) or not all(np.all(np.isfinite(value))
            for value in (body_state, eta, pi, acceleration))):
        raise ValueError('finite body_state(13), eta/Pi(P4 surface nodes), acceleration(6) required')
    if not np.isfinite(step) or step <= 0:
        raise ValueError('positive finite directional difference step required')
    if abs(np.linalg.norm(body_state[3:7])-1) > 1e-10:
        raise ValueError('body_state requires a unit wxyz quaternion')
    rotation = Rotation.from_quat(body_state[3:7], scalar_first=True)
    rest_rotation = Rotation.from_quat(geometry.domain.rest.orientation, scalar_first=True)
    q = np.r_[body_state[:3], (rotation*rest_rotation.inv()).as_rotvec(), eta]
    size = len(q)
    surface_nodes = face_nodes(geometry.mesh, geometry.free_faces)
    outer = (geometry.mesh.points[:, :2].min(axis=0), geometry.mesh.points[:, :2].max(axis=0))

    # The only base current operator. Surface eta_factor is small and independent
    # of the volume LU; retaining it below does not retain the volume factor.
    data = prepare(geometry, q)
    clearance = _rim_bound(geometry, data['state'].points, q[:3], data['orientation'],
                           surface_nodes, outer)[0]
    eta_factor = data['eta_factor']
    psi = eta_factor.solve(pi/1000, trans='T')
    velocity = body_state[7:13].copy()
    canonical = encode(data, psi, velocity)
    canonical[-count:] = pi
    free = evaluate(geometry, data, canonical[size:], derivative_step=step)
    qdot = free['qdot'].copy()
    eta_rate = qdot[6:].copy()
    pi_rate = free['rate'][size+6:].copy()
    free_body_momentum_rate = free['rate'][size:size+6].copy()

    # Pi=rho*Ceta.T*psi. Derivative Ceta uses ONLY actual geometry and P4
    # boundary integration; there are no implicit fluid solves in this chain.
    plus_map = surface_maps(geometry, _configuration_state(geometry, q+step*qdot),
                            order=data['boundary_order'])[0]
    minus_map = surface_maps(geometry, _configuration_state(geometry, q-step*qdot),
                             order=data['boundary_order'])[0]
    map_rate = (plus_map-minus_map)/(2*step)
    psi_rate = eta_factor.solve(pi_rate/1000-map_rate.T@psi, trans='T')

    # Use the prepared quaternion representative, so q and -q inputs share the
    # same normalized offset path. Omega is world: qdot=.5*(0,Omega)*q.
    quat = data['orientation'].copy()
    quaternion_rate = .5*np.r_[-velocity[3:]@quat[1:],
        quat[0]*velocity[3:]+np.cross(velocity[3:], quat[1:])]
    physical = np.r_[q[:3], quat, velocity, eta, psi]
    physical_rate = np.r_[velocity[:3], quaternion_rate, acceleration, eta_rate, psi_rate]
    rigid_mass = data['rigid_mass'].copy()
    hamiltonian, fluid_volume = free['hamiltonian'], free['fluid_volume']
    normal_residual = free['full_weak_normal_velocity_relative_residual']
    work_residual = free['fluid_work_identity_absolute_residual']
    # Crucial ownership: data is the ONLY reference to base volume LU. The
    # returned field has arrays/scalars, not the prepared data/factor. Drop both
    # before preparing either offset; likewise release each offset immediately.
    del data, free, eta_factor, plus_map, minus_map, map_rate
    offsets = []
    offset_clearances = []
    for sign in (1, -1):
        moved = physical+sign*step*physical_rate
        moved[3:7] /= np.linalg.norm(moved[3:7])
        relative = Rotation.from_quat(moved[3:7], scalar_first=True)*rest_rotation.inv()
        moved_q = np.r_[moved[:3], relative.as_rotvec(), moved[13:13+count]]
        if np.linalg.norm(moved_q[3:6]-q[3:6]) > np.pi:
            raise ValueError('inverse-momentum difference crossed the local rotation chart')
        offset = prepare(geometry, moved_q)
        offset_clearances.append(_rim_bound(geometry, offset['state'].points,
            moved_q[:3], offset['orientation'], surface_nodes, outer)[0])
        offsets.append(encode(offset, moved[13+count:], moved[7:13])[size:].copy())
        del offset
    momentum_rate = (offsets[0]-offsets[1])/(2*step)
    controller = np.linalg.solve(body_map(q[3:6]).T,
                                 momentum_rate[:6]-free_body_momentum_rate)
    inertia = rigid_mass[3:, 3:]
    rigid_load = rigid_mass@acceleration
    rigid_load[3:] += np.cross(velocity[3:], inertia@velocity[3:])
    weight = np.array([0., 0., -geometry.domain.body.mass*9.81, 0., 0., 0.])
    hydrodynamic = rigid_load-weight-controller
    if not np.all(np.isfinite(np.r_[hydrodynamic, controller, eta_rate, pi_rate, psi_rate])):
        raise FloatingPointError('nonfinite P4 prescribed canoe reaction')
    return dict(eta_rate=eta_rate, pi_rate=pi_rate, psi_rate=psi_rate,
        hydrodynamic_load=hydrodynamic, controller_load=controller,
        physical_state=physical, physical_rate=physical_rate,
        fluid_momentum_pushforward_error=momentum_rate[6:]-pi_rate,
        hamiltonian=hamiltonian, fluid_volume=fluid_volume,
        minimum_rim_clearance_lower_bound=clearance,
        minimum_offset_rim_clearance_lower_bound=min(offset_clearances),
        full_weak_normal_velocity_relative_residual=normal_residual,
        fluid_work_identity_absolute_residual=work_residual)
