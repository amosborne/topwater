"""Current restricted-action RHS with analytic geometry and momentum tangents.

Finite-difference geometry/field controls live in independent audit modules.
This closure uses one current operator and its reverse shape gradient.
"""
from dataclasses import dataclass
from time import perf_counter

import numpy as np

from topwater import InitialState
from topwater.frames import rotation_matrix
from topwater.constants import GRAVITY, WATER_DENSITY
from .boat_ale_pressure import normal_slip
from .boat_action_tangent import action_terms
from .boat_incident_tangent import incident_rate
from .boat_parent_action import port_seam_load
from .body_pressure import evaluate_body_pressure


@dataclass
class ActionResult:
    acceleration: np.ndarray
    hydrodynamic_load: np.ndarray
    eta_rate: np.ndarray
    psi_rate: np.ndarray
    surface: dict
    diagnostics: dict
    ale: dict


def evaluate(body, arrays, current, assembly, state, values, time, geometry_tangent):
    """Evaluate the exact-current action closure on its owned semantic chart.

    geometry_tangent() returns the current refined geometry and its coordinate
    and free-trace Jacobians. Its primal arrays must own the supplied operator.
    """
    from .boat_momentum_tangent import field_momentum_rate
    from .cuda_reverse_gradient import vertex_gradient

    started = perf_counter()
    pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
    twist = state[7:13]
    before = perf_counter()
    gradient = vertex_gradient(arrays['points'], arrays['faces'], current['potential'],
        assembly=assembly, order=assembly.quadrature_order, density=WATER_DENSITY)
    timings = dict(gradient_seconds=perf_counter()-before)
    before = perf_counter()
    fine, dx, ds = geometry_tangent()
    for key in ('points', 'faces', 'free_nodes', 'free_parent_indices', 'free_parent_weights'):
        if not np.array_equal(fine[key], arrays[key]):
            raise ValueError('analytic geometry does not own this current operator/trace')
    timings['analytic_geometry_seconds'] = perf_counter()-before
    before = perf_counter()
    chain = action_terms(arrays, current, gradient, values, twist, time, dx, ds)
    timings['analytic_action_terms_seconds'] = perf_counter()-before
    eta_rate, psi_rate = chain['eta_rate'], chain['psi_rate']
    seam_difference = float(np.max(abs(port_seam_load(arrays, current, values, time, order=6)
        -port_seam_load(arrays, current, values, time, order=12))))
    if not np.isfinite(np.r_[eta_rate, psi_rate]).all():
        raise FloatingPointError('nonfinite action surface rates')
    rotation = rotation_matrix(pose.orientation)
    inertia = rotation@body.inertia_cg@rotation.T
    rigid = np.zeros((6, 6)); rigid[:3, :3] = body.mass*np.eye(3); rigid[3:, 3:] = inertia
    added = WATER_DENSITY*current['body_load'].T@current['acceleration_potentials']
    reciprocity = np.linalg.norm(added-added.T)/max(np.linalg.norm(added), 1e-30)
    if not np.isfinite(added).all() or reciprocity > 1e-10:
        raise ValueError('invalid action added mass')
    np.linalg.cholesky(rigid+added)
    weight_gyro = np.r_[0., 0., -body.mass*GRAVITY, -np.cross(twist[3:], inertia@twist[3:])]
    before = perf_counter()
    velocity = np.einsum('nid,d->ni', dx, np.r_[twist, eta_rate])
    incident = incident_rate(arrays, velocity, values, time)
    timings['incident_rate_seconds'] = perf_counter()-before
    exact = field_momentum_rate(arrays, assembly, current, values, twist, time,
        eta_rate, psi_rate, dx, ds, chain['object_derivatives']['body_load_jacobian'], incident,
        body_pose_load=chain['body_pose_load'])
    timings['momentum_tangent_seconds'] = exact['diagnostics']['total_seconds']
    load = exact['known_load']
    acceleration = np.linalg.solve(rigid+added, load+weight_gyro)
    rate, velocity = exact['potential_rate_ale'], exact['mesh_velocity']
    ale = dict(representation='action with exact current momentum tangent',
        diagnostics=normal_slip(arrays, velocity, pose.position, twist))
    pressure = evaluate_body_pressure(arrays['points'], arrays['faces'], arrays['body_faces'],
        pose.position, body.mass, inertia, twist, current['potential'], rate, velocity,
        current['acceleration_potentials'], density=WATER_DENSITY, gravity=GRAVITY)
    if not np.isfinite(np.r_[load, acceleration]).all():
        raise FloatingPointError('nonfinite action body load/acceleration')
    hydro = load-added@acceleration
    balance = rigid@acceleration-hydro-weight_gyro
    diagnostics = dict(formulation='experimental restricted action',
        shape_gradient='reverse', momentum_derivative='exact current',
        geometry_derivatives='analytic owned chart',
        incident_port_derivatives='Exact cancellation before geometry differencing; full incident forcing remains in fluid solves.',
        added_mass_reciprocity_relative=float(reciprocity), rigid_body_balance_world=balance.tolist(),
        known_load_world=load.tolist(), pressure_known_load_world=pressure.known_pressure_load.tolist(),
        action_minus_pressure_known_load_world=(load-pressure.known_pressure_load).tolist(),
        geometry_evaluations=1, momentum_tangent=exact['diagnostics'],
        timings=dict(timings, total_seconds=perf_counter()-started))
    surface = dict(formulation='experimental restricted action',
        parent_mass_condition=chain['parent_mass_condition'],
        maximum_seam_order6_12_difference=seam_difference)
    return ActionResult(acceleration, hydro, eta_rate, psi_rate, surface, diagnostics, ale)
