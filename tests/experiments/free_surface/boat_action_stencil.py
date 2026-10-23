"""Experimental full action RHS using one current gradient and owned offsets.

No operator reuse or material-chart fallback. Both h/h2 controls are retained;
choosing h/2 here is diagnostic and does not establish an accuracy bound.
"""
from dataclasses import dataclass
from time import perf_counter

import numpy as np
from scipy.spatial.transform import Rotation

from topwater import InitialState
from topwater.frames import rotation_matrix
from tests.reference_boat import GRAVITY, WATER_DENSITY
from .boat_ale_pressure import differentiate
from .boat_body_action import body_derivatives, known_load
from .boat_geometry_template import InvalidTemplate
from .boat_parent_action import action_rates, port_transport_covector, port_seam_load
from .body_pressure import evaluate_body_pressure


@dataclass
class ActionResult:
    acceleration: np.ndarray
    hydrodynamic_load: np.ndarray
    diagnostics: dict


def evaluate(body, arrays, current, assembly, state, values, time, geometry, field, *, geometry_tangent=None):
    """Return the moving driver's stencil tuple, with a fresh action closure.

    geometry(pose, eta, psi, label) must stay on the selected base chart;
    field(arrays, time) returns a fresh field and energy on each owned offset.
    """
    from .cuda_shape_gradient import vertex_gradient
    started = perf_counter()
    size = len(arrays['parent_nodes'])
    pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
    twist, eta, psi = state[7:13], state[13:13+size], state[13+size:]
    before = perf_counter()
    time_geometry_evaluations = 0
    gradient = vertex_gradient(arrays['points'], arrays['faces'], current['potential'],
        assembly=assembly, order=assembly.quadrature_order, density=WATER_DENSITY)
    timings = dict(gradient_seconds=perf_counter()-before)
    def fixed_geometry(current_pose, height):
        return geometry(current_pose, height, psi, 'action-shape')
    if geometry_tangent is None:
        before = perf_counter()
        covector = port_transport_covector(arrays, current['potential'], values, time)
        timings['port_covector_seconds'] = perf_counter()-before
        before = perf_counter()
        parent = action_rates(arrays, current, gradient, values, twist, time, fixed_geometry,
            combine_incident_ports=True, port_covector=covector)
        timings['parent_geometry_load_seconds'] = perf_counter()-before
        before = perf_counter()
        derivatives = body_derivatives(arrays, current, gradient, values, twist, time,
            parent['eta_rate'], fixed_geometry, combine_incident_ports=True, port_covector=covector)
        timings['body_geometry_load_seconds'] = perf_counter()-before
        eta_rate, psi_rate = parent['eta_rate'], parent['psi_rate'][1]
        shape_evaluations = parent['geometry_evaluations']+derivatives['geometry_evaluations']
        parent_difference = float(np.max(abs(parent['psi_rate'][0]-parent['psi_rate'][1])))
        condition = parent['parent_mass_condition']
        seam_difference = float(np.max(abs(parent['port_seam_weak_order6']-parent['port_seam_weak'])))
    else:
        from .boat_action_tangent import action_terms
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
        # One exact current pose derivative serves both TIME-offset controls.
        derivatives = dict(pose_load=np.broadcast_to(chain['body_pose_load'], (2, 6)))
        shape_evaluations, parent_difference, condition = 1, None, chain['parent_mass_condition']
        seam_difference = float(np.max(abs(port_seam_load(arrays, current, values, time, order=6)
            -port_seam_load(arrays, current, values, time, order=12))))
    if not np.isfinite(np.r_[eta_rate, psi_rate]).all():
        raise FloatingPointError('nonfinite action surface rates')
    before = perf_counter()
    for halving in range(13):
        h = 1e-5/2**halving
        offsets = {}
        try:
            for dt in (h, -h, h/2, -h/2):
                q = (Rotation.from_rotvec(dt*twist[3:])*Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
                time_geometry_evaluations += 1
                moved = geometry(InitialState(pose.position+dt*twist[:3], q), eta+dt*eta_rate,
                    psi+dt*psi_rate, f'action-time-{dt:g}')
                for key in ('faces', 'free_nodes', 'free_faces', 'body_faces', 'wall_faces', 'bottom_faces', 'node_chart_key'):
                    if not np.array_equal(moved[key], arrays[key]):
                        raise ValueError('action time derivative leaves owned chart/basis')
                offsets[dt] = moved
        except InvalidTemplate:
            continue
        break
    else:
        raise InvalidTemplate('no common action time stencil')
    timings['time_geometry_seconds'] = perf_counter()-before
    before = perf_counter()
    fields, energies = {}, {}
    for dt, moved in offsets.items():
        fields[dt], energies[dt] = field(moved, time+dt)
    timings['offset_field_seconds'] = perf_counter()-before
    rotation = rotation_matrix(pose.orientation)
    inertia = rotation@body.inertia_cg@rotation.T
    rigid = np.zeros((6, 6)); rigid[:3, :3] = body.mass*np.eye(3); rigid[3:, 3:] = inertia
    added = WATER_DENSITY*current['body_load'].T@current['acceleration_potentials']
    reciprocity = np.linalg.norm(added-added.T)/max(np.linalg.norm(added), 1e-30)
    if not np.isfinite(added).all() or reciprocity > 1e-10:
        raise ValueError('invalid action added mass')
    np.linalg.cholesky(rigid+added)
    weight_gyro = np.r_[0., 0., -body.mass*GRAVITY, -np.cross(twist[3:], inertia@twist[3:])]
    loads, accelerations, pressure_loads = [], [], []
    for level, step in enumerate((h, h/2)):
        load, _, _ = known_load(derivatives, current, twist, fields[step], fields[-step], step, level)
        loads.append(load)
        accelerations.append(np.linalg.solve(rigid+added, load+weight_gyro))
        rate, velocity, ale = differentiate(arrays, offsets[step], offsets[-step],
            fields[step]['potential'], fields[-step]['potential'], step, pose.position, twist)
        pressure = evaluate_body_pressure(arrays['points'], arrays['faces'], arrays['body_faces'],
            pose.position, body.mass, inertia, twist, current['potential'], rate, velocity,
            current['acceleration_potentials'], density=WATER_DENSITY, gravity=GRAVITY)
        pressure_loads.append(pressure.known_pressure_load)
    if not np.isfinite(np.r_[loads[1], accelerations[1]]).all():
        raise FloatingPointError('nonfinite action body load/acceleration')
    hydro = loads[1]-added@accelerations[1]
    balance = rigid@accelerations[1]-hydro-weight_gyro
    diagnostics = dict(formulation='experimental restricted action',
        geometry_derivatives='analytic owned chart' if geometry_tangent is not None else 'h/h2 finite differences',
        incident_port_derivatives='Exact cancellation before geometry differencing; full incident forcing remains in fluid solves.',
        added_mass_reciprocity_relative=float(reciprocity), rigid_body_balance_world=balance.tolist(),
        known_load_world=loads[1].tolist(), h_h2_known_load_difference=(loads[0]-loads[1]).tolist(),
        h_h2_acceleration_difference=(accelerations[0]-accelerations[1]).tolist(),
        pressure_known_load_world=pressure_loads[1].tolist(),
        action_minus_pressure_known_load_world=(loads[1]-pressure_loads[1]).tolist(),
        geometry_evaluations=shape_evaluations+time_geometry_evaluations,
        timings=dict(timings, total_seconds=perf_counter()-started))
    surface = dict(formulation='experimental restricted action',
        maximum_h_h2_psi_rate_difference=parent_difference,
        parent_mass_condition=condition, maximum_seam_order6_12_difference=seam_difference)
    result = ActionResult(accelerations[1], hydro, diagnostics)
    return (result, eta_rate, psi_rate, surface, [arrays, offsets[h/2], offsets[-h/2]],
        [assembly.energy, energies[h/2], energies[-h/2]],
        dict(representation='action with owned nodal offsets', step_seconds=h/2, diagnostics=ale))
