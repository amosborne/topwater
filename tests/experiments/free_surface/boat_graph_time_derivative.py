"""Moving-hull pressure/ALE derivative gate, not a coupled trajectory.

Prescribed six-DOF motion clips the input hull against one fixed nonflat graph's
tangent plane. Only exterior Dirichlet data are supplied; hidden TOTAL trace is
filled algebraically, then the FULL source trace is subtracted. The graph and
its imposed potential are a forced boundary, not an evolved free surface.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation

from topwater import InitialState, SimulationInputs
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.frames import rotation_matrix
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from .boat_graph_bvp import GraphBody
from .boat_graph_fields import body_velocity, actual_graph_values, actual_field_values, source_velocity_bulk
from .boat_transport_pressure import prescribed_pressure
from .boat_graph_images import image_height
from .boat_hos_continuation import tilted_geometry, tangent_graph
from .boat_spatial_work import subdivide
from .boat_parent_clip_proof import coplanar_groups
from .cuda_graph import ExteriorTrace, body_maps, harmonic_velocity, cp
from .panel_gradient_reference import checked_influence
from .graph_geometry_reference import fourier_points
from .boat_frozen_field import frozen_field
from .boat_coupled_pressure import coupled_pressure, body_dynamics


def prescribed_pose(initial, time):
    """Fixed world rotation axis; all rates agree exactly with the pose curve."""
    initial_rotation = Rotation.from_quat(initial.orientation, scalar_first=True)
    omega = initial_rotation.apply(initial.angular_velocity)
    path_time = time+.15*time*time
    orientation = Rotation.from_rotvec(omega*path_time)*initial_rotation
    return InitialState(initial.position+initial.linear_velocity*path_time,
        orientation.as_quat(scalar_first=True), initial.linear_velocity*(1+.3*time),
        orientation.inv().apply(omega*(1+.3*time)))


def constant_twist_pose(initial, time):
    """Zero world linear/angular acceleration for the coupled known RHS."""
    initial_rotation = Rotation.from_quat(initial.orientation, scalar_first=True)
    omega = initial_rotation.apply(initial.angular_velocity)
    orientation = Rotation.from_rotvec(omega*time)*initial_rotation
    return InitialState(initial.position+initial.linear_velocity*time,
        orientation.as_quat(scalar_first=True), initial.linear_velocity,
        orientation.inv().apply(omega))


def wet_snapshot(inputs, pose, plane_origin, frame, grid_xy, length, refinement):
    orientation = Rotation.from_matrix(frame)*Rotation.from_quat(pose.orientation, scalar_first=True)
    state = InitialState(frame@(pose.position-plane_origin), orientation.as_quat(scalar_first=True))
    wet = hydrostatics(inputs.body, state, density=WATER_DENSITY, gravity=GRAVITY)
    triangles = wet.wetted_triangles@frame+plane_origin
    for _ in range(refinement):
        triangles = subdivide(triangles)
    footprint = ConvexHull((wet.waterline_segments@frame+plane_origin)[..., :2].reshape(-1, 2))
    distance = grid_xy@footprint.equations[:, :2].T+footprint.equations[:, 2]
    inside = np.all(distance <= 64*np.finfo(float).eps*length, axis=1)
    relative = (inputs.body.mesh.vertices-inputs.body.cg)@rotation_matrix(pose.orientation).T
    signs = np.sign(((relative+pose.position-plane_origin)@frame.T)[:, 2])
    if np.any(signs == 0):
        raise ValueError('stencil needs persistent clipping recipes, no on-plane input vertex')
    return TrianglePanels(triangles, -1), inside, wet.source_faces, signs


class BulkField:
    def __init__(self, context, points, length, origin):
        shifted = points.copy()
        shifted[:, 2] -= context.reference_height
        size = context.graph.padded_shape[0]
        self.potential, unused = body_maps(shifted, np.zeros_like(shifted), size, length/size, origin)
        del unused
        self.size, self.length = size, length
        self.single = np.empty((len(points), len(context.panels)))
        self.gradient = np.empty((*self.single.shape, 3))
        for first in range(0, len(points), 128):
            last = min(first+128, len(points))
            reflected = points[first:last].copy()
            reflected[:, 2] = 2*context.image_height-reflected[:, 2]
            real = checked_influence(points[first:last], context.panels)
            image = checked_influence(reflected, context.panels)
            self.single[first:last] = real.single_layer-image.single_layer
            self.gradient[first:last] = (real.gradient_single_layer
                -image.gradient_single_layer*np.array([1., 1., -1.]))

    def value(self, density, spectrum):
        return cp.asnumpy((self.potential@spectrum).real)+self.single@density

    def velocity(self, density, spectrum):
        return cp.asnumpy(harmonic_velocity(self.potential, spectrum, self.size, self.length)) \
            +np.einsum('tni,n->ti', self.gradient, density)


def run(inputs, args, report):
    started = perf_counter()
    mesh = inputs.body.mesh
    coplanar_groups(mesh.vertices, mesh.faces, mesh.artificial_faces)
    initial = inputs.initial_state
    cg = initial.position
    _, level, slope = tilted_geometry(inputs)
    normal = np.r_[-slope, 1.]
    normal /= np.linalg.norm(normal)
    tangent = np.array([1., 0., slope[0]])
    tangent /= np.linalg.norm(tangent)
    frame = np.vstack((tangent, np.cross(normal, tangent), normal))
    plane_origin = np.r_[cg[:2], level]
    origin = cg[:2]-args.length/2
    x, y = np.meshgrid(origin[0]+np.arange(args.grid)*args.length/args.grid,
                       origin[1]+np.arange(args.grid)*args.length/args.grid)
    grid_xy = np.column_stack((x.ravel(), y.ravel()))
    eta, _, outer = tangent_graph(inputs, x, y, level, slope)
    if np.any(outer >= args.length/2):
        raise ValueError('fixed graph support must fit strictly inside the box')
    times = [0., args.step, -args.step, args.step/2, -args.step/2, args.step/4, -args.step/4]
    snapshots = {}
    for time in times:
        pose = prescribed_pose(initial, time)
        snapshots[time] = (pose, *wet_snapshot(inputs, pose, plane_origin, frame, grid_xy,
                                              args.length, args.body_refinement))
    _, center_panels, inside, source, signs = snapshots[0.]
    for pose, panels, mask, owners, pattern in snapshots.values():
        if (len(panels) != len(center_panels) or not np.array_equal(owners, source)
                or not np.array_equal(pattern, signs)):
            raise ValueError('density DOF clipping recipes changed across the stencil')
        if not np.array_equal(mask, inside):
            raise ValueError('physical exterior labels changed; state transfer is not implemented')
    extension = ExteriorTrace(inside.reshape(args.grid, args.grid))
    trace_nodes = extension.exterior[np.argsort(np.linalg.norm(
        grid_xy[extension.exterior]-cg[:2], axis=1))[:16]]
    shifts = args.length/args.grid*np.array([[.5, 0.], [-.5, 0.], [0., .5], [0., -.5]])
    off_xy = np.unique(np.concatenate([grid_xy[trace_nodes]+shift for shift in shifts]), axis=0)
    off_height = fourier_points(eta, off_xy, args.length, origin)[0]
    off_points = np.column_stack((off_xy, off_height))
    valid = np.ones(len(off_points), dtype=bool)
    for pose, *_ in snapshots.values():
        hull = ConvexHull((inputs.body.mesh.vertices-inputs.body.cg)@rotation_matrix(pose.orientation).T
                          +pose.position)
        valid &= np.max(off_points@hull.equations[:, :3].T+hull.equations[:, 3], axis=1) > 64*np.finfo(float).eps*args.length
    off_points = off_points[valid]
    off_points = off_points[np.argsort(np.linalg.norm(off_points[:, :2]-cg[:2], axis=1))[:16]]
    if not len(off_points):
        raise ValueError('no safe physical-exterior half-cell trace probes across the stencil')
    image = max(image_height(p.triangles, eta, [args.length]*2, origin)[0]
                for _, p, *_ in snapshots.values())
    true_image_bound = max(float(np.max(p.triangles[..., 2]+level
        +(p.triangles[..., :2]-cg[:2])@slope))/2 for _, p, *_ in snapshots.values())
    image = max(image, float(np.nextafter(true_image_bound, np.inf)))
    minimum_image = image
    if args.image_height is not None:
        if not np.isfinite(args.image_height) or args.image_height < minimum_image:
            raise ValueError('fixed image height must satisfy every stencil image/graph bound')
        image = args.image_height
    height = max(float(p.triangles[..., 2].max()) for _, p, *_ in snapshots.values())
    depth = level+(center_panels.centroids[:, :2]-cg[:2])@slope-center_panels.centroids[:, 2]
    beam = float(np.ptp(inputs.body.mesh.vertices[:, 1]))
    offset = np.minimum(.03*beam, .25*depth/(1+np.linalg.norm(slope)))
    if np.any(offset <= 0):
        raise ValueError('strictly submerged center probes required')
    probes = np.concatenate([center_panels.centroids-fraction*offset[:, None]*center_panels.normals
                             for fraction in [1., .5]])
    material = (probes-cg)@rotation_matrix(initial.orientation)
    height = max(height, float(probes[:, 2].max()), max(float((pose.position
        +material@rotation_matrix(pose.orientation).T)[:, 2].max())
        for pose, *_ in snapshots.values()))
    minimum_height = height
    if args.reference_height is not None:
        if not np.isfinite(args.reference_height) or args.reference_height < minimum_height:
            raise ValueError('fixed expansion height must cover every stencil source and target')
        height = args.reference_height
    omega = rotation_matrix(initial.orientation)@initial.angular_velocity
    probe_speed = initial.linear_velocity+np.cross(omega, probes-cg)
    k = 2*np.pi*np.array([3., 1.])/args.length
    wave_frequency = float(np.sqrt(GRAVITY*np.linalg.norm(k)))
    report.update(input_sha256=sha256(inputs.to_json().encode()).hexdigest(),
        body_panels=len(center_panels), grid_size=args.grid, box_length_m=args.length,
        hos_order=args.order, body_refinement=args.body_refinement,
        image_height_m=image, reference_height_m=height, probe_count=len(probes),
        minimum_sufficient_image_height_m=minimum_image,
        minimum_sufficient_reference_height_m=minimum_height,
        minimum_image_gap_above_true_plane_m=min(float(np.min(2*image-p.triangles[..., 2]
            -level-(p.triangles[..., :2]-cg[:2])@slope)) for _, p, *_ in snapshots.values()),
        exterior_nodes=len(extension.exterior), hidden_nodes=len(extension.hidden),
        frozen_density_source_faces=source.tolist(), frozen_vertex_plane_signs=signs.tolist(),
        imposed_periodic_trace=dict(amplitude_m2_per_s=.1, wavevector_per_m=k.tolist(),
                                    angular_frequency_per_s=wave_frequency),
        probe_offsets_m=offset.tolist())
    center = None
    records = {}
    for time in times:
        pose, panels, *_ = snapshots[time]
        context = GraphBody(panels, eta, args.length, origin, height, image, args.order)
        ext = .1*np.cos(grid_xy[extension.exterior]@k-wave_frequency*time)
        trace = extension.extend(ext)
        angular = rotation_matrix(pose.orientation)@pose.angular_velocity
        speed = pose.linear_velocity+np.cross(angular, panels.centroids-pose.position)
        right = context.right_hand_side(trace, np.sum(speed*panels.normals, axis=1))
        density, diagnostics = context.solve_density(right)
        spectrum = context.spectrum(trace-context.trace@cp.asarray(density))
        body = context.body_potential(density, spectrum)
        moving = pose.position+material@rotation_matrix(pose.orientation).T
        # Both fixed and moving targets must remain strictly in physical fluid.
        hull = ConvexHull((inputs.body.mesh.vertices-inputs.body.cg)@rotation_matrix(pose.orientation).T
                          +pose.position)
        for points in [probes, moving]:
            if (np.any(np.max(points@hull.equations[:, :3].T+hull.equations[:, 3], axis=1) <= 0)
                    or np.any(points[:, 2] >= level+(points[:, :2]-cg[:2])@slope)):
                raise ValueError('bulk stencil target crossed hull or true graph')
        targets = np.vstack((panels.centroids, probes, moving))
        plane_targets = targets.copy()
        plane_targets[:, 2] = 0.
        eta_map, unused = body_maps(plane_targets, np.zeros_like(plane_targets),
                                    args.grid, args.length/args.grid, origin)
        represented_eta = cp.asnumpy((eta_map@cp.fft.fft2(cp.asarray(eta)).ravel()).real)
        represented_depth = float(np.min(represented_eta-targets[:, 2]))
        del eta_map, unused
        if represented_depth <= 0:
            raise ValueError('bulk or body-centroid target crossed represented Fourier graph')
        fixed_field = BulkField(context, probes, args.length, origin)
        fixed = fixed_field.value(density, spectrum)
        actual_trace, maximum_trace_exponent = actual_graph_values(context, density, spectrum, trace_nodes)
        if time == 0:
            center, center_solution = context, (density, spectrum)
            center_field = fixed_field
            velocity = fixed_field.velocity(density, spectrum)
            record = dict(right=right, fixed=fixed, moving=fixed, body=body,
                          frozen_fixed=fixed, frozen_moving=fixed, frozen_body=body)
        else:
            frozen_spectrum = context.spectrum(trace-context.trace@cp.asarray(center_solution[0]))
            frozen_fixed = fixed_field.value(center_solution[0], frozen_spectrum)
            del fixed_field
            moving_field = BulkField(context, moving, args.length, origin)
            record = dict(right=right, fixed=fixed,
                body=body, frozen_body=context.body_potential(center_solution[0], frozen_spectrum),
                moving=moving_field.value(density, spectrum),
                frozen_fixed=frozen_fixed,
                frozen_moving=moving_field.value(center_solution[0], frozen_spectrum),
                operator_on_center=context.apply(center_solution[0]))
            del moving_field
            if abs(time) == args.step/4:
                action = frozen_field(panels, eta, args.length, origin, height, image,
                                      args.order, trace, center_solution[0])
                action['residual'] = np.sum(speed*panels.normals, axis=1)-action['normal']
                action['diagnostics'].update(
                    full_context_residual_max_difference=float(np.max(abs(action['residual']
                        -(right-record['operator_on_center'])))),
                    full_context_potential_max_difference=float(np.max(abs(action['potential']
                        -record['frozen_body']))))
                record['frozen_action'] = action
        records[time] = record
        record['actual_graph_values'] = actual_trace
        record['desired_graph_values'] = cp.asnumpy(trace[trace_nodes])
        record['actual_trace_defect'] = actual_trace-record['desired_graph_values']
        off_actual, _ = actual_field_values(context, density, spectrum, off_points)
        off_desired = fourier_points(cp.asnumpy(trace).reshape(context.shape), off_points[:, :2],
                                    args.length, origin)[0]
        record['off_actual'] = off_actual
        record['off_desired'] = off_desired
        record['off_defect'] = off_actual-off_desired
        report['snapshots'].append(dict(time_seconds=time, **diagnostics,
            context_seconds=context.setup_seconds,
            minimum_target_depth_below_represented_graph_m=represented_depth,
            maximum_actual_trace_raw_exponent=maximum_trace_exponent,
            panel_sha256=sha256(panels.triangles.tobytes()).hexdigest()))
        print(json.dumps(report['snapshots'][-1]), flush=True)
        if time != 0:
            del context
    convection = np.sum(probe_speed*velocity, axis=1)
    surface_velocity = body_velocity(center, *center_solution, args.length)
    body_speed = initial.linear_velocity+np.cross(omega, center_panels.centroids-cg)
    normal_error = np.sum((surface_velocity-body_speed)*center_panels.normals, axis=1)
    weighted = center_panels.areas[:, None]*np.column_stack((center_panels.normals,
        np.cross(center_panels.centroids-cg, center_panels.normals)))
    bulk_errors = [float(np.max(np.linalg.norm(
            velocity[first:first+len(center_panels)]-surface_velocity, axis=1)))
            for first in [0, len(center_panels)]]
    report['surface_trace'] = dict(normal_boundary_condition_max_error_m_per_s=float(np.max(abs(normal_error))))
    limit_rows = [dict(offset_fraction=fraction, velocity_max_error_m_per_s=error)
        for fraction, error in zip([1., .5],
            bulk_errors, strict=True)]
    for fraction in [.25, .125, .0625]:
        points = center_panels.centroids-fraction*offset[:, None]*center_panels.normals
        field = BulkField(center, points, args.length, origin)
        error = np.linalg.norm(field.velocity(*center_solution)-surface_velocity, axis=1)
        limit_rows.append(dict(offset_fraction=fraction, velocity_max_error_m_per_s=float(error.max())))
        del field
    report['surface_trace']['one_sided_bulk_limit'] = limit_rows
    chi_exterior = .1*wave_frequency*np.sin(grid_xy[extension.exterior]@k)
    direct_pressure = prescribed_pressure(center, *center_solution, initial,
        np.r_[.3*initial.linear_velocity, .3*omega], extension, chi_exterior)
    report['current_domain_pressure'] = dict(**direct_pressure['diagnostics'],
        body_pressure_pa=direct_pressure['pressure'].tolist(),
        pressure_functionals=direct_pressure['pressure_functionals'].tolist())
    transport_solution = direct_pressure['transport_solution']
    actual_transport, _ = actual_graph_values(center, transport_solution.density, transport_solution.spectrum, trace_nodes)
    actual_chi = actual_transport-direct_pressure['known_source_transport'][trace_nodes]
    desired_chi = .1*wave_frequency*np.sin(grid_xy[trace_nodes]@k)
    acceleration_trace_defect = actual_chi-desired_chi
    report['actual_graph_trace'] = dict(node_indices=trace_nodes.tolist(),
        positions_world_m=center.surface_points[trace_nodes].tolist(),
        current_flow_trace_defect=records[0.]['actual_trace_defect'].tolist(),
        current_acceleration_trace_defect=acceleration_trace_defect.tolist(),
        scope='Actual finite harmonic sum plus analytic layer trace at fixed exterior grid nodes; no formal Taylor surface evaluation')
    off_transport, _ = actual_field_values(center, transport_solution.density, transport_solution.spectrum, off_points)
    off_source_velocity = source_velocity_bulk(center, center_solution[0], off_points)
    off_rigid = initial.linear_velocity+np.cross(omega, off_points-cg)
    off_chi = off_transport-np.sum(off_rigid*off_source_velocity, axis=1)
    off_desired_chi = fourier_points(cp.asnumpy(extension.extend(chi_exterior)).reshape(center.shape),
                                    off_points[:, :2], args.length, origin)[0]
    off_acceleration_defect = off_chi-off_desired_chi
    report['off_grid_graph_trace'] = dict(positions_world_m=off_points.tolist(),
        current_flow_defect=records[0.]['off_defect'].tolist(),
        current_acceleration_defect=off_acceleration_defect.tolist(),
        prescribed_data_interpolation_error=(records[0.]['off_desired']-.1*np.cos(off_points[:, :2]@k)).tolist(),
        scope='Half-cell probes strictly outside the full hull throughout the stencil; desired state is its owned full-trace Fourier interpolant')
    for step in [args.step, args.step/2, args.step/4]:
        plus, minus = records[step], records[-step]
        defect_rate = (plus['actual_trace_defect']-minus['actual_trace_defect'])/(2*step)
        graph_fd_chi = (plus['actual_graph_values']-minus['actual_graph_values'])/(2*step)
        bdot = (plus['right']-minus['right'])/(2*step)
        adot_sigma = (plus['operator_on_center']-minus['operator_on_center'])/(2*step)
        density_dot, diagnostics = center.solve_density(bdot-adot_sigma)
        response_spectrum = -center.spectrum(center.trace@cp.asarray(density_dot))
        response = center_field.value(density_dot, response_spectrum)
        body_response = center.body_potential(density_dot, response_spectrum)
        implicit_fixed = (plus['frozen_fixed']-minus['frozen_fixed'])/(2*step)+response
        implicit_material = (plus['frozen_moving']-minus['frozen_moving'])/(2*step)+response
        resolved_fixed = (plus['fixed']-minus['fixed'])/(2*step)
        resolved_material = (plus['moving']-minus['moving'])/(2*step)
        body_chart_rate = (plus['frozen_body']-minus['frozen_body'])/(2*step)+body_response
        resolved_body_rate = (plus['body']-minus['body'])/(2*step)
        chart_velocity = (snapshots[step][1].centroids-snapshots[-step][1].centroids)/(2*step)
        chart_convection = np.sum(chart_velocity*surface_velocity, axis=1)
        body_phi_t = body_chart_rate-chart_convection
        body_pressure = -WATER_DENSITY*(body_phi_t+.5*np.sum(surface_velocity**2, axis=1)
                                      +GRAVITY*center_panels.centroids[:, 2])
        resolved_body_pressure = -WATER_DENSITY*(resolved_body_rate-chart_convection
            +.5*np.sum(surface_velocity**2, axis=1)+GRAVITY*center_panels.centroids[:, 2])
        incomplete_density, _ = center.solve_density(bdot)
        incomplete_spectrum = -center.spectrum(center.trace@cp.asarray(incomplete_density))
        incomplete_response = center_field.value(incomplete_density, incomplete_spectrum)
        incomplete_body_response = center.body_potential(incomplete_density, incomplete_spectrum)
        bernoulli = -WATER_DENSITY*(implicit_material-convection+.5*np.sum(velocity*velocity, axis=1)
                                   +GRAVITY*probes[:, 2])
        dynamic_pressure_difference = (bernoulli+WATER_DENSITY*GRAVITY*probes[:, 2]).reshape(2, -1) \
            -body_pressure-WATER_DENSITY*GRAVITY*center_panels.centroids[:, 2]
        row = dict(step_seconds=step, density_derivative_solve=diagnostics,
            implicit_vs_resolved_fixed_phi_t_max_m2_per_s2=float(np.max(abs(implicit_fixed-resolved_fixed))),
            implicit_vs_resolved_material_phi_rate_max_m2_per_s2=float(np.max(abs(implicit_material-resolved_material))),
            resolved_ale_vs_fixed_phi_t_max_m2_per_s2=float(np.max(abs(resolved_material-convection-resolved_fixed))),
            implicit_ale_vs_fixed_phi_t_max_m2_per_s2=float(np.max(abs(implicit_material-convection-implicit_fixed))),
            omitting_operator_derivative_phi_t_max_change_m2_per_s2=float(np.max(abs(incomplete_response-response))),
            omitting_ale_correction_pressure_max_change_pa=float(WATER_DENSITY*np.max(abs(convection))),
            bulk_pressure_min_max_pa=[float(bernoulli.min()), float(bernoulli.max())],
            surface_chart_normal_velocity_max_error_m_per_s=float(np.max(abs(
                np.sum((chart_velocity-body_speed)*center_panels.normals, axis=1)))),
            implicit_vs_resolved_body_pressure_max_pa=float(np.max(abs(body_pressure-resolved_body_pressure))),
            using_rigid_velocity_for_chart_pressure_max_change_pa=float(WATER_DENSITY*np.max(abs(
                np.sum((chart_velocity-body_speed)*surface_velocity, axis=1)))),
            omitting_operator_derivative_body_pressure_max_change_pa=float(WATER_DENSITY*np.max(abs(
                incomplete_body_response-body_response))),
            bulk_to_surface_dynamic_pressure_max_error_pa=np.max(abs(dynamic_pressure_difference), axis=1).tolist(),
            body_pressure_force_and_moment=(weighted.T@body_pressure).tolist(),
            current_domain_minus_fd_pressure_max_pa=float(np.max(abs(direct_pressure['pressure']-body_pressure))),
            current_domain_minus_fd_pressure_functionals=(direct_pressure['pressure_functionals']
                -weighted.T@body_pressure).tolist(),
            actual_graph_trace_diagnostics=dict(flow_defect_rate=defect_rate.tolist(),
                acceleration_defect=acceleration_trace_defect.tolist(),
                fd_minus_current_acceleration=(graph_fd_chi-actual_chi).tolist()),
            off_grid_trace_diagnostics=dict(
                flow_defect_rate=((plus['off_defect']-minus['off_defect'])/(2*step)).tolist(),
                acceleration_defect=off_acceleration_defect.tolist(),
                fd_minus_current_acceleration=((plus['off_actual']-minus['off_actual'])/(2*step)-off_chi).tolist()),
            implicit_minus_resolved_body_load=(weighted.T@(body_pressure-resolved_body_pressure)).tolist(),
            body_pressure_pa=body_pressure.tolist(),
            bulk_pressure_pa=bernoulli.tolist(), fixed_phi_t_m2_per_s2=implicit_fixed.tolist())
        report['rows'].append(row)
        if step == args.step/4:
            actions = [record['frozen_action'] for record in (plus, minus)]
            action_density_dot, action_solve = center.solve_density(
                (actions[0]['residual']-actions[1]['residual'])/(2*step))
            action_response_spectrum = -center.spectrum(center.trace@cp.asarray(action_density_dot))
            action_response = center.body_potential(action_density_dot, action_response_spectrum)
            action_body_rate = (actions[0]['potential']-actions[1]['potential'])/(2*step)+action_response
            action_phi_t = action_body_rate-chart_convection
            action_pressure = -WATER_DENSITY*(action_phi_t+.5*np.sum(surface_velocity**2, axis=1)
                                              +GRAVITY*center_panels.centroids[:, 2])
            report['frozen_field_pressure'] = dict(step_seconds=step,
                offsets=[a['diagnostics'] for a in actions], density_derivative_solve=action_solve,
                implicit_full_context_pressure_max_difference_pa=float(np.max(abs(action_pressure-body_pressure))),
                implicit_full_context_six_functional_difference=(weighted.T@(action_pressure-body_pressure)).tolist(),
                scope='Two frozen-density field actions plus one CURRENT derivative solve; full contexts remain independent verification references. No sixDOF evolution or end-to-end timing claim.')
            # Genuine zero-world-acceleration stencils own the coupled known RHS.
            # Subtracting W*a from an accelerating stencil leaves O(h^2) geometry terms.
            modes = weighted/center_panels.areas[:, None]
            prescribed_acceleration = .3*np.r_[initial.linear_velocity, omega]
            constant_actions, constant_panels = [], []
            for offset_time in [step, -step]:
                pose = constant_twist_pose(initial, offset_time)
                panels, mask, owners, pattern = wet_snapshot(inputs, pose, plane_origin, frame,
                    grid_xy, args.length, args.body_refinement)
                if (len(panels) != len(center_panels) or not np.array_equal(owners, source)
                        or not np.array_equal(pattern, signs) or not np.array_equal(mask, inside)):
                    raise ValueError('zero-acceleration stencil changed clipping/exterior ownership')
                plane_bound = float(np.max(panels.triangles[..., 2]+level
                    +(panels.triangles[..., :2]-cg[:2])@slope))/2
                required_image = image_height(panels.triangles, eta, [args.length]*2, origin)[0]
                actual_depth = fourier_points(eta, panels.centroids[:, :2], args.length, origin)[0]-panels.centroids[:, 2]
                if (required_image > image or plane_bound >= image
                        or panels.triangles[..., 2].max() > height or np.min(actual_depth) <= 0):
                    raise ValueError('zero-acceleration stencil failed image/continuation/physical-depth guard')
                trace = extension.extend(.1*np.cos(grid_xy[extension.exterior]@k-wave_frequency*offset_time))
                action = frozen_field(panels, eta, args.length, origin, height, image,
                                      args.order, trace, center_solution[0])
                angular = rotation_matrix(pose.orientation)@pose.angular_velocity
                speed = pose.linear_velocity+np.cross(angular, panels.centroids-pose.position)
                action['residual'] = np.sum(speed*panels.normals, axis=1)-action['normal']
                constant_actions.append(action)
                constant_panels.append(panels)
            normal_rhs = (constant_actions[0]['residual']-constant_actions[1]['residual'])/(2*step)
            constant_chart_velocity = (constant_panels[0].centroids-constant_panels[1].centroids)/(2*step)
            frozen_phi_t = (constant_actions[0]['potential']-constant_actions[1]['potential'])/(2*step) \
                -np.sum(constant_chart_velocity*surface_velocity, axis=1)
            coupled = coupled_pressure(center, inputs.body, initial, normal_rhs,
                                       frozen_phi_t, surface_velocity)
            known_density, _ = center.solve_density(normal_rhs)
            known_spectrum = -center.spectrum(center.trace@cp.asarray(known_density))
            known_phi_t = frozen_phi_t+center.body_potential(known_density, known_spectrum)
            known_pressure = -WATER_DENSITY*(known_phi_t+.5*np.sum(surface_velocity**2, axis=1)
                                            +GRAVITY*center_panels.centroids[:, 2])
            responses = [center.solve(np.zeros(args.grid**2), modes[:, axis]) for axis in range(6)]
            pressure_modes = np.column_stack([solution.potential for solution in responses])
            fluid_mass = WATER_DENSITY*weighted.T@pressure_modes
            mass, external = body_dynamics(inputs.body, initial)
            modal_acceleration = np.linalg.solve(mass+fluid_mass, weighted.T@known_pressure+external)
            modal_pressure = known_pressure-WATER_DENSITY*pressure_modes@modal_acceleration
            recovered_prescribed = known_pressure-WATER_DENSITY*pressure_modes@prescribed_acceleration
            report['coupled_current_pressure'] = dict(diagnostics=coupled['diagnostics'],
                known_rhs_stencil='Genuine zero-world-linear/angular-acceleration poses; explicit clipping/label/image/height/physical-depth guards',
                zero_acceleration_offsets=[a['diagnostics'] for a in constant_actions],
                acceleration_world=coupled['acceleration'].tolist(),
                explicit_modal_acceleration_world=modal_acceleration.tolist(),
                acceleration_difference=(coupled['acceleration']-modal_acceleration).tolist(),
                pressure_max_difference_pa=float(np.max(abs(coupled['pressure']-modal_pressure))),
                six_functional_difference=(weighted.T@(coupled['pressure']-modal_pressure)).tolist(),
                recovered_prescribed_pressure_max_difference_pa=float(np.max(abs(recovered_prescribed-body_pressure))),
                prescribed_recovery_scope='Zero-acceleration and accelerating centered paths differ at O(h^2); this is a stencil convergence comparison, not roundoff parity.',
                body_mass_matrix_world=mass.tolist(), fluid_mass_matrix_world=fluid_mass.tolist(),
                fluid_mass_max_asymmetry=float(np.max(abs(fluid_mass-fluid_mass.T))),
                scope='One matrix-free rank-six coupled solve vs seven independent current modal solves; no matrix symmetry repair, no evolved trajectory.')
    report.update(complete=True, total_run_seconds=perf_counter()-started,
        cupy_pool_used_bytes=cp.get_default_memory_pool().used_bytes(),
        cupy_pool_reserved_bytes=cp.get_default_memory_pool().total_bytes())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--grid', type=int, default=129)
    parser.add_argument('--length', type=float, default=64.)
    parser.add_argument('--order', type=int, choices=range(2, 9), default=4)
    parser.add_argument('--body-refinement', type=int, choices=[1, 2], default=1)
    parser.add_argument('--step', type=float, default=.04)
    parser.add_argument('--image-height', type=float,
                        help='hold a valid image height across spatial-resolution controls')
    parser.add_argument('--reference-height', type=float,
                        help='hold a valid expansion height across derivative-stencil controls')
    args = parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists()
            or args.grid < 17 or args.grid % 2 != 1 or not np.isfinite(args.length) or args.length <= 0
            or not np.isfinite(args.step) or args.step <= 0):
        parser.error('fresh output/archive, odd grid>=17 and positive finite length/step required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True)
    hashes = {}
    sources = [*Path('src/topwater').glob('*.py'), *Path('tests').glob('*.py'),
               *Path(__file__).parent.glob('*.py'), *Path(__file__).parent.glob('*.cu'),
               *Path(__file__).parent.glob('*.cuh')]
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    report = dict(complete=False, qualified=False, scope=__doc__, producer_sha256=hashes,
        snapshots=[], rows=[], limits=[
            'Prescribed small-pose stencil, fixed forced graph/Dirichlet data; not free-surface evolution or sixDOF advancement',
            'Centroid surface-pressure integration and bulk probes; no spatial/load-convergence or physical hull-force qualification',
            'Implicit derivative and offset solves share the elliptic discretization; this is consistency, not PDE accuracy',
            'Current-domain transport pressure is PDE-consistent, not an exact derivative of the finite source/HOS span; spatial controls are required',
            'Exterior labels and clipping recipes must persist; newly exposed state transfer is not implemented',
            'Image/expansion heights are fixed across the entire stencil; no representation-switch derivative',
            'Aperiodic source/periodic HOS discrepancy and earlier finite-frequency radiation defect remain open'])
    started = perf_counter()
    try:
        run(SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8')), args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        report['peak_process_rss_bytes'] = peak_rss_bytes()
        report['memory_scope'] = 'Process peak RSS; ending CUDA pool is not total or peak VRAM.'
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
