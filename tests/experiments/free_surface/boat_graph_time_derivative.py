"""Moving-hull bulk pressure/ALE derivative gate, not a trajectory or hull load.

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
from .boat_graph_images import image_height
from .boat_hos_continuation import tilted_geometry, tangent_graph
from .boat_spatial_work import subdivide
from .cuda_graph import ExteriorTrace, body_maps, cp
from .panel_gradient_reference import checked_influence


def prescribed_pose(initial, time):
    """Fixed world rotation axis; all rates agree exactly with the pose curve."""
    initial_rotation = Rotation.from_quat(initial.orientation, scalar_first=True)
    omega = initial_rotation.apply(initial.angular_velocity)
    path_time = time+.15*time*time
    orientation = Rotation.from_rotvec(omega*path_time)*initial_rotation
    return InitialState(initial.position+initial.linear_velocity*path_time,
        orientation.as_quat(scalar_first=True), initial.linear_velocity*(1+.3*time),
        orientation.inv().apply(omega*(1+.3*time)))


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
        kx, ky = cp.meshgrid(2*np.pi*cp.fft.fftfreq(size, d=length/size),
                             2*np.pi*cp.fft.fftfreq(size, d=length/size))
        self.derivatives = (1j*kx.ravel(), 1j*ky.ravel(), cp.hypot(kx, ky).ravel())
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
        return np.column_stack([cp.asnumpy((self.potential@(spectrum*d)).real)
                                for d in self.derivatives])+np.einsum('tni,n->ti', self.gradient, density)


def run(inputs, args, report):
    started = perf_counter()
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
    image = max(image_height(p.triangles, eta, [args.length]*2, origin)[0]
                for _, p, *_ in snapshots.values())
    true_image_bound = max(float(np.max(p.triangles[..., 2]+level
        +(p.triangles[..., :2]-cg[:2])@slope))/2 for _, p, *_ in snapshots.values())
    image = max(image, float(np.nextafter(true_image_bound, np.inf)))
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
    omega = rotation_matrix(initial.orientation)@initial.angular_velocity
    probe_speed = initial.linear_velocity+np.cross(omega, probes-cg)
    k = 2*np.pi*np.array([3., 1.])/args.length
    wave_frequency = float(np.sqrt(GRAVITY*np.linalg.norm(k)))
    report.update(input_sha256=sha256(inputs.to_json().encode()).hexdigest(),
        body_panels=len(center_panels), grid_size=args.grid, box_length_m=args.length,
        hos_order=args.order, body_refinement=args.body_refinement,
        image_height_m=image, reference_height_m=height, probe_count=len(probes),
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
        if time == 0:
            center, center_solution = context, (density, spectrum)
            center_field = fixed_field
            velocity = fixed_field.velocity(density, spectrum)
            record = dict(right=right, fixed=fixed, moving=fixed,
                          frozen_fixed=fixed, frozen_moving=fixed)
        else:
            frozen_spectrum = context.spectrum(trace-context.trace@cp.asarray(center_solution[0]))
            frozen_fixed = fixed_field.value(center_solution[0], frozen_spectrum)
            del fixed_field
            moving_field = BulkField(context, moving, args.length, origin)
            record = dict(right=right, fixed=fixed,
                moving=moving_field.value(density, spectrum),
                frozen_fixed=frozen_fixed,
                frozen_moving=moving_field.value(center_solution[0], frozen_spectrum),
                operator_on_center=context.apply(center_solution[0]))
            del moving_field
        records[time] = record
        report['snapshots'].append(dict(time_seconds=time, **diagnostics,
            context_seconds=context.setup_seconds,
            minimum_target_depth_below_represented_graph_m=represented_depth,
            panel_sha256=sha256(panels.triangles.tobytes()).hexdigest()))
        print(json.dumps(report['snapshots'][-1]), flush=True)
        if time != 0:
            del context
    convection = np.sum(probe_speed*velocity, axis=1)
    for step in [args.step, args.step/2, args.step/4]:
        plus, minus = records[step], records[-step]
        bdot = (plus['right']-minus['right'])/(2*step)
        adot_sigma = (plus['operator_on_center']-minus['operator_on_center'])/(2*step)
        density_dot, diagnostics = center.solve_density(bdot-adot_sigma)
        response_spectrum = -center.spectrum(center.trace@cp.asarray(density_dot))
        response = center_field.value(density_dot, response_spectrum)
        implicit_fixed = (plus['frozen_fixed']-minus['frozen_fixed'])/(2*step)+response
        implicit_material = (plus['frozen_moving']-minus['frozen_moving'])/(2*step)+response
        resolved_fixed = (plus['fixed']-minus['fixed'])/(2*step)
        resolved_material = (plus['moving']-minus['moving'])/(2*step)
        incomplete_density, _ = center.solve_density(bdot)
        incomplete_response = center_field.value(incomplete_density,
            -center.spectrum(center.trace@cp.asarray(incomplete_density)))
        bernoulli = -WATER_DENSITY*(implicit_material-convection+.5*np.sum(velocity*velocity, axis=1)
                                   +GRAVITY*probes[:, 2])
        row = dict(step_seconds=step, density_derivative_solve=diagnostics,
            implicit_vs_resolved_fixed_phi_t_max_m2_per_s2=float(np.max(abs(implicit_fixed-resolved_fixed))),
            implicit_vs_resolved_material_phi_rate_max_m2_per_s2=float(np.max(abs(implicit_material-resolved_material))),
            resolved_ale_vs_fixed_phi_t_max_m2_per_s2=float(np.max(abs(resolved_material-convection-resolved_fixed))),
            implicit_ale_vs_fixed_phi_t_max_m2_per_s2=float(np.max(abs(implicit_material-convection-implicit_fixed))),
            omitting_operator_derivative_phi_t_max_change_m2_per_s2=float(np.max(abs(incomplete_response-response))),
            omitting_ale_correction_pressure_max_change_pa=float(WATER_DENSITY*np.max(abs(convection))),
            bulk_pressure_min_max_pa=[float(bernoulli.min()), float(bernoulli.max())],
            bulk_pressure_pa=bernoulli.tolist(), fixed_phi_t_m2_per_s2=implicit_fixed.tolist())
        report['rows'].append(row)
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
            'Bulk pressures only; no surface-gradient limit, wetted pressure integration or physical hull-force qualification',
            'Implicit derivative and offset solves share the elliptic discretization; this is consistency, not PDE accuracy',
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
