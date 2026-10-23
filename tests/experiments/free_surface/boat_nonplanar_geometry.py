"""Actual-sea curved graph versus local P1 wet-hull clipping; no fluid solve.

The Fourier graph is a windowed observed-sea snapshot. Local P1 parents sample
that SAME interpolant. Geometric six-mode closure and symmetric mask difference
expose disagreement without fitting a fluid result or choosing a runtime mesh.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.spatial import ConvexHull

from topwater import InitialState, SimulationInputs, WaveComponents, sea_states
from topwater.constants import GRAVITY
from topwater.frames import rotation_matrix
from topwater.hydrostatics import FloodingError
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from .boat_hos_continuation import window
from .boat_world_boundary import construct
from .boat_parent_clip_proof import coplanar_groups
from .incident_baseline_benchmark import coefficients, fused_cpu
from .graph_geometry_reference import tensor_fourier, fourier_points, structured_P1


def local_geometry_inputs(inputs):
    """One static XY translation; physical axes and wave realization are unchanged."""
    pose, waves = inputs.initial_state, inputs.waves
    origin = np.r_[pose.position[:2], 0.]
    local_pose = InitialState(pose.position-origin, pose.orientation,
                              pose.linear_velocity, pose.angular_velocity)
    phase_shift = waves.angular_frequencies**2/GRAVITY*(
        origin[0]*np.cos(waves.directions)+origin[1]*np.sin(waves.directions))
    local_waves = WaveComponents(waves.amplitudes, waves.angular_frequencies,
        waves.directions, waves.phases+phase_shift, seed=waves.seed,
        provenance={**waves.provenance, 'local_frame_origin_world_xy_m': str(origin[:2].tolist())})
    return origin, local_pose, local_waves


def local_axis(radius, count):
    """Evaluate integer-indexed fractions before scaling; preserve zero/half steps."""
    return radius*((2*np.arange(count+1)-count)/count)


def triangulation_failure(error):
    """Inspect the actual failed post-split fan, without bypassing validation."""
    result = {}
    trace = error.__traceback__
    while trace is not None:
        frame = trace.tb_frame
        if frame.f_code.co_name == 'construct' and 'world' in frame.f_locals:
            world, roof = frame.f_locals['world'], frame.f_locals['roof']
            center = (min(world[:, 1])+max(world[:, 1]))/2
            nearest = min(set(roof[:, 1]), key=lambda y: abs(y-center))
            result['nearest_parent_y_minus_exact_hull_bbox_midpoint_m'] = float(nearest-center)
        if frame.f_code.co_name == 'triangulate' and 'inverse' in frame.f_locals:
            values = frame.f_locals
            points, faces = values['points'], values['faces']
            samples = []
            for index, face in enumerate(faces):
                xyz = np.array([points[i] for i in face], dtype=float)
                ratios = []
                for anchor in range(3):
                    u, v = np.roll(xyz, -anchor, axis=0)[1:]-xyz[anchor]
                    denominator = np.linalg.norm(u)*np.linalg.norm(v)
                    ratios.append(float(np.linalg.norm(np.cross(u, v))/denominator) if denominator else 0.)
                if ratios[0] > 64*np.finfo(float).eps:
                    continue
                from .boat_parent_clip_proof import sub, cross, dot
                u, v = sub(points[face[1]], points[face[0]]), sub(points[face[2]], points[face[0]])
                normal = cross(u, v)
                samples.append(dict(face_index=index, kind=values['kinds'][index],
                    source=int(values['source'][index]), exact_twice_area_squared=float(dot(normal, normal)),
                    float_cross_relative_at_cyclic_anchors=ratios,
                    coordinates=xyz.tolist()))
            result['failed_faces'] = samples
        trace = trace.tb_next
    return result


def rigid_surface_flux(triangles, cg):
    """Exact six rigid-motion fluxes for consistently oriented flat triangles."""
    panels = TrianglePanels(triangles, 1)
    return np.sum(panels.areas[:, None]*np.column_stack((panels.normals,
        np.cross(panels.centroids-cg, panels.normals))), axis=0)


def hidden_mask(height, x, y, equations):
    """Classify graph points against the full convex SOLID, never its wet hull."""
    maximum = np.full(height.shape, -np.inf)
    for a, b, c, d in equations:
        maximum = np.maximum(maximum, a*x[None, :]+b*y[:, None]+c*height+d)
    tolerance = 64*np.finfo(float).eps*max(1., float(np.max(abs(height))),
                                        float(np.max(abs(x))), float(np.max(abs(y))))
    return maximum <= 0, abs(maximum) <= tolerance


def graph_flux(height, slope, x, y, mask, measure, cg):
    n = np.concatenate((-slope, np.ones((*height.shape, 1))), axis=-1)
    rx, ry, rz = x[None, :]-cg[0], y[:, None]-cg[1], height-cg[2]
    modes = np.stack((n[..., 0], n[..., 1], n[..., 2],
        ry*n[..., 2]-rz*n[..., 1], rz*n[..., 0]-rx*n[..., 2],
        rx*n[..., 1]-ry*n[..., 0]), axis=-1)
    return np.sum(modes*mask[..., None], axis=(0, 1))*measure


def fourier_dry_certificate(inputs, world, eta, length, origin):
    """Sufficient Lipschitz dry-rim/cap bounds, with explicit roundoff allowance."""
    size = eta.shape[0]
    k = 2*np.pi*np.fft.fftfreq(size, d=length/size)
    spectrum = np.fft.fft2(eta)/eta.size
    derivative = np.array([np.sum(abs(spectrum*k[None, :])),
                           np.sum(abs(spectrum*k[:, None]))])
    allowance = 64*np.finfo(float).eps*eta.size*max(1., float(np.sum(abs(spectrum))), float(abs(world).max()))
    rim = inputs.body.mesh.opening_rim
    lower, sampled = [], []
    for i, j in zip(rim, np.roll(rim, -1)):
        a, b = world[i], world[j]
        points = a+((np.arange(64)+.5)/64)[:, None]*(b-a)
        height = fourier_points(eta, points[:, :2], length, origin)[0]
        gaps = points[:, 2]-height
        bound = .5/64*(abs(b[2]-a[2])+derivative@abs(b[:2]-a[:2]))+allowance
        sampled.append(float(gaps.min()))
        lower.append(float(gaps.min()-bound))
    if sampled and min(sampled) <= 0:
        raise FloodingError('actual Fourier graph reaches a sampled physical rim point')
    if lower and min(lower) <= 0:
        raise ValueError('cannot certify dry rim against Fourier graph; not proven flooding')
    cap_bound = None
    if len(inputs.body.mesh.artificial_faces):
        cap = world[np.unique(inputs.body.mesh.faces[inputs.body.mesh.artificial_faces])]
        x = np.linspace(cap[:, 0].min(), cap[:, 0].max(), 129)
        y = np.linspace(cap[:, 1].min(), cap[:, 1].max(), 65)
        height = tensor_fourier(eta, x, y, length, origin)[0]
        error = .5*derivative@np.array([np.diff(x).max(), np.diff(y).max()])+allowance
        cap_bound = float(cap[:, 2].min()-height.max()-error)
        if cap_bound <= 0:
            raise ValueError('cannot certify artificial closure dry against Fourier graph; not proven flooding')
    return dict(rim_sample_minimum_m=min(sampled) if sampled else None,
        rim_clearance_lower_bound_m=min(lower) if lower else None,
        artificial_closure_clearance_lower_bound_m=cap_bound,
        derivative_absolute_sum_bounds=derivative.tolist(), rounding_allowance_m=allowance,
        scope='Sufficient analytic Lipschitz bounds with engineering floating allowance; not interval arithmetic')


def run(inputs, args, report):
    started = perf_counter()
    body, pose = inputs.body, inputs.initial_state
    cg = pose.position
    beam = float(np.ptp(body.mesh.vertices[:, 1]))
    coplanar_groups(body.mesh.vertices, body.mesh.faces, body.mesh.artificial_faces)
    sea_states.validate_envelope(inputs.waves, beam_m=beam)
    relative = (body.mesh.vertices-body.cg)@rotation_matrix(pose.orientation).T
    world = relative+cg
    origin = cg[:2]-args.length/2
    x, y = np.meshgrid(origin[0]+np.arange(args.grid)*args.length/args.grid,
                       origin[1]+np.arange(args.grid)*args.length/args.grid)
    inner = np.max(abs(world[:, :2]-cg[:2]), axis=0)+beam
    outer = inner+2*beam
    if np.any(outer >= args.length/2):
        raise ValueError('window support must fit strictly inside Fourier box')
    incident = fused_cpu(np.column_stack((x.ravel(), y.ravel(), np.zeros(x.size))),
                         coefficients(inputs.waves), 0.)[:, 0].reshape(x.shape)
    eta = incident*window(x-cg[0], inner[0], outer[0])*window(y-cg[1], inner[1], outer[1])
    report.update(input_sha256=sha256(inputs.to_json().encode()).hexdigest(),
        graph_sha256=sha256(eta.tobytes()).hexdigest(), graph_grid_size=args.grid,
        graph_box_length_m=args.length, dry_fourier_certificate=fourier_dry_certificate(
            inputs, world, eta, args.length, origin))
    np.savez(args.output.with_name(args.output.stem+'-graph.npz'), eta=eta, origin=origin, length=args.length)
    radius = np.max(abs(relative[:, :2]), axis=0)+beam/2
    low, high = cg[:2]-radius, cg[:2]+radius
    geometry_origin, geometry_pose, local_waves = local_geometry_inputs(inputs)
    equations = ConvexHull(world).equations
    report['local_patch_xy_bounds_m'] = [low.tolist(), high.tolist()]
    report['local_patch_definition'] = 'Integer-indexed local axes about input CG; clipping stays in that XY-translated frame'
    # Independently refined midpoint quadrature. No area estimate from FFT mask counts.
    quadratures = []
    for count in args.quadrature:
        qx = low[0]+(np.arange(count)+.5)*(high[0]-low[0])/count
        qy = low[1]+(np.arange(count//2)+.5)*(high[1]-low[1])/(count//2)
        height, slope = tensor_fourier(eta, qx, qy, args.length, origin)
        inside, ambiguous = hidden_mask(height, qx, qy, equations)
        area = (high[0]-low[0])*(high[1]-low[1])/(count*(count//2))
        flux = graph_flux(height, slope, qx, qy, inside, area, cg)
        reference = dict(midpoint_shape_yx=list(height.shape), cell_area_m2=area,
            hidden_projected_area_m2=float(inside.sum()*area), six_mode_flux=flux.tolist(),
            ambiguous_classification_area_m2=float(ambiguous.sum()*area))
        report['fourier_quadrature'].append(reference)
        quadratures.append((qx, qy, height, slope, inside, area, flux))
    for refinement in args.roof_refinement:
        nx, ny = 4*refinement, 2*refinement
        local_x, local_y = local_axis(radius[0], nx), local_axis(radius[1], ny)
        nodesx, nodesy = cg[0]+local_x, cg[1]+local_y
        heights = tensor_fourier(eta, nodesx, nodesy, args.length, origin)[0]
        xx, yy = np.meshgrid(nodesx, nodesy)
        xy = np.column_stack((xx.ravel(), yy.ravel()))
        local_xx, local_yy = np.meshgrid(local_x, local_y)
        local_xy = np.column_stack((local_xx.ravel(), local_yy.ravel()))
        parents = np.array([tri for j in range(ny) for i in range(nx)
            for a in [j*(nx+1)+i] for tri in [(a, a+1, a+nx+2), (a, a+nx+2, a+nx+1)]])
        geometry_path = args.output.with_name(args.output.stem+f'-roof{refinement}-geometry.json')
        try:
            _, arrays, geometry = construct(geometry_path, body=body, pose=geometry_pose, waves=local_waves,
                parent_xy=local_xy, parents=parents, parent_eta=heights.ravel(), parent_psi=np.zeros(len(xy)))
        except ValueError as error:
            report['triangulation_failure'] = triangulation_failure(error)
            raise
        faces = arrays['faces'][arrays['body_faces']]
        arrays['geometry_frame_origin_world_m'] = geometry_origin.copy()
        np.savez(geometry_path.with_suffix('.npz'), **arrays)
        points = arrays['points']+geometry_origin
        triangles = points[faces]  # ALREADY fluid-outward.
        TrianglePanels(triangles, 1)  # The world view must also pass the downstream kernel guard.
        body_flux = rigid_surface_flux(arrays['points'][faces], geometry_pose.position)
        all_roof = np.column_stack((local_xy, heights.ravel()))[parents]
        hidden_p1_flux = rigid_surface_flux(all_roof, geometry_pose.position)-rigid_surface_flux(
            arrays['points'][arrays['faces'][arrays['free_faces']]], geometry_pose.position)
        contact = np.intersect1d(np.unique(faces), arrays['free_nodes'])
        edges = np.sort(np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]])), axis=1)
        edge, count = np.unique(edges, axis=0, return_counts=True)
        boundary = edge[count == 1]
        if not np.isin(boundary, contact).all():
            raise ValueError('wet-body boundary contains an edge outside the P1 waterline')
        contact_points = np.vstack((points[contact],
            .5*(points[boundary[:, 0]]+points[boundary[:, 1]])))
        true_contact = fourier_points(eta, contact_points[:, :2], args.length, origin)[0]
        vertices = np.unique(triangles.reshape(-1, 3), axis=0)
        true_vertices = fourier_points(eta, vertices[:, :2], args.length, origin)[0]
        centroid = triangles.mean(axis=1)
        true_centroids = fourier_points(eta, centroid[:, :2], args.length, origin)[0]
        row = dict(roof_refinement=refinement, parents=len(parents), body_panels=len(triangles),
            geometry_frame_origin_world_m=geometry_origin.tolist(),
            body_triangles_coordinate_frame='world',
            local_body_triangles_sha256=sha256(arrays['points'][faces].tobytes()).hexdigest(),
            geometry_report=str(geometry_path), geometry_report_sha256=sha256(geometry_path.read_bytes()).hexdigest(),
            geometry_npz_sha256=sha256(geometry_path.with_suffix('.npz').read_bytes()).hexdigest(),
            geometry_seconds=geometry['seconds'], geometry_planes=geometry['planes'],
            body_triangles_sha256=sha256(triangles.tobytes()).hexdigest(),
            physical_group_count=int(arrays['source_physical_group_count']),
            fluid_outward_body_six_mode_flux=body_flux.tolist(),
            p1_hidden_six_mode_flux=hidden_p1_flux.tolist(),
            p1_self_closure_absolute_error=(body_flux-hidden_p1_flux).tolist(),
            contact_fourier_minus_p1_height_range_m=[float((true_contact-contact_points[:, 2]).min()),
                                                   float((true_contact-contact_points[:, 2]).max())],
            maximum_wet_vertex_above_fourier_graph_m=float(np.max(vertices[:, 2]-true_vertices)),
            minimum_body_centroid_depth_below_fourier_graph_m=float(np.min(true_centroids-centroid[:, 2])),
            quadratures=[])
        for qx, qy, height, slope, inside, measure, flux in quadratures:
            p1_height, p1_slope = structured_P1(local_x, local_y, heights, qx-cg[0], qy-cg[1])
            p1_inside, ambiguous = hidden_mask(p1_height, qx, qy, equations)
            union = inside | p1_inside
            if not np.any(union):
                raise ValueError('quadrature missed the hidden region; refine before using its metrics')
            dh, ds = p1_height-height, np.linalg.norm(p1_slope-slope, axis=-1)
            row['quadratures'].append(dict(midpoint_shape_yx=list(height.shape),
                body_minus_fourier_six_mode_flux=(body_flux-flux).tolist(),
                p1_hidden_projected_area_m2=float(p1_inside.sum()*measure),
                mask_symmetric_difference_area_m2=float(np.count_nonzero(p1_inside != inside)*measure),
                ambiguous_p1_classification_area_m2=float(ambiguous.sum()*measure),
                hidden_union_height_rms_m=float(np.sqrt(np.mean(dh[union]**2))),
                hidden_union_height_max_m=float(np.max(abs(dh[union]))),
                hidden_union_slope_rms=float(np.sqrt(np.mean(ds[union]**2))),
                hidden_union_slope_max=float(ds[union].max())))
        report['rows'].append(row)
        print(json.dumps({k:v for k,v in row.items() if k != 'quadratures'}), flush=True)
    report.update(complete=True, total_run_seconds=perf_counter()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--grid', type=int, default=129)
    parser.add_argument('--length', type=float, default=64.)
    parser.add_argument('--roof-refinement', type=int, nargs='+', default=[1, 2])
    parser.add_argument('--quadrature', type=int, nargs='+', default=[256, 512, 1024])
    args = parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists()
            or args.grid < 17 or args.grid % 2 != 1 or not np.isfinite(args.length) or args.length <= 0
            or any(level not in (1, 2) for level in args.roof_refinement)
            or any(n < 32 or n % 2 for n in args.quadrature)):
        parser.error('fresh output/archive, odd graph>=17, positive length, roof1/2 and even quadrature>=32 required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True)
    sources = [*Path('src/topwater').glob('*.py'), *Path('tests').glob('*.py'),
               *Path(__file__).parent.glob('*.py'), *Path(__file__).parent.glob('*.cu'),
               *Path(__file__).parent.glob('*.cuh')]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    report = dict(complete=False, qualified=False, scope=__doc__, producer_sha256=hashes,
        rows=[], fourier_quadrature=[], limits=[
            'Geometry only: windowed observed-sea snapshot, not a prepared evolving fluid or trajectory',
            'P1 and Fourier surfaces differ between nodes; wet-hull panels approximate a different waterline',
            'Midpoint mask quadrature is independently refined, not an interval-certified area integral',
            'Contact and slope diagnostics are samples; cap/rim Lipschitz certificates carry explicit floating allowance',
            'Reference clipping constructs boundary bookkeeping with existing plane/time/bed limits, not a volume mesh',
            'Only physical wet body triangles can feed later BVP; no pressure on artificial caps'])
    started = perf_counter()
    try:
        run(SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8')), args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        report['peak_process_rss_bytes'] = peak_rss_bytes()
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
