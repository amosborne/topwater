"""Nonflat HOS continuation falsifier on actual reference-boat targets.

The true graph equals the observed sea's local tangent plane over the entire
hull, and tapers smoothly outside it. The hull is clipped to that true plane.
Known paired image layers lie outside both true and represented graphs.
This tests a necessary continuation component, not a body BVP or trajectory.
Pressure metrics are same-point centroid functionals, not converged integrals.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.spatial.transform import Rotation

from topwater import SimulationInputs, InitialState
from topwater.constants import WATER_DENSITY, GRAVITY
from topwater.frames import rotation_matrix
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from .boat_graph_images import image_height
from .boat_spatial_work import subdivide
from .boat_radiation_time import harmonic_maps
from .hos_graph import graph_field
from .incident_baseline_benchmark import coefficients, fused_cpu
from .panel_gradient_reference import checked_influence


def tilted_geometry(inputs):
    cg = inputs.initial_state.position
    incident = fused_cpu(np.array([[*cg[:2], 0.]]), coefficients(inputs.waves), 0.)[0]
    level, slope = incident[0], incident[1:3]
    normal = np.r_[-slope, 1.]
    normal /= np.linalg.norm(normal)
    tangent = np.array([1., 0., slope[0]])
    tangent /= np.linalg.norm(tangent)
    transform = np.vstack((tangent, np.cross(normal, tangent), normal))
    origin = np.r_[cg[:2], level]
    orientation = (Rotation.from_matrix(transform)
        *Rotation.from_quat(inputs.initial_state.orientation, scalar_first=True)).as_quat(scalar_first=True)
    state = InitialState(position=transform@(cg-origin), orientation=orientation)
    wet = hydrostatics(inputs.body, state, density=WATER_DENSITY, gravity=GRAVITY)
    triangles = wet.wetted_triangles@transform+origin
    return subdivide(triangles), float(level), slope


def window(coordinate, inner, outer):
    u = (np.abs(coordinate)-inner)/(outer-inner)
    result = np.ones_like(u)
    result[u >= 1] = 0
    middle = (u > 0) & (u < 1)
    a, b = np.exp(-1/u[middle]), np.exp(-1/(1-u[middle]))
    result[middle] = b/(a+b)
    return result


def image_field(points, panels, density, heights):
    """Known harmonic field from two opposite EXTERIOR image layers."""
    value = np.empty((len(points), density.shape[1]))
    gradient = np.empty((len(points), density.shape[1], 3))
    for first in range(0, len(points), 128):
        last = min(first+128, len(points))
        one, two = points[first:last].copy(), points[first:last].copy()
        one[:, 2], two[:, 2] = 2*heights[0]-one[:, 2], 2*heights[1]-two[:, 2]
        a, b = checked_influence(one, panels), checked_influence(two, panels)
        value[first:last] = (a.single_layer-b.single_layer)@density
        gradient[first:last] = np.einsum('tnj,nm->tmj',
            (a.gradient_single_layer-b.gradient_single_layer)*np.array([1., 1., -1.]), density)
    return value, gradient


def run(inputs, grid_cases, orders, prescribed_image_height=None):
    started = perf_counter()
    triangles, level, slope = tilted_geometry(inputs)
    panels = TrianglePanels(triangles, -1)
    cg = inputs.initial_state.position
    beam = float(np.ptp(inputs.body.mesh.vertices[:, 1]))
    world = (inputs.body.mesh.vertices-inputs.body.cg)@rotation_matrix(inputs.initial_state.orientation).T+cg
    inner = np.max(np.abs(world[:, :2]-cg[:2]), axis=0)+beam
    outer = inner+2*beam
    expansion_height = float(triangles[..., 2].max())
    modes = np.column_stack((panels.normals, np.cross(panels.centroids-cg, panels.normals)))
    density = modes[:, [2, 3]]
    load = -WATER_DENSITY*modes.T*panels.areas
    grids = []
    for length, size in grid_cases:
        if np.any(outer >= length/2):
            raise ValueError('smooth graph support must fit strictly inside the box')
        origin = cg[:2]-length/2
        x, y = np.meshgrid(origin[0]+np.arange(size)*length/size,
                           origin[1]+np.arange(size)*length/size)
        xy = np.stack((x-cg[0], y-cg[1]), axis=-1)
        eta = (level+xy@slope)*window(x-cg[0], inner[0], outer[0])*window(y-cg[1], inner[1], outer[1])
        height, bound = image_height(triangles, eta, [length, length], origin)
        grids.append((length, size, origin, x, y, eta, height, bound))
    true_plane_max = float(np.max(level+(triangles[..., :2]-cg[:2])@slope))
    safe_height = max([item[6] for item in grids]
                      +[.5*(expansion_height+true_plane_max), expansion_height])
    if prescribed_image_height is not None:
        if not np.isfinite(prescribed_image_height) or prescribed_image_height < safe_height:
            raise ValueError('prescribed manufactured image height must meet all graph bounds')
        safe_height = prescribed_image_height
    rows = []
    for length, size, origin, x, y, eta, _, bound in grids:
        padding = 2+1/(4*size)  # ceil gives odd2N+1: no Nyquist ambiguity at hull targets.
        padded = int(np.ceil(padding*size))
        targets = panels.centroids.copy()
        targets[:, 2] -= expansion_height
        p, n, _ = harmonic_maps(targets, panels.normals, padded, length/padded, origin)
        plane_points = panels.centroids.copy()
        plane_points[:, 2] = 0
        eta_map, _, _ = harmonic_maps(plane_points, np.zeros_like(plane_points), size, length/size, origin)
        represented_eta = (eta_map@np.fft.fft2(eta).ravel()).real
        if np.any(represented_eta <= panels.centroids[:, 2]):
            raise ValueError('hull continuation targets must lie below represented auxiliary graph')
        graph_points = np.column_stack((x.ravel(), y.ravel(), eta.ravel()))
        near = (np.abs(x-cg[0]) <= inner[0]) & (np.abs(y-cg[1]) <= inner[1])
        kx, ky = np.meshgrid(2*np.pi*np.fft.fftfreq(size, d=length/size),
                             2*np.pi*np.fft.fftfreq(size, d=length/size))
        eta_spectrum = np.fft.fft2(eta)
        slopes = np.stack((np.fft.ifft2(1j*kx*eta_spectrum).real,
                           np.fft.ifft2(1j*ky*eta_spectrum).real), axis=-1).reshape(-1, 2)
        flat_points = graph_points.copy()
        flat_points[:, 2] = expansion_height
        for clearance in [0., .1*beam]:
            heights = safe_height+clearance, safe_height+clearance+.25*beam
            true_plane_at_sources = level+(triangles[..., :2]-cg[:2])@slope
            image_gaps = [float(np.min(2*h-triangles[..., 2]-true_plane_at_sources)) for h in heights]
            if min(image_gaps) <= 0:
                raise ValueError('manufactured image layers must lie strictly above the true graph')
            exact, gradient = image_field(panels.centroids, panels, density, heights)
            normal_exact = np.einsum('tmj,tj->tm', gradient, panels.normals)
            psi, graph_gradient = image_field(graph_points, panels, density, heights)
            flat_psi, _ = image_field(flat_points, panels, density, heights)
            flat = np.column_stack([(p@np.fft.fft2(graph_field(np.zeros_like(eta),
                flat_psi[:, col].reshape(size, size), [length, length], order=1,
                padding=padding).flat_terms[0]).ravel()).real for col in range(2)])
            body_norm = lambda a: np.sqrt(np.sum(panels.areas[:, None]*a*a, axis=0))
            for order in orders:
                before = perf_counter()
                answer, normal_answer, graph_error, near_error, periodic_flux, reference_flux = [], [], [], [], [], []
                for col in range(2):
                    result = graph_field(eta-expansion_height, psi[:, col].reshape(size, size),
                                         [length, length], order=order, padding=padding)
                    spectrum = np.fft.fft2(sum(result.flat_terms)).ravel()
                    answer.append((p@spectrum).real)
                    normal_answer.append((n@spectrum).real)
                    difference = result.velocity.reshape(-1, 3)-graph_gradient[:, col]
                    graph_error.append(float(np.linalg.norm(difference)/np.linalg.norm(graph_gradient[:, col])))
                    near_error.append(float(np.linalg.norm(difference[near.ravel()])
                        /np.linalg.norm(graph_gradient[near.ravel(), col])))
                    periodic_flux.append(float(result.normal_flux.sum()*(length/size)**2))
                    q = graph_gradient[:, col, 2]-np.sum(slopes*graph_gradient[:, col, :2], axis=1)
                    reference_flux.append(float(q.sum()*(length/size)**2))
                answer, normal_answer = np.column_stack(answer), np.column_stack(normal_answer)
                row = dict(box_length_m=length, grid_size=size, padded_size=padded, hos_order=order,
                    image_clearance_from_common_bound_m=clearance, image_heights_m=list(heights),
                    minimum_image_gap_above_true_plane_m=image_gaps,
                    local_image_bound=bound,
                    body_potential_relative_l2=(body_norm(answer-exact)/body_norm(exact)).tolist(),
                    body_normal_relative_l2=(body_norm(normal_answer-normal_exact)/body_norm(normal_exact)).tolist(),
                    flat_control_body_potential_relative_l2=(body_norm(flat-exact)/body_norm(exact)).tolist(),
                    graph_gradient_relative_sample_l2=graph_error,
                    near_hull_graph_gradient_relative_sample_l2=near_error,
                    auxiliary_periodic_flux_integral=periodic_flux,
                    aperiodic_reference_graph_flux_integral=reference_flux,
                    centroid_pressure_functional_error=(load@(answer-exact)).tolist(),
                    exact_centroid_pressure_functionals=(load@exact).tolist(),
                    represented_minus_true_graph_at_hull_centroids_max_m=float(np.max(np.abs(
                        represented_eta-(level+(panels.centroids[:, :2]-cg[:2])@slope)))),
                    minimum_target_depth_below_represented_graph_m=float(np.min(represented_eta-panels.centroids[:, 2])),
                    hos_and_continuation_seconds=perf_counter()-before)
                rows.append(row)
                print(json.dumps(row), flush=True)
    return dict(complete=True, qualified=False, scope=__doc__,
        input_sha256=sha256(inputs.to_json().encode()).hexdigest(),
        body_panel_count=len(panels), panel_sha256=sha256(triangles.tobytes()).hexdigest(),
        tangent_plane_height_at_cg_m=level, tangent_plane_slopes=slope.tolist(),
        expansion_plane_height_m=expansion_height, common_safe_image_height_m=safe_height,
        manufactured_field_normalization='Acceleration potentials in m^2/s^2: density patterns multiply unit1m/s^2 translation or1rad/s^2 rotation; gradients are acceleration fields. No physical time derivative was solved.',
        source_field_columns=['heave-normal pattern', 'world-x-rotation-normal pattern'],
        centroid_functional_rows=['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'],
        centroid_functional_units=['N']*3+['Nm']*3,
        rows=rows, total_seconds=perf_counter()-started, peak_process_rss_bytes=peak_rss_bytes(),
        limits=['CPU component reference, not GPU or complete-simulation throughput',
                'Known external images test HOS continuation; no unknown body density or rigid-body motion is solved',
                'Centroid functionals compare identical sample quadrature, not independently converged surface-pressure integrals',
                'Finite periodic continuation of aperiodic image fields retains domain bias; flat controls expose part of it',
                'Centroid depth guards do not certify whole wet panels below the represented Fourier graph; true wet cut is to the tangent plane',
                'Image heights fixed across grids; HOS expansion order and spatial resolution varied independently'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--length', type=float, nargs='+', default=[64., 64., 64.])
    parser.add_argument('--grid', type=int, nargs='+', default=[65, 129, 257])
    parser.add_argument('--order', type=int, nargs='+', choices=range(2, 9), default=[2, 4, 6, 8])
    parser.add_argument('--image-height', type=float)
    args = parser.parse_args()
    if (args.output.exists() or len(args.length) != len(args.grid)
            or any(not np.isfinite(x) or x <= 0 for x in args.length)
            or any(x < 17 or x % 2 != 1 for x in args.grid)):
        parser.error('new output, paired positive lengths and odd grids>=17 required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    sources = [*Path('src/topwater').glob('*.py'), Path(__file__),
        *[Path('tests')/name for name in ['laplace_reference.py', 'benchmark_rankine.py', 'reference_hull.py']],
        *[Path(__file__).with_name(name) for name in ['boat_graph_images.py', 'boat_spatial_work.py',
          'boat_radiation_time.py', 'hos_graph.py', 'incident_baseline_benchmark.py', 'panel_gradient_reference.py']]]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    result = run(SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8')),
                 list(zip(args.length, args.grid)), args.order, args.image_height)
    result['producer_sha256'] = hashes
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
