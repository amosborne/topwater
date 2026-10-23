"""Flat actual-hull gate for an elevated-image/current-graph formulation.

R uses opposite images ABOVE the fluid. H restores the ACTUAL plane trace:
phi = H(-R|z=0) + R for prescribed body velocity and zero plane potential.
Changing the artificial image height must not change the continuum answer.
This is a static precursor, not a nonlinear solver or a runtime alternative.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.spatial import ConvexHull

from topwater import SimulationInputs
from topwater.bem import assemble_dirichlet_halfspace
from topwater.constants import WATER_DENSITY, GRAVITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels, panel_influence
from tests.benchmark_rankine import peak_rss_bytes
from .boat_spatial_work import source_matrix, subdivide
from .boat_radiation_time import harmonic_maps
from .boat_graph_images import image_height
from .panel_gradient_reference import checked_influence


def _potential_block(points, panels, height):
    reflected = points.copy()
    reflected[:, 2] = 2*height-reflected[:, 2]
    return (panel_influence(points, panels).single_layer
            -panel_influence(reflected, panels).single_layer)


def elevated_potential(points, panels, height):
    """Parallel target blocks; unchanged guarded CPU integration, no worker CUDA."""
    value = np.empty((len(points), len(panels)))
    starts = range(0, len(points), 128)
    with ThreadPoolExecutor(max_workers=4) as pool:
        blocks = pool.map(lambda first: _potential_block(points[first:first+128], panels, height), starts)
        for first, block in zip(starts, blocks, strict=True):
            value[first:first+len(block)] = block
    return value


def elevated_trace(points, panels, height):
    """Unperiodized paired-source potential and vertical derivative, blocked."""
    value = np.empty((len(points), len(panels)))
    vertical = np.empty_like(value)
    for first in range(0, len(points), 128):
        last = min(first+128, len(points))
        real = checked_influence(points[first:last], panels)
        reflected = points[first:last].copy()
        reflected[:, 2] = 2*height-reflected[:, 2]
        image = checked_influence(reflected, panels)
        value[first:last] = real.single_layer-image.single_layer
        vertical[first:last] = real.gradient_single_layer[:, :, 2]+image.gradient_single_layer[:, :, 2]
    return value, vertical


def exterior_probes(segments, footprint, beam):
    """Fixed physical targets, independent of resolution, at each segment."""
    center = footprint.points[footprint.vertices].mean(axis=0)
    points = []
    for a, b in segments:
        middle = .37*a+.63*b
        edge = b[:2]-a[:2]
        normal = np.array([edge[1], -edge[0]])/np.linalg.norm(edge)
        if normal@(middle[:2]-center) < 0:
            normal *= -1
        for fraction in [.01, .03]:
            points.append(middle+np.array([*(fraction*beam*normal), 0.]))
    points = np.asarray(points)
    if not np.all(np.any(points[:, :2]@footprint.equations[:, :2].T
                         +footprint.equations[:, 2] > 0, axis=1)):
        raise ValueError('probe must lie strictly on the physical exterior plane')
    return points


def run(inputs, reference, grids):
    started = perf_counter()
    from .cuda_operator import cp
    wet = hydrostatics(inputs.body, inputs.initial_state, density=WATER_DENSITY, gravity=GRAVITY)
    triangles = subdivide(wet.wetted_triangles)
    panels = TrianglePanels(triangles, -1)
    if np.any(panels.centroids[:, 2] >= 0):
        raise ValueError('strictly submerged body centroids required')
    modes = np.column_stack((panels.normals,
        np.cross(panels.centroids-inputs.initial_state.position, panels.normals)))
    right = modes[:, [2, 3]]
    weighted_modes = panels.areas[:, None]*modes
    zero_single = assemble_dirichlet_halfspace(panels).single_layer
    zero_matrix = source_matrix(panels)
    zero_potential = zero_single@np.linalg.solve(zero_matrix, right)
    zero_inertia = WATER_DENSITY*weighted_modes.T@zero_potential
    panel_hash = sha256(triangles.tobytes()).hexdigest()
    refs = [row['reference'] for row in reference['rows']
            if row['panel_sha256'] == panel_hash and row['reference']['angular_frequency'] == 'infinite'
            and row['reference']['density'] == WATER_DENSITY
            and row['reference']['gravity'] == GRAVITY
            and row['reference']['water_depth'] == 'infinite'
            and row['reference']['radiating_modes'] == ['Heave', 'Roll']
            and row['reference']['influenced_modes'] == ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']]
    if len(refs) != 1:
        raise ValueError('one exact-panel infinite-frequency reference required')
    independent = np.asarray(refs[0]['added_mass'])
    footprint = ConvexHull(wet.waterline_segments[:, :, :2].reshape(-1, 2))
    if not np.isclose(footprint.volume, wet.waterplane_area, rtol=1e-12, atol=1e-12):
        raise ValueError('convex physical waterplane required')
    beam = float(np.ptp(inputs.body.mesh.vertices[:, 1]))
    rows = []
    for length, size in grids:
        before = perf_counter()
        spacing = length/size
        origin = inputs.initial_state.position[:2]-length/2
        if np.any(np.abs(footprint.points-inputs.initial_state.position[:2]) >= length/2):
            raise ValueError('box must strictly contain the actual waterplane')
        potential, normal, wave = harmonic_maps(panels.centroids, panels.normals, size, spacing, origin)
        potential, normal, wave = cp.asarray(potential), cp.asarray(normal), cp.asarray(wave)
        x, y = np.meshgrid(origin[0]+np.arange(size)*spacing,
                           origin[1]+np.arange(size)*spacing)
        points = np.column_stack((x.ravel(), y.ravel(), np.zeros(size*size)))
        exterior = np.any(points[:, :2]@footprint.equations[:, :2].T+footprint.equations[:, 2] > 0, axis=1)
        probes = exterior_probes(wet.waterline_segments, footprint, beam)
        cp.cuda.get_current_stream().synchronize()
        map_seconds = perf_counter()-before
        local_height, local_bound = image_height(triangles, np.zeros((size, size)),
                                                [length, length], origin)
        for plane_kind, height in [('local graph bound', local_height),
                                   ('quarter beam', .25*beam), ('half beam', .5*beam)]:
            case_started = perf_counter()
            shifted = TrianglePanels(triangles-np.array([0., 0., height]), -1)
            matrix = source_matrix(shifted)
            single = assemble_dirichlet_halfspace(shifted).single_layer
            before = perf_counter()
            trace, vertical = elevated_trace(points, panels, height)
            source_trace_seconds = perf_counter()-before
            before = perf_counter()
            spectrum = cp.fft.fft2(cp.asarray(trace).reshape(size, size, len(panels)),
                                  axes=(0, 1)).reshape(-1, len(panels))
            p_lift, n_lift = potential@spectrum, normal@spectrum
            imaginary = float(cp.max(cp.abs(p_lift.imag)).item())
            imaginary_normal = float(cp.max(cp.abs(n_lift.imag)).item())
            effective_single = single-cp.asnumpy(p_lift.real)
            effective_matrix = matrix-cp.asnumpy(n_lift.real)
            cp.cuda.get_current_stream().synchronize()
            coupling_seconds = perf_counter()-before
            before = perf_counter()
            density = lu_solve(lu_factor(effective_matrix), right)
            body_potential = effective_single@density
            inertia = WATER_DENSITY*weighted_modes.T@body_potential
            solve_seconds = perf_counter()-before
            before = perf_counter()
            combined = spectrum@cp.asarray(density)
            h_vertical = cp.fft.ifft2((wave.ravel()[:, None]*combined).reshape(size, size, 2),
                                      axes=(0, 1)).real.reshape(-1, 2)
            flux = cp.asarray(vertical)@cp.asarray(density)-h_vertical
            surface_flux = spacing**2*cp.sum(flux[cp.asarray(exterior)], axis=0)
            direct, _ = elevated_trace(probes, panels, height)
            probe_field = direct@density
            for first in range(0, len(probes), 8):
                last = min(first+8, len(probes))
                p, _, _ = harmonic_maps(probes[first:last],
                    np.tile([0., 0., 1.], (last-first, 1)), size, spacing, origin)
                probe_field[first:last] -= cp.asnumpy((cp.asarray(p)@combined).real)
            cp.cuda.get_current_stream().synchronize()
            audit_seconds = perf_counter()-before
            difference = body_potential-zero_potential
            relative = np.sqrt(np.sum(panels.areas[:, None]*difference**2, axis=0)
                               /np.sum(panels.areas[:, None]*zero_potential**2, axis=0))
            row = dict(box_length_m=length, grid_size=size, spacing_m=spacing, plane_kind=plane_kind,
                image_plane_height_m=height, local_image_height_bound=local_bound,
                minimum_image_point_height_m=2*height-float(np.max(triangles[:, :, 2])),
                added_mass=inertia.tolist(), zero_plane_added_mass_error=(inertia-zero_inertia).tolist(),
                capytaine_added_mass_error=(inertia-independent).tolist(),
                body_potential_relative_l2=relative.tolist(),
                effective_single_relative_frobenius_error=float(np.linalg.norm(effective_single-zero_single)/np.linalg.norm(zero_single)),
                effective_normal_relative_frobenius_error=float(np.linalg.norm(effective_matrix-zero_matrix)/np.linalg.norm(zero_matrix)),
                offgrid_plane_trace_max=np.max(np.abs(probe_field), axis=0).tolist(),
                offgrid_plane_trace_relative_to_body_peak=(np.max(np.abs(probe_field), axis=0)
                    /np.max(np.abs(zero_potential), axis=0)).tolist(),
                exterior_surface_flux_rectangle=cp.asnumpy(surface_flux).tolist(),
                prescribed_body_flux=(panels.areas@right).tolist(),
                mean_flux_closure_rectangle=cp.asnumpy(surface_flux+cp.asarray(panels.areas@right)).tolist(),
                matrix_condition_number_2=float(np.linalg.cond(effective_matrix)),
                centroid_neumann_relative_residual=float(np.linalg.norm(effective_matrix@density-right)/np.linalg.norm(right)),
                max_imaginary_potential_lift=imaginary, max_imaginary_normal_lift=imaginary_normal,
                shared_grid_map_seconds=map_seconds, source_trace_seconds=source_trace_seconds,
                coupling_seconds=coupling_seconds, solve_seconds=solve_seconds,
                audit_seconds=audit_seconds, case_seconds=perf_counter()-case_started)
            rows.append(row)
            print(json.dumps(row), flush=True)
            # These maps are owned by this image-height case. Releasing them
            # before the next FFT avoids retaining two full source spectra.
            del spectrum, combined, flux, h_vertical, p_lift, n_lift, trace, vertical
        del potential, normal, wave
    return dict(complete=True, qualified=False, scope=__doc__,
        input_sha256=sha256(inputs.to_json().encode()).hexdigest(), panel_sha256=panel_hash,
        body_panel_count=len(panels), driven_modes=['Heave', 'Roll'],
        fixed_probe_distances_m=[.01*beam, .03*beam],
        force_rows=['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'],
        coefficient_units=[['kg', 'kg*m']]*3+[['kg*m', 'kg*m^2']]*3,
        zero_plane_added_mass=zero_inertia.tolist(), capytaine_added_mass=independent.tolist(),
        zero_plane_capytaine_error=(zero_inertia-independent).tolist(),
        image_bound='Local source projection; Fourier derivative bounds and reported floating allowance. The flat graph here is eta=0.',
        device=cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)['name'].decode(),
        cupy_version=cp.__version__, peak_process_rss_bytes=peak_rss_bytes(),
        cupy_pool_used_bytes=cp.get_default_memory_pool().used_bytes(),
        cupy_pool_reserved_bytes=cp.get_default_memory_pool().total_bytes(),
        memory_scope='Process peak RSS; ending CUDA pool is not peak or total VRAM.',
        timing_scope='Run includes imports/setup/maps/all cases; shared_grid_map_seconds counted once per grid, not once per image height. Source archiving and final JSON outside run.',
        rows=rows, total_seconds=perf_counter()-started,
        limits=['Flat graph only; no HOS truncation, moving hull, nonlinear pressure or time evolution',
                'Full auxiliary source trace is sampled under hull too; no hidden physical states and no harmonic fill of the source trace',
                'Unperiodized paired sources plus periodic Poisson lift require independent domain control',
                'Rectangle surface flux excludes hidden nodes but does not integrate cut cells; not a certified mass-conservation error',
                'Centroid Neumann residual is algebraic, not off-collocation physical boundary accuracy'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--length', type=float, nargs='+', default=[64., 64., 64., 128., 128.])
    parser.add_argument('--grid', type=int, nargs='+', default=[129, 257, 513, 257, 513])
    args = parser.parse_args()
    if (args.output.exists() or len(args.length) != len(args.grid)
            or any(not np.isfinite(x) or x <= 0 for x in args.length)
            or any(x < 17 or x % 2 != 1 for x in args.grid)):
        parser.error('new output, paired positive lengths and odd grids>=17 required')
    inputs = SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8'))
    reference = json.loads(args.reference.read_text(encoding='utf-8'))
    if reference['input_sha256'] != sha256(inputs.to_json().encode()).hexdigest():
        parser.error('reference physical inputs differ')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    sources = [*Path('src/topwater').glob('*.py'), Path(__file__),
        *[Path('tests')/name for name in ['laplace_reference.py', 'benchmark_rankine.py', 'reference_hull.py']],
        *[Path(__file__).with_name(name) for name in ['boat_spatial_work.py', 'boat_radiation_time.py',
          'panel_gradient_reference.py', 'boat_graph_images.py', 'cuda_operator.py', 'cuda_pair.cuh', 'cuda_energy.cu']]]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    shutil.copyfile(args.reference, archive/'reference.json')
    result = run(inputs, reference, list(zip(args.length, args.grid)))
    result['producer_sha256'] = hashes
    result['reference_sha256'] = sha256(args.reference.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
