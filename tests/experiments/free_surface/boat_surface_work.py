"""Compact PHYSICAL surface trace / hull reciprocal-work falsifier.

The aperiodic trace vanishes outside a disk strictly beyond the actual convex
waterplane. Its work integral therefore needs no truncated outer boundary.
The FFT lift is a separate resolution control, not an alternate runtime path.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import lu_factor, lu_solve

from topwater import SimulationInputs
from topwater.bem import assemble_dirichlet_halfspace
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from .boat_hidden_extension import poisson_lift
from .boat_spatial_work import source_matrix, source_correction, subdivide


def trace(xy, center, radius):
    squared = np.sum((xy-center)**2, axis=-1)/radius**2
    values = np.zeros(squared.shape)
    inside = squared < 1
    values[inside] = np.exp(1-1/(1-squared[inside]))
    return values


def disk_quadrature(center, radius, order):
    """Independent polar area quadrature of the declared physical trace."""
    nodes, weights = leggauss(order)
    r = radius*(nodes+1)/2
    angle = np.arange(2*order)*np.pi/order
    xy = center+np.stack((r[:, None]*np.cos(angle),
                         r[:, None]*np.sin(angle)), axis=-1).reshape(-1, 2)
    area = np.repeat(radius*weights/2*r*np.pi/order, 2*order)
    return np.column_stack((xy, np.zeros(len(xy)))), area*trace(xy, center, radius)


def fourier_lift(points, normals, center, radius, size, length, origin):
    """Bounded direct sums; independent of the evolving driver's maps."""
    x, y = np.meshgrid(origin[0]+np.arange(size)*length/size,
                       origin[1]+np.arange(size)*length/size)
    values = trace(np.stack((x, y), axis=-1), center, radius)
    spectrum = np.fft.fft2(values).ravel()/(size*size)
    kx, ky = np.meshgrid(2*np.pi*np.fft.fftfreq(size, d=length/size),
                         2*np.pi*np.fft.fftfreq(size, d=length/size))
    kx, ky = kx.ravel(), ky.ravel()
    wave = np.hypot(kx, ky)
    potential, normal = np.empty(len(points)), np.empty(len(points))
    for first in range(0, len(points), 8):
        last = min(first+8, len(points))
        p, n = points[first:last], normals[first:last]
        phase = (p[:, 0, None]-origin[0])*kx+(p[:, 1, None]-origin[1])*ky
        basis = np.exp(1j*phase+p[:, 2, None]*wave)*spectrum
        potential[first:last] = basis.sum(axis=1).real
        normal[first:last] = (basis*(1j*n[:, 0, None]*kx
            +1j*n[:, 1, None]*ky+n[:, 2, None]*wave)).sum(axis=1).real
    return potential, normal, int(np.count_nonzero(values))


def run(inputs, refinements, orders):
    started = perf_counter()
    wet = hydrostatics(inputs.body, inputs.initial_state, density=WATER_DENSITY, gravity=GRAVITY)
    waterline = wet.waterline_segments[:, :, :2].reshape(-1, 2)
    radius = .5*float(np.min(np.ptp(waterline, axis=0)))
    center = waterline.mean(axis=0)
    center[1] = np.max(waterline[:, 1])+1.5*radius
    if radius <= 0 or center[1]-radius <= np.max(waterline[:, 1]):
        raise ValueError('strictly physical exterior disk required')
    rows = []
    for level in refinements:
        triangles = wet.wetted_triangles
        for _ in range(level):
            triangles = subdivide(triangles)
        panels = TrianglePanels(triangles, -1)
        modes = np.column_stack((panels.normals,
            np.cross(panels.centroids-inputs.initial_state.position, panels.normals)))
        single = assemble_dirichlet_halfspace(panels).single_layer
        matrix = source_matrix(panels)
        factor = lu_factor(matrix)
        density = lu_solve(factor, modes)
        radiated_potential = single@density
        weighted_modes = panels.areas[:, None]*modes
        previous = None
        for order in orders:
            potential, gradient = poisson_lift(panels.centroids, center, radius, order)
            normal = np.einsum('ij,ij->i', gradient, panels.normals)
            returned = potential-single@lu_solve(factor, normal)
            body_work = weighted_modes.T@returned
            green_work = radiated_potential.T@(panels.areas*normal)-weighted_modes.T@potential
            points, weights = disk_quadrature(center, radius, order)
            _, gradient_q = source_correction(points, panels, density, surface=True)
            surface_work = weights@gradient_q[:, :, 2]
            scale = np.maximum(np.abs(body_work), np.abs(surface_work))
            residual = body_work+surface_work
            row = dict(refinement=level, panel_count=len(panels), quadrature_order=order,
                panel_sha256=sha256(triangles.tobytes()).hexdigest(),
                body_work=body_work.tolist(), surface_work=surface_work.tolist(),
                green_work=green_work.tolist(),
                surface_minus_green=(surface_work-green_work).tolist(),
                pressure_minus_green=(-body_work-green_work).tolist(),
                residual=residual.tolist(), relative_to_mode_work=np.divide(residual, scale,
                    out=np.zeros_like(scale), where=scale != 0).tolist(),
                centroid_neumann_relative_residual=float(np.linalg.norm(matrix@density-modes)/np.linalg.norm(modes)))
            if previous is not None:
                row['quadrature_change'] = dict(body_work=(body_work-previous[0]).tolist(),
                                              surface_work=(surface_work-previous[1]).tolist())
            previous = body_work, surface_work
            # One body discretization isolates plane representation error. The
            # fixed 256m box matches the current time-domain diagnostic only.
            if level == refinements[0] and order == orders[-1]:
                row['fourier_controls'] = []
                length = 256.
                origin = inputs.initial_state.position[:2]-length/2
                for size in [129, 257, 513, 1025, 2049]:
                    if np.any(np.abs(center-origin-length/2)+radius >= length/2):
                        raise ValueError('physical trace support must fit strictly within FFT box')
                    h, n, active = fourier_lift(panels.centroids, panels.normals,
                        center, radius, size, length, origin)
                    work = weighted_modes.T@(h-single@lu_solve(factor, n))
                    x, y = np.meshgrid(origin[0]+np.arange(size)*length/size,
                                       origin[1]+np.arange(size)*length/size)
                    xy = np.column_stack((x.ravel(), y.ravel()))
                    values = trace(xy, center, radius)
                    xy, values = xy[values != 0], values[values != 0]
                    _, grid_gradient = source_correction(
                        np.column_stack((xy, np.zeros(len(xy)))), panels, density, surface=True)
                    grid_work = (length/size)**2*(values@grid_gradient[:, :, 2])
                    row['fourier_controls'].append(dict(grid_size=size, box_length_m=length,
                        nonzero_trace_nodes=active, body_work=work.tolist(),
                        harmonic_fill_equivalence_certified_by_support_gap=bool(.5*radius > length/size),
                        surface_rectangle_work=grid_work.tolist(),
                        surface_rectangle_change_from_disk=(grid_work-surface_work).tolist(),
                        body_work_change_from_aperiodic=(work-body_work).tolist(),
                        body_potential_relative_error=float(np.linalg.norm(h-potential)/np.linalg.norm(potential)),
                        body_normal_relative_error=float(np.linalg.norm(n-normal)/np.linalg.norm(normal))))
            rows.append(row)
            print(json.dumps(row), flush=True)
    return dict(complete=True, qualified=False,
        scope='Static aperiodic physical-exterior compact trace; no trajectory or throughput qualification',
        input_sha256=sha256(inputs.to_json().encode()).hexdigest(),
        modes=['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'],
        work_units=['m^2']*3+['m^3']*3,
        normalization='Normalized operator coefficients, not physical energy: dimensionless unit trace and rigid normal modes (dimensionless translations, metre rotations). No density multiplier.',
        identity='integral_body(q * total_fixed_body_potential) + integral_surface(trace * radiated_normal_velocity) = 0',
        bump_center_xy_m=center.tolist(), bump_radius_m=radius, support_y_clearance_m=.5*radius,
        elapsed_seconds=perf_counter()-started, rows=rows,
        limits=['P0 centroid hull traces and loads retain their raw collocation/quadrature error',
                'FFT controls change sampled plane representation and include periodic-image error',
                'FFT controls use the exact full compact trace, not runtime harmonic fill; a support gap greater than one grid spacing certifies identical zero hidden and fill-boundary samples',
                'A compact verification trace is not an incident sea-state component'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--refinements', type=int, nargs='+', default=[1, 2, 3])
    parser.add_argument('--orders', type=int, nargs='+', default=[24, 48])
    args = parser.parse_args()
    if (args.output.exists() or not args.refinements or not args.orders
            or any(x < 0 or x > 4 for x in args.refinements)
            or any(x < 4 for x in args.orders)):
        parser.error('new output, refinements0..4 and orders>=4 required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    sources = [*Path('src/topwater').glob('*.py'), Path('tests/laplace_reference.py'), Path(__file__),
        *[Path(__file__).with_name(name) for name in
          ['boat_hidden_extension.py', 'boat_spatial_work.py', 'panel_gradient_reference.py']]]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    result = run(SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8')),
                 args.refinements, args.orders)
    result['producer_sha256'] = hashes
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
