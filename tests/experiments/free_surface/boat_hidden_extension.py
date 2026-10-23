"""Aperiodic hidden-plane cancellation on supplied physical convex hulls.

The smooth plane trace is supported strictly inside the actual waterplane.
It changes no physical free-surface data. Its harmonic lift plus the fixed-body
Neumann correction must therefore vanish in the exterior fluid. This is a
static spatial falsifier, not an incident sea state or a dynamics benchmark.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.spatial import ConvexHull
from scipy.linalg import lu_factor, lu_solve
from scipy.linalg.lapack import dgecon

from topwater import SimulationInputs
from topwater.bem import assemble_dirichlet_halfspace
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from .boat_spatial_work import source_matrix, source_correction, subdivide


def poisson_lift(points, center, radius, order):
    """Analytic value/gradient kernels, integrated over a compact smooth disk."""
    if np.any((points[:, 2] == 0) & (np.linalg.norm(points[:, :2]-center, axis=1) <= radius)):
        raise ValueError('plane targets inside the hidden trace require a separate boundary limit')
    nodes, gauss = leggauss(order)
    r, radial_weights = radius*(nodes+1)/2, radius*gauss/2
    angles = np.arange(2*order)*np.pi/order
    sources = center+np.stack((r[:, None]*np.cos(angles),
                              r[:, None]*np.sin(angles)), axis=-1).reshape(-1, 2)
    bump = np.exp(1-1/(1-(r/radius)**2))
    weights = np.repeat(radial_weights*r*bump/(2*order), 2*order)
    value, gradient = np.empty(len(points)), np.empty((len(points), 3))
    for first in range(0, len(points), 32):
        last = min(first+32, len(points))
        xy = points[first:last, None, :2]-sources
        h = -points[first:last, 2, None]
        r2 = np.sum(xy*xy, axis=2)+h*h
        if np.any(h < 0) or np.any(r2 == 0):
            raise ValueError('Poisson evaluation requires z<=0 away from the trace support')
        inverse3 = r2**-1.5
        inverse5 = inverse3/r2
        value[first:last] = (h*inverse3)@weights
        gradient[first:last, :2] = np.einsum('tsj,ts,s->tj', xy, -3*h*inverse5, weights)
        gradient[first:last, 2] = (-inverse3+3*h*h*inverse5)@weights
    return value, gradient


def errors(lift, correction, weights):
    total = lift+correction
    norm = lambda a: float(np.sqrt(np.sum(weights*a*a)))
    scale = norm(lift)
    if scale == 0:
        raise ValueError('nonzero lift required for a cancellation control')
    return dict(relative_weighted_l2=norm(total)/scale,
                max_absolute=float(np.max(np.abs(total))), lift_weighted_l2=scale,
                correction_weighted_l2=norm(correction),
                max_relative_to_peak_lift=float(np.max(np.abs(total))/np.max(np.abs(lift))))


def targets(footprint, segments):
    lo = footprint.points.min(axis=0)-2.
    hi = footprint.points.max(axis=0)+2.
    x, y = np.meshgrid(np.linspace(lo[0], hi[0], 49), np.linspace(lo[1], hi[1], 33))
    xy = np.column_stack((x.ravel(), y.ravel()))
    planes = footprint.equations
    # These are bulk-field samples, not trace-limit quadrature. Keep a declared
    # 5 mm gap so grid alignment never silently places a target on a hull edge.
    distances_to_planes = (xy@planes[:, :2].T+planes[:, 2])/np.linalg.norm(planes[:, :2], axis=1)
    outside = np.any(distances_to_planes > .005, axis=1)
    grid = np.column_stack((xy[outside], np.zeros(np.sum(outside))))
    # Distances do not shrink with body refinement; include contact-near failure.
    distances = np.array([.01, .03, .1, .3])
    near = []
    center = footprint.points[footprint.vertices].mean(axis=0)
    for a, b in segments:
        point = .63*a+.37*b
        edge = b[:2]-a[:2]
        normal = np.array([edge[1], -edge[0]])/np.linalg.norm(edge)
        if normal@(point[:2]-center) < 0:
            normal *= -1
        for d in distances:
            near.append(point+np.array([*(d*normal), 0.]))
    return grid, np.asarray(near), distances


def run(inputs, refinements, orders):
    started = perf_counter()
    wet = hydrostatics(inputs.body, inputs.initial_state, density=WATER_DENSITY, gravity=GRAVITY)
    footprint = ConvexHull(wet.waterline_segments[:, :, :2].reshape(-1, 2))
    if not np.isclose(footprint.volume, wet.waterplane_area, rtol=1e-12, atol=1e-12):
        raise ValueError('waterplane is not the required convex footprint')
    center = footprint.points[footprint.vertices].mean(axis=0)
    plane_norm = np.linalg.norm(footprint.equations[:, :2], axis=1)
    clearance = np.min(-(footprint.equations[:, :2]@center+footprint.equations[:, 2])/plane_norm)
    radius = .6*clearance
    if radius <= 0:
        raise ValueError('strictly interior disk required')
    if np.any(footprint.equations[:, :2]@center+footprint.equations[:, 2]+radius*plane_norm >= 0):
        raise ValueError('hidden trace support crosses the actual waterline')
    grid, near, distances = targets(footprint, wet.waterline_segments)
    surface = np.vstack((grid, near))
    if np.any(np.linalg.norm(surface[:, :2]-center, axis=1) <= radius):
        raise ValueError('physical surface targets must be outside the hidden trace support')
    surface_lifts = [poisson_lift(surface, center, radius, order)[1][:, 2] for order in orders]
    rows = []
    for level in refinements:
        triangles = wet.wetted_triangles
        for _ in range(level):
            triangles = subdivide(triangles)
        panels = TrianglePanels(triangles, -1)
        before = perf_counter()
        op = assemble_dirichlet_halfspace(panels)
        matrix = source_matrix(panels)
        body_lifts = [poisson_lift(panels.centroids, center, radius, order) for order in orders]
        q = np.column_stack([-np.einsum('tj,tj->t', grad, panels.normals) for _, grad in body_lifts])
        factor, pivots = lu_factor(matrix)
        reciprocal_condition, info = dgecon(factor, np.linalg.norm(matrix, 1))
        if info:
            raise RuntimeError('condition estimate failed')
        density = lu_solve((factor, pivots), q)
        potential = op.single_layer@density
        modes = np.column_stack((panels.normals,
            np.cross(panels.centroids-inputs.initial_state.position, panels.normals)))
        _, gradient = source_correction(surface, panels, density, surface=True)
        flux = gradient[:, :, 2]
        for col, order in enumerate(orders):
            body_phi = body_lifts[col][0]
            lifted_flux = surface_lifts[col]
            grid_weights = np.ones(len(grid))
            near_weights = np.ones(len(near))
            row = dict(refinement=level, panel_count=len(panels), lift_order=order,
                panel_sha256=sha256(triangles.tobytes()).hexdigest(),
                setup_and_evaluation_seconds=perf_counter()-before,
                matrix_reciprocal_condition_1=float(reciprocal_condition),
                max_abs_prescribed_normal_velocity=float(np.max(np.abs(q[:, col]))),
                max_abs_density=float(np.max(np.abs(density[:, col]))),
                potential_moments=dict(
                    modes=['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'],
                    units=['kg*m/s']*3+['kg*m^2/s']*3,
                    meaning='rho integral of fluid-normal rigid mode times potential; impulse units, not static force',
                    lift=(WATER_DENSITY*modes.T@(panels.areas*body_phi)).tolist(),
                    correction=(WATER_DENSITY*modes.T@(panels.areas*potential[:, col])).tolist(),
                    total=(WATER_DENSITY*modes.T@(panels.areas*(body_phi+potential[:, col]))).tolist()),
                body_potential=errors(body_phi, potential[:, col], panels.areas),
                exterior_grid_flux=errors(lifted_flux[:len(grid)], flux[:len(grid), col], grid_weights),
                contact_near_flux=errors(lifted_flux[len(grid):], flux[len(grid):, col], near_weights),
                contact_distances_m=distances.tolist(),
                contact_flux_residual_by_segment=(lifted_flux[len(grid):]+flux[len(grid):, col]).reshape(-1, len(distances)).tolist(),
                linear_relative_residual=float(np.linalg.norm(matrix@density[:, col]-q[:, col])/np.linalg.norm(q[:, col])))
            if col:
                row['lift_change_from_previous_order'] = dict(
                    body_potential_relative=float(np.linalg.norm(body_phi-body_lifts[col-1][0])/np.linalg.norm(body_phi)),
                    body_normal_derivative_relative=float(np.linalg.norm(q[:, col]-q[:, col-1])/np.linalg.norm(q[:, col])),
                    exterior_flux_relative=float(np.linalg.norm(lifted_flux-surface_lifts[col-1])/np.linalg.norm(lifted_flux)))
            rows.append(row)
            print(json.dumps(row), flush=True)
    return dict(complete=True, qualified=False,
        scope='Static hidden-plane gauge/extension falsifier on actual supplied hull; no dynamics or throughput claim.',
        exact_exterior_total_field='zero potential and zero velocity',
        input_sha256=sha256(inputs.to_json().encode()).hexdigest(),
        bump_center_xy_m=center.tolist(), bump_radius_m=radius,
        support_clearance_m=clearance-radius, exterior_grid_target_count=len(grid),
        exterior_grid_minimum_plane_gap_m=.005,
        near_rim_target_count=len(near), elapsed_seconds=perf_counter()-started, rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--refinements', type=int, nargs='+', default=[0, 1, 2, 3])
    parser.add_argument('--orders', type=int, nargs='+', default=[24, 48])
    args = parser.parse_args()
    if args.output.exists() or any(r < 0 or r > 4 for r in args.refinements) or any(o < 4 for o in args.orders):
        parser.error('new output, refinement 0..4 and lift order >=4 required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    sources = [*Path('src/topwater').glob('*.py'), Path('tests/laplace_reference.py'), Path(__file__),
        Path(__file__).with_name('boat_spatial_work.py'), Path(__file__).with_name('panel_gradient_reference.py')]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    inputs = SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8'))
    result = run(inputs, args.refinements, args.orders)
    result['producer_sha256'] = hashes
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
