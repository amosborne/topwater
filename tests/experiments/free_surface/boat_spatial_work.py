"""Actual-hull Green-work falsifier for a spatial surface/body split.

Static verification only: no evolving water, operator reuse or timing claim.
A harmonic plane-wave lift and a zero-plane-potential body correction share
the physical wetted hull. The exposed plane, sides and bottom close one finite
fluid volume. Raw P0 collocation is retained, including its reciprocity error.
"""
import argparse
from functools import partial
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.spatial import ConvexHull

from topwater import SimulationInputs
from topwater.bem import assemble_dirichlet_halfspace
from topwater.boundary import BoundarySystem
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.laplace_reference import point_source
from .panel_gradient_reference import checked_influence


def subdivide(triangles):
    """Coplanar refinement, independent of physical hull topology."""
    a, b, c = triangles.transpose(1, 0, 2)
    ab, bc, ca = (a+b)/2, (b+c)/2, (c+a)/2
    return np.concatenate([np.stack(t, axis=1) for t in
                           ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))])


def clip_polygon(polygon, plane, inside):
    """Float clipping for integration cells only; never solver contact geometry."""
    result = []
    if not len(polygon):
        return np.empty((0, 2))
    values = polygon@plane[:2]+plane[2]
    if not inside:
        values = -values
    for a, b, va, vb in zip(polygon, np.roll(polygon, -1, axis=0),
                            values, np.roll(values, -1)):
        if va <= 0:
            result.append(a)
        if (va < 0 < vb) or (vb < 0 < va):
            result.append(a+(b-a)*(va/(va-vb)))
    return np.asarray(result).reshape(-1, 2)


def box_boundary(footprint, half_width, depth):
    """Disjoint convex exposed-plane pieces plus four sides and the bottom."""
    h = half_width
    square = np.array([[-h, -h], [h, -h], [h, h], [-h, h]])
    if np.any(np.abs(footprint.points[footprint.vertices]) >= h):
        raise ValueError('integration box must enclose the actual waterplane')
    remaining = square
    pieces = []
    for plane in footprint.equations:
        outside = clip_polygon(remaining, plane, False)
        if len(outside) >= 3:
            pieces.append(np.column_stack((outside, np.zeros(len(outside)))))
        remaining = clip_polygon(remaining, plane, True)
    surface = [np.array([p[0], p[j], p[j+1]]) for p in pieces
               for j in range(1, len(p)-1)]
    cap = np.column_stack((square, np.full(4, -depth)))
    quads = [cap]
    normals = [np.array([0., 0., -1.])]
    for a, b in zip(square, np.roll(square, -1, axis=0)):
        quads.append(np.array([[*a, 0.], [*b, 0.], [*b, -depth], [*a, -depth]]))
        edge = b-a
        normals.append(np.array([edge[1], -edge[0], 0.])/np.linalg.norm(edge))
    outer = np.array([q[list(ids)] for q in quads for ids in ((0, 1, 2), (0, 2, 3))])
    return np.array(surface), outer, np.repeat(normals, 2, axis=0)


def lift(points, wavevector):
    """Exact harmonic cos/sin modes, independent of the numerical BEM."""
    k = np.linalg.norm(wavevector)
    phase = points[:, :2]@wavevector
    scale = np.exp(k*points[:, 2])
    values = scale[:, None]*np.column_stack((np.cos(phase), np.sin(phase)))
    tangential = scale[:, None]*np.column_stack((-np.sin(phase), np.cos(phase)))
    gradients = np.empty((len(points), 2, 3))
    gradients[:, :, :2] = tangential[:, :, None]*wavevector
    gradients[:, :, 2] = k*values
    return values, gradients


def correction(points, panels, q, potential, *, surface=False):
    """Contract blocked physical layer potentials without storing the full grid."""
    values = np.zeros((len(points), q.shape[1]))
    gradients = np.empty((len(points), q.shape[1], 3))
    for first in range(0, len(points), 128):
        last = min(first+128, len(points))
        targets = points[first:last]
        real = checked_influence(targets, panels)
        if surface:
            gradients[first:last, :, :2] = 0.
            gradients[first:last, :, 2] = 2*(real.gradient_single_layer[:, :, 2]@q
                -real.gradient_double_layer[:, :, 2]@potential)
        else:
            mirror = targets*np.array([1., 1., -1.])
            image = checked_influence(mirror, panels)
            values[first:last] = ((real.single_layer-image.single_layer)@q
                                 -(real.double_layer-image.double_layer)@potential)
            gradients[first:last] = (np.einsum('tnj,nm->tmj',
                real.gradient_single_layer-image.gradient_single_layer*np.array([1., 1., -1.]), q)
                -np.einsum('tnj,nm->tmj', real.gradient_double_layer
                    -image.gradient_double_layer*np.array([1., 1., -1.]), potential))
    return values, gradients


def source_matrix(panels):
    """Target-normal derivative of S_D, including the fluid-side +1/2 jump."""
    count = len(panels)
    matrix = np.empty((count, count))
    for row in range(count):
        other = np.arange(count) != row
        subset = TrianglePanels(panels.triangles[other], panels.normal_sign[other])
        real = checked_influence(panels.centroids[row:row+1], subset)
        matrix[row, other] = real.gradient_single_layer[0]@panels.normals[row]
        matrix[row, row] = .5  # Real-panel principal normal derivative is zero.
    mirror = panels.centroids*np.array([1., 1., -1.])
    for first in range(0, count, 128):
        last = min(first+128, count)
        image = checked_influence(mirror[first:last], panels)
        matrix[first:last] -= np.einsum('tnj,tj->tn', image.gradient_single_layer,
                            panels.normals[first:last]*np.array([1., 1., -1.]))
    return matrix


def source_correction(points, panels, density, *, surface=False):
    """Continuous odd-image single-layer field; P0 unknowns are densities."""
    values = np.zeros((len(points), density.shape[1]))
    gradients = np.zeros((len(points), density.shape[1], 3))
    for first in range(0, len(points), 128):
        last = min(first+128, len(points))
        targets = points[first:last]
        real = checked_influence(targets, panels)
        if surface:
            gradients[first:last, :, 2] = 2*real.gradient_single_layer[:, :, 2]@density
        else:
            image = checked_influence(targets*np.array([1., 1., -1.]), panels)
            values[first:last] = (real.single_layer-image.single_layer)@density
            gradients[first:last] = np.einsum('tnj,nm->tmj', real.gradient_single_layer
                -image.gradient_single_layer*np.array([1., 1., -1.]), density)
    return values, gradients


def integrate(triangles, order):
    points, weights = triangle_samples(triangles, order)
    area = np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0],
                                   triangles[:, 2]-triangles[:, 0]), axis=1)/2
    return points.reshape(-1, 3), (area[:, None]*weights).ravel()


def contact_profile(segments, footprint, field):
    """Resolve the actual contact edge, away from triangle vertices."""
    segment = segments[np.argmax(np.linalg.norm(segments[:, 1]-segments[:, 0], axis=1))]
    base = .63*segment[0]+.37*segment[1]
    edge = segment[1, :2]-segment[0, :2]
    outward = np.array([edge[1], -edge[0]])/np.linalg.norm(edge)
    center = footprint.points[footprint.vertices].mean(axis=0)
    if outward@(base[:2]-center) < 0:
        outward = -outward
    distances = np.array([.01, .003, .001, .0003, .0001, .00003, .00001])
    points = np.broadcast_to(base, (len(distances), 3)).copy()
    points[:, :2] += distances[:, None]*outward
    _, gradient = field(points, surface=True)
    flux = gradient[:, :, 2]
    return dict(edge_point=base.tolist(), outward_xy=outward.tolist(),
                distance_m=distances.tolist(), plane_flux=flux.tolist(),
                distance_times_flux=(distances[:, None]*flux).tolist())


def run(inputs, *, refinements, orders, half_widths, depth, frequency, heading):
    started = perf_counter()
    wet = hydrostatics(inputs.body, inputs.initial_state, density=WATER_DENSITY, gravity=GRAVITY)
    footprint = ConvexHull(wet.waterline_segments[:, :, :2].reshape(-1, 2))
    wavevector = (2*np.pi*frequency)**2/GRAVITY*np.array([np.cos(heading), np.sin(heading)])
    rows = []
    for refinement in refinements:
        triangles = wet.wetted_triangles
        for _ in range(refinement):
            triangles = subdivide(triangles)
        panels = TrianglePanels(triangles, -1)
        modes = np.column_stack((panels.normals[:, 2],
            np.cross(panels.centroids-inputs.initial_state.position, panels.normals)[:, 0]))
        p, gradient = lift(panels.centroids, wavevector)
        n = np.einsum('fmj,fj->fm', gradient, panels.normals)
        before = perf_counter()
        op = assemble_dirichlet_halfspace(panels)
        system = BoundarySystem(op, dirichlet=np.zeros(len(panels), bool),
            areas=panels.areas, length_scale=float(np.ptp(triangles[:, :, 0])), bounded=False)
        solution = system.solve(normal_velocity=np.column_stack((modes, -n)))
        setup_direct = perf_counter()-before
        before = perf_counter()
        matrix = source_matrix(panels)
        rhs = np.column_stack((modes, -n))
        density = np.linalg.solve(matrix, rhs)
        source_phi = op.single_layer@density
        source_residual = float(np.linalg.norm(matrix@density-rhs)/np.linalg.norm(rhs))
        exact_phi, exact_grad = point_source(panels.centroids, wet.center_of_buoyancy)
        mirror_phi, mirror_grad = point_source(panels.centroids,
            wet.center_of_buoyancy*np.array([1., 1., -1.]))
        control_q = np.einsum('tmj,tj->tm', exact_grad-mirror_grad, panels.normals)
        control_density = np.linalg.solve(matrix, control_q)
        exact_phi -= mirror_phi
        control_error = float(np.linalg.norm(op.single_layer@control_density-exact_phi)
                              /np.linalg.norm(exact_phi))
        setup_source = perf_counter()-before+op.diagnostics['assembly_seconds']
        trials = [
            ('direct_p0_potential', solution.potential,
             partial(correction, panels=panels, q=modes, potential=solution.potential[:, :2]),
             setup_direct, solution.diagnostics['original_bie_residual']),
            ('indirect_p0_density', source_phi,
             partial(source_correction, panels=panels, density=density[:, :2]),
             setup_source, source_residual)]
        for name, phi, field, setup_seconds, linear_residual in trials:
            body = phi[:, :2].T@(panels.areas[:, None]*n)-modes.T@(panels.areas[:, None]*p)
            inertia = WATER_DENSITY*modes.T@(panels.areas[:, None]*phi[:, :2])
            profile = contact_profile(wet.waterline_segments, footprint, field)
            for half_width in half_widths:
                if (np.any(np.abs(triangles[:, :, :2]) >= half_width)
                        or np.any(triangles[:, :, 2] <= -depth)):
                    raise ValueError('integration side/bottom must strictly enclose the wetted hull')
                surface, outer, normals = box_boundary(footprint, half_width, depth)
                for order in orders:
                    before = perf_counter()
                    points, weights = integrate(surface, order)
                    trace, _ = lift(points, wavevector)
                    _, rg = field(points, surface=True)
                    work = rg[:, :, 2].T@(weights[:, None]*trace)
                    top_area = float(weights.sum())
                    points, weights = integrate(outer, order)
                    trace, hg = lift(points, wavevector)
                    rp, rg = field(points)
                    normal = np.repeat(normals, order*order, axis=0)
                    hn = np.einsum('tmj,tj->tm', hg, normal)
                    rn = np.einsum('tmj,tj->tm', rg, normal)
                    tail = rp.T@(weights[:, None]*hn)-rn.T@(weights[:, None]*trace)
                    residual = work-body-tail
                    scale = np.abs(work)+np.abs(body)+np.abs(tail)
                    rows.append(dict(representation=name, refinement=refinement, panel_count=len(panels), order=order,
                        panel_sha256=sha256(triangles.tobytes()).hexdigest(),
                        contact_profile=profile,
                        half_width_m=half_width, depth_m=depth, setup_seconds=setup_seconds,
                        integration_seconds=perf_counter()-before, exposed_area_m2=top_area,
                        expected_exposed_area_m2=4*half_width**2-footprint.volume,
                        surface_work=work.tolist(), body_work=body.tolist(), outer_work=tail.tolist(),
                        residual=residual.tolist(), relative_residual=(residual/np.maximum(scale, 1e-14)).tolist(),
                        added_mass_heave_roll=inertia.tolist(),
                        body_potential_from_lift_rms=np.sqrt(np.mean((p+phi[:, 2:])**2, axis=0)).tolist(),
                    linear_residual=linear_residual,
                    source_manufactured_potential_relative_error=control_error))
                    print(json.dumps(rows[-1]), flush=True)
    return dict(complete=True, qualified=False, scope='Static actual-hull finite-box Green identity; no dynamics or speed qualification',
        input_sha256=sha256(inputs.to_json().encode()).hexdigest(), frequency_hz=frequency,
        heading_radians=heading, modes=['Heave', 'Roll'], lift_modes=['cos', 'sin'],
        identity='surface_work = body_work + outer_work; fluid-outward hull normals',
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=perf_counter()-started, rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--refinements', type=int, nargs='+', default=[0, 1, 2])
    parser.add_argument('--orders', type=int, nargs='+', default=[8, 16])
    parser.add_argument('--half-widths', type=float, nargs='+', default=[12., 24.])
    parser.add_argument('--depth', type=float, default=8.)
    parser.add_argument('--frequency', type=float, default=.1)
    parser.add_argument('--heading', type=float, default=.7)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output already exists')
    if (any(r < 0 or r > 3 for r in args.refinements) or any(o < 2 for o in args.orders)
            or not np.isfinite(args.depth) or args.depth <= 0
            or not np.isfinite(args.heading) or not 0 < args.frequency <= .6
            or any(not np.isfinite(h) or h <= 0 for h in args.half_widths)):
        parser.error('refinement 0..3, order >=2, positive depth and frequency <=.6 Hz required')
    inputs = SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8'))
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    sources = [*Path('src/topwater').glob('*.py'), Path(__file__),
               Path(__file__).with_name('panel_gradient_reference.py'), Path('tests/laplace_reference.py')]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    result = run(inputs, refinements=args.refinements, orders=args.orders,
                 half_widths=args.half_widths, depth=args.depth,
                 frequency=args.frequency, heading=args.heading)
    result['producer_sha256'] = hashes
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
