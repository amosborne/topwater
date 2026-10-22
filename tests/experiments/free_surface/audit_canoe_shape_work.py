"""Geometry-only kinetic-energy shape audit of the nonlinear mixed canoe solve.

Continuum fixed-material-trace identity:
  .5*rho*phi.T*E_dot*phi = rho*integral_boundary(
      .5*|u|^2*w.n - (w.u)*(u.n)) dS.
The left side is isolated from fresh geometry-only field solves, subtracting
trace variation. Physical-boundary quadrature uses reconstructed velocities.
A discrepancy also contains spatial/FEM interior and artificial-interface work;
it does not by itself identify a pressure sign error or an incorrect FD step.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from tests.laplace_reference import triangle_samples
from .audit_canoe_energy import kinetic_and_load
from .moving_canoe import make_domain, initial_state, field, surface, moved_orientation

RHO = 1000.


def boundary_work(domain, current, position, psi, velocity, height_rate, mesh_velocity):
    bary, weights = triangle_samples(np.eye(3)[None], 6)
    bary = bary[0]
    points = current['points'][domain.local_mapping]
    faces = domain.contact.faces[domain.contact.body_faces]
    triangles = points[faces]
    normal_area = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    twice_area = np.linalg.norm(normal_area, axis=1)
    normals = normal_area/twice_area[:, None]
    opposite = np.roll(triangles, -2, axis=1)-np.roll(triangles, -1, axis=1)
    gradient = np.cross(normals[:, None], opposite)/twice_area[:, None, None]
    tangent = np.einsum('fci,fc->fi', gradient,
        current['potential'][domain.local_mapping][faces])
    samples = np.einsum('qc,fci->fqi', bary, triangles)
    rigid = velocity[:3]+np.cross(velocity[3:], samples-position)
    normal_velocity = np.einsum('fqi,fi->fq', rigid, normals)
    fluid = tangent[:, None]+normal_velocity[:, :, None]*normals[:, None]
    moving = np.einsum('qc,fci->fqi', bary, mesh_velocity[domain.local_mapping][faces])
    moving_normal = np.einsum('fqi,fi->fq', moving, normals)
    integrand = .5*np.sum(fluid*fluid, axis=2)*moving_normal-np.sum(moving*fluid, axis=2)*normal_velocity
    hull = RHO*np.einsum('f,q,fq->', twice_area/2, weights, integrand)

    faces = domain.free_faces
    triangles = current['points'][faces]
    normal_area = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    area = normal_area[:, 2]/2
    slope = -normal_area[:, :2]/normal_area[:, 2, None]
    opposite = np.roll(triangles[:, :, :2], -2, axis=1)-np.roll(triangles[:, :, :2], -1, axis=1)
    gradient = np.stack((-opposite[:, :, 1], opposite[:, :, 0]), axis=-1)/(2*area[:, None, None])
    lookup = np.full(len(current['points']), -1, dtype=int)
    lookup[domain.free_nodes] = np.arange(len(domain.free_nodes))
    indices = lookup[faces]
    grad_psi = np.einsum('fci,fc->fi', gradient, psi[indices])
    moving = np.einsum('qc,fci->fqi', bary, mesh_velocity[faces])
    eta_rate = np.einsum('qc,fc->fq', bary, height_rate[indices])
    # This is the P1 weak-flux recovery used by the undamped evolution, not
    # an independent pointwise Neumann trace or an exact BEM gradient limit.
    graph_flux = eta_rate-np.einsum('fqi,fi->fq', moving[:, :, :2], slope)
    vertical = (graph_flux+np.sum(slope*grad_psi, axis=1)[:, None])/(1+np.sum(slope*slope, axis=1))[:, None]
    horizontal = grad_psi[:, None]-vertical[:, :, None]*slope[:, None]
    fluid = np.concatenate((horizontal, vertical[:, :, None]), axis=2)
    integrand = graph_flux*(.5*np.sum(fluid*fluid, axis=2)-np.sum(moving*fluid, axis=2))
    free = RHO*np.einsum('f,q,fq->', area, weights, integrand)
    return float(hull), float(free)


def fem_volume_shape(domain, current, mesh_velocity):
    """Exact P1 volume shape term at fixed nodal potential, not boundary work."""
    tets = domain.exterior_tets
    coordinates = current['points'][tets]
    edges = coordinates[:, 1:]-coordinates[:, :1]
    gradients = np.empty((len(tets), 4, 3))
    gradients[:, 1:] = np.linalg.inv(edges).transpose(0, 2, 1)
    gradients[:, 0] = -gradients[:, 1:].sum(axis=1)
    u = np.einsum('tci,tc->ti', gradients, current['potential'][tets])
    dw = np.einsum('tcj,tci->tij', gradients, mesh_velocity[tets])
    volume = np.abs(np.linalg.det(edges))/6
    return float(RHO*np.sum(volume*(.5*np.sum(u*u, axis=1)*np.trace(dw, axis1=1, axis2=2)-
        np.einsum('ti,tij,tj->t', u, dw, u))))


def run(n=12, wet_subdivisions=2, initial_scale=10., gpu=True):
    domain = make_domain(n, gpu, wet_subdivisions)
    state = initial_state(domain, scale=initial_scale)
    size = len(domain.free_nodes)
    position, orientation, velocity = state[:3], state[3:7], state[7:13]
    eta, psi = state[13:13+size], state[13+size:]
    current = field(domain, position, orientation, eta, psi, velocity)
    height_rate, _, moving = surface(domain, current, psi, velocity)
    kinetic, load = kinetic_and_load(domain, current, position, velocity)
    hull, free = boundary_work(domain, current, position, psi, velocity, height_rate, moving)
    fe_shape = fem_volume_shape(domain, current, moving)
    rows = []
    for h in (2e-4, 1e-4):
        solutions = []
        energies = []
        for dt in (h, -h):
            moved_position = position+dt*velocity[:3]
            solution = field(domain, moved_position, moved_orientation(orientation, velocity[3:], dt),
                             eta+dt*height_rate, psi, velocity)
            solutions.append(solution)
            energies.append(kinetic_and_load(domain, solution, moved_position, velocity)[0])
        kinetic_rate = (energies[0]-energies[1])/(2*h)
        trace_rate = (solutions[0]['potential']-solutions[1]['potential'])/(2*h)
        trace_work = RHO*trace_rate@load
        shape = kinetic_rate-trace_work
        rows.append(dict(directional_step_seconds=h, kinetic_geometry_derivative_watts=float(kinetic_rate),
            trace_variation_work_watts=float(trace_work), operator_shape_work_watts=float(shape),
            physical_boundary_shape_work_watts=hull+free,
            shape_minus_boundary_watts=float(shape-hull-free),
            bem_operator_shape_work_watts=float(shape-fe_shape)))
    return dict(qualified=False, n=n, wet_subdivisions=wet_subdivisions, initial_scale=initial_scale,
        gpu=gpu, density=RHO, kinetic_energy_joules=kinetic,
        hull_boundary_shape_work_watts=hull, free_surface_boundary_shape_work_watts=free,
        fem_volume_shape_work_watts=fe_shape, results=rows,
        scope='Undamped actual geometry-rate direction; psi and world twist held fixed in offsets. '
              'Boundary u uses exact rigid normal velocity on hull plus P1 tangential gradient, '
              'and the evolution weak-flux reconstruction on FS. Physical rigid walls have zero '
              'normal velocity/flux. FEM interior shape work and artificial BEM/FEM interface '
              'trace/gradient mismatch can contribute; a residual does not alone identify a '
              'pressure bug. Finite-difference and spatial convergence remain unqualified.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=12)
    parser.add_argument('--wet-subdivisions', type=int, default=2)
    parser.add_argument('--initial-scale', type=float, default=10.)
    parser.add_argument('--cpu', action='store_true')
    args = vars(parser.parse_args()); args['gpu'] = not args.pop('cpu')
    with threadpool_limits(1):
        report = run(**args)
    path = Path('artifacts'); path.mkdir(exist_ok=True)
    suffix = f"n{args['n']}-wet{args['wet_subdivisions']}-initial{args['initial_scale']:g}"
    (path/f'canoe-shape-work-{suffix}.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
