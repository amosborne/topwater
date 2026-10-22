"""Static P1-BEM / tetra-P4 interface accuracy falsifier on the same flat box.

No body, nonlinear motion, absorber, time stepping, or performance qualification.
The full-depth pocket is exactly pocket.grid/pocket.solve_case's partition.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from .operator import assemble_energy
from .pocket import grid, boundary, field
from . import p4_tetra


def shared_faces(left, right):
    identities = {tuple(sorted(face)) for face in right}
    mask = np.array([tuple(sorted(face)) in identities for face in left])
    return mask


def boundary_terms(points, faces, nodes, values, weights, samples):
    """Independent physical quadrature of q and consistent trace mass."""
    panels = TrianglePanels(points[faces], 1)
    _, gradient = field(samples)
    normal_flux = np.einsum('fqi,fi->fq', gradient, panels.normals)
    top = np.all(points[faces, 2] == 0, axis=1)
    local_load = panels.areas[:, None]*np.einsum('q,qi,fq->fi', weights, values, normal_flux)
    outer_load = np.zeros(len(points)); exact_flux = np.zeros(len(points))
    for i in range(nodes.shape[1]):
        np.add.at(outer_load, nodes[~top, i], local_load[~top, i])
        np.add.at(exact_flux, nodes[top, i], local_load[top, i])
    local_mass = panels.areas[top, None, None]*np.einsum('q,qi,qj->ij', weights, values, values)
    top_nodes = nodes[top]
    mass = coo_matrix((local_mass.ravel(),
                       (np.repeat(top_nodes, nodes.shape[1], axis=1).ravel(),
                        np.tile(top_nodes, (1, nodes.shape[1])).ravel())),
                      shape=(len(points), len(points))).tocsr()
    return mass, outer_load, exact_flux, dict(
        nodes=nodes, values=values, weights=weights, samples=samples,
        areas=panels.areas, normal_flux=normal_flux, top=top)


def solve_case(n, coupled, *, quadrature_order=12):
    points, tets, near, _ = grid(n, 0., 0)
    selected = tets[~near] if coupled else tets
    mesh = p4_tetra.make_mesh(points, selected)
    fem_energy = p4_tetra.stiffness(mesh)
    fem_faces = boundary(points, selected)
    local_nodes = local_faces = local_energy = interface = None
    local_load = np.zeros(len(mesh.points))
    if coupled:
        near_faces = boundary(points, tets[near])
        interface_mask = shared_faces(near_faces, fem_faces)
        interface = near_faces[interface_mask]
        local_nodes, inverse = np.unique(near_faces, return_inverse=True)
        local_faces = inverse.reshape(-1, 3)
        local_energy, bem_diagnostics = assemble_energy(
            local_faces, TrianglePanels(points[local_nodes][local_faces], 1), 3)
        restriction, retained = p4_tetra.restrict_trace(mesh, interface)
        physical_fem = fem_faces[~shared_faces(fem_faces, interface)]
        physical_local = near_faces[~interface_mask]
    else:
        restriction = eye(len(mesh.points), format='csr')
        retained = np.arange(len(mesh.points))
        physical_fem = fem_faces
        bem_diagnostics = None

    barycentric, weights = triangle_samples(np.eye(3)[None], quadrature_order)
    barycentric = barycentric[0]
    p4_values, _ = p4_tetra.basis(barycentric)
    fem_samples = barycentric@mesh.points[physical_fem]
    mass, outer_load, exact_flux, fem_part = boundary_terms(
        mesh.points, physical_fem, p4_tetra.face_nodes(mesh, physical_fem),
        p4_values, weights, fem_samples)
    parts = [fem_part]
    if coupled:
        extra_mass, local_load, extra_flux, local_part = boundary_terms(
            mesh.points, physical_local, physical_local, barycentric, weights,
            barycentric@mesh.points[physical_local])
        mass += extra_mass; exact_flux += extra_flux
        outer_load += local_load
        parts.append(local_part)

    energy_fem = (restriction.T@fem_energy@restriction).tocsr()
    energy = energy_fem.copy()
    if coupled:
        local_map = restriction[local_nodes]
        energy += (local_map.T@coo_matrix(local_energy)@local_map).tocsr()
    mass = (restriction.T@mass@restriction).tocsr()
    load = restriction.T@outer_load
    exact_flux = restriction.T@exact_flux
    active = np.flatnonzero(energy.diagonal() > 0)
    free = np.flatnonzero(mass.diagonal() > 0)
    unknown = np.setdiff1d(active, free)
    potential = np.zeros(len(retained))
    potential[free] = field(mesh.points[retained[free]])[0]
    right = load[unknown]-energy[unknown][:, free]@potential[free]
    potential[unknown] = splu(energy[unknown][:, unknown].tocsc()).solve(right)
    residual = energy@potential-load
    surface_mass = mass[free][:, free].tocsc()
    mass_factor = splu(surface_mass)
    flux = residual[free]
    flux_error = flux-exact_flux[free]
    projected_q = np.zeros(len(retained))
    projected_q[free] = mass_factor.solve(flux)
    weak_error = np.sqrt(flux_error@mass_factor.solve(flux_error)
                         / (exact_flux[free]@mass_factor.solve(exact_flux[free])))
    full_potential = restriction@potential
    full_q = restriction@projected_q
    potential_error = potential_norm = flux_squared_error = flux_norm = 0.
    for part in parts:
        values = part['values']; nodes = part['nodes']
        actual = np.einsum('qa,fa->fq', values, full_potential[nodes])
        exact = field(part['samples'])[0]
        physical_weights = part['areas'][:, None]*weights
        potential_error += np.sum(physical_weights*(actual-exact)**2)
        potential_norm += np.sum(physical_weights*exact**2)
        top = part['top']
        actual_q = np.einsum('qa,fa->fq', values, full_q[nodes])[top]
        exact_q = part['normal_flux'][top]
        flux_squared_error += np.sum(physical_weights[top]*(actual_q-exact_q)**2)
        flux_norm += np.sum(physical_weights[top]*exact_q**2)

    # |grad exp(kz)cos(kx*x+ky*y+phase)|² = k² exp(2kz).
    exact_energy = 8*np.pi**2*.7/2*(1-np.exp(-8*.7))
    computed_energy = float(potential@(energy@potential))
    trace_error = balance_error = None
    if coupled:
        interface_nodes = p4_tetra.face_nodes(mesh, interface)
        coarse_values = np.einsum('qa,fa->fq', barycentric, full_potential[interface])
        fine_values = np.einsum('qa,fa->fq', p4_values, full_potential[interface_nodes])
        trace_error = float(np.max(abs(coarse_values-fine_values)))
        # At shared free-surface nodes, the remaining load is physical FS flux;
        # balance the unknown interface rows, including known bottom Neumann data.
        interface_vertices = np.unique(interface)
        interface_unknown = np.intersect1d(unknown, interface_vertices)
        local_balance = local_map.T@(local_energy@(local_map@potential)) - restriction.T@local_load
        fem_balance = energy_fem@potential-restriction.T@(outer_load-local_load)
        balance_error = float(np.linalg.norm((local_balance+fem_balance)[interface_unknown])
                              / np.linalg.norm(local_balance[interface_unknown]))
    return dict(n=n, kind='P1-BEM/tetra-P4' if coupled else 'full tetra-P4',
        qualified=False, potential_dofs=len(active), unknown_dofs=len(unknown),
        free_surface_dofs=len(free), tetrahedra=len(selected),
        local_bem_faces=0 if not coupled else len(local_faces),
        eliminated_interface_trace_dofs=len(mesh.points)-len(retained),
        boundary_potential_relative_rms=float(np.sqrt(potential_error/potential_norm)),
        whole_free_surface_flux_relative_rms=float(np.sqrt(flux_squared_error/flux_norm)),
        projected_weak_free_surface_flux_relative_error=float(weak_error),
        energy_relative_error=float(abs(computed_energy/exact_energy-1)),
        energy=computed_energy, exact_energy=float(exact_energy),
        relative_unknown_residual=float(np.linalg.norm(residual[unknown])/np.linalg.norm(right)),
        shared_interface_potential_max_abs=trace_error,
        relative_shared_weak_flux_residual=balance_error,
        quadrature_order=quadrature_order, bem_diagnostics=bem_diagnostics,
        scope='Static harmonic accuracy only; no canoe, nonlinear evolution, outgoing '
              'waves, interface-scattering or performance qualification.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, nargs='+', default=[6, 12], choices=[6, 12])
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-pocket.json'))
    args = parser.parse_args()
    cases = []
    with threadpool_limits(1):
        for n in args.n:
            for coupled in (False, True):
                case = solve_case(n, coupled)
                cases.append(case)
                print(json.dumps(case), flush=True)
    args.output.parent.mkdir(exist_ok=True, parents=True)
    args.output.write_text(json.dumps(cases, indent=2)+'\n')
