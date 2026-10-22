"""Current whole-box tetra-P4 canoe inertia, with an unrestrained shared trace.

One static finite-domain calculation. No wave, nonlinear, convergence, speed,
or Capytaine qualification follows from this experiment.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from topwater.rankine import TrianglePanels
from .moving_canoe import make_domain, geometry, field
from .pocket import boundary
from .p4_canoe import tetrahedralize, harmonic_quartic, loads, diagnostics
from . import p4_tetra


def face_keys(faces):
    return {tuple(sorted(map(int, face))) for face in faces}


def make_volume(domain):
    points, _, _ = geometry(domain, domain.rest.position, domain.rest.orientation,
                             np.zeros(len(domain.free_nodes)))
    local = points[domain.local_mapping]
    markers = np.full(len(domain.contact.faces), 3)
    markers[domain.contact.body_faces] = 1
    markers[domain.contact.free_faces] = 2
    lp, lt, lf, lm, _ = tetrahedralize(local, domain.contact.faces, markers)
    # No coordinate tolerance or nearest-node matching: retain the exact source
    # vertex prefix and map its existing topological IDs into the global mesh.
    if not np.array_equal(lp[:len(local)], local):
        raise RuntimeError('TetGen changed the input vertex prefix')
    appended = np.arange(len(points), len(points)+len(lp)-len(local))
    mapping = np.r_[domain.local_mapping, appended]
    points = np.vstack((points, lp[len(local):]))
    local_tets, local_faces = mapping[lt], mapping[lf]
    interface = local_faces[lm == 3]
    exterior_faces = boundary(points, domain.exterior_tets)
    if not face_keys(interface).issubset(face_keys(exterior_faces)):
        raise RuntimeError('local/exterior triangular interface does not match exactly')
    tets = np.vstack((domain.exterior_tets, local_tets))
    faces = boundary(points, tets)
    expected = face_keys(exterior_faces)-face_keys(interface)
    expected |= face_keys(local_faces[lm != 3])
    if face_keys(faces) != expected:
        raise RuntimeError('merged physical boundary differs from the current domain')
    determinant = np.linalg.det(points[tets[:, 1:]]-points[tets[:, :1]])
    if not np.all(np.isfinite(determinant)) or np.any(determinant == 0):
        raise RuntimeError('invalid merged tetrahedral Jacobian')
    panels = TrianglePanels(points[faces], 1)
    boundary_volume = np.einsum('f,fi,fi->', panels.areas, panels.centroids, panels.normals)/3
    tetra_volume = np.sum(abs(determinant))/6
    if boundary_volume <= 0 or abs(tetra_volume/boundary_volume-1) > 1e-10:
        raise RuntimeError('merged tetrahedra disagree with the physical boundary volume')
    mesh = p4_tetra.make_mesh(points, tets)
    if len(mesh.points) > 100000:
        raise RuntimeError('merged mesh exceeded the 100000-P4-node preflight cap')
    interface_nodes = set(p4_tetra.face_nodes(mesh, interface).ravel())
    near_nodes = set(mesh.elements[len(domain.exterior_tets):].ravel())
    far_nodes = set(mesh.elements[:len(domain.exterior_tets)].ravel())
    if near_nodes & far_nodes != interface_nodes:
        raise RuntimeError('local/exterior P4 trace nodes are not exactly shared')
    free_faces = faces[np.all(points[faces, 2] == 0, axis=1)]
    body_faces = local_faces[lm == 1]
    preflight = dict(local_tetrahedra=len(lt), exterior_tetrahedra=len(domain.exterior_tets),
        interface_triangles=len(interface), shared_p4_trace_nodes=len(interface_nodes),
        eliminated_interface_modes=0, tetra_volume=float(tetra_volume),
        boundary_volume=float(boundary_volume),
        relative_volume_difference=float(abs(tetra_volume/boundary_volume-1)),
        minimum_absolute_tet_determinant=float(np.min(abs(determinant))))
    return mesh, faces, free_faces, body_faces, preflight


def run():
    domain = make_domain(12, False, 2)
    mesh, faces, free_faces, body_faces, preflight = make_volume(domain)
    operator = p4_tetra.stiffness(mesh)
    _, patch = loads(mesh, faces, domain.rest.position)
    polynomial, _ = harmonic_quartic(mesh.points, domain.rest.position)
    patch_error = float(np.linalg.norm(operator@polynomial-patch)/np.linalg.norm(patch))
    if patch_error > 1e-9:
        raise RuntimeError(f'merged-domain quartic patch failed: {patch_error:.12g}')
    load, _ = loads(mesh, body_faces, domain.rest.position)
    active = np.unique(mesh.elements)
    fixed = np.unique(p4_tetra.face_nodes(mesh, free_faces))
    unknown = np.setdiff1d(active, fixed)
    factor = splu(operator[unknown][:, unknown].tocsc())
    potential = np.zeros_like(load)
    potential[unknown] = factor.solve(load[unknown])
    actual = 1000*load.T@potential
    residual = float(np.linalg.norm((operator@potential-load)[unknown])/np.linalg.norm(load[unknown]))
    if residual > 1e-9:
        raise RuntimeError(f'whole-domain P4 residual failed: {residual:.12g}')
    energy = 1000*potential.T@(operator@potential)
    # Current unchanged nonlinear experiment, linear instantaneous auxiliary
    # problem on exactly the same finite tank/physical canoe boundary.
    count = len(domain.free_nodes)
    mixed = field(domain, domain.rest.position, domain.rest.orientation,
                  np.zeros(count), np.zeros(count), np.zeros(6))
    control = mixed['added_mass']
    d = np.array([1., 1., 1., 3.5, 3.5, 3.5])
    denominator = 1000*3.5**3*d[:, None]*d[None, :]
    current = mixed['points'][domain.local_mapping]
    wet = current[domain.contact.faces[domain.contact.body_faces]][:, ::-1]
    return dict(qualified=False, case='whole current finite tank, shared tetra-P4 canoe inertia',
        n=12, wet_subdivisions=2, density=1000., depth=4.,
        rotation_center=domain.rest.position.tolist(), tetrahedra=len(mesh.tets),
        active_p4_nodes=len(active), allocated_p4_nodes=len(mesh.points),
        orphan_nodes=len(mesh.points)-len(active), free_surface_p4_nodes=len(fixed),
        p4_unknowns=len(unknown), physical_boundary_faces=len(faces), wetted_faces=len(body_faces),
        current_wet_triangles_sha256=sha256(np.asarray(wet, dtype='<f8').tobytes()).hexdigest(),
        preflight=preflight, independent_harmonic_quartic_relative_patch_error=patch_error,
        relative_system_residual=residual,
        relative_energy_work_identity_error=float(np.linalg.norm(energy-actual)/np.linalg.norm(actual)),
        p4_inertia=actual.tolist(), mixed_p1_bem_fem_inertia=control.tolist(),
        p4_diagnostics=diagnostics(actual), mixed_diagnostics=diagnostics(control),
        diagonal_relative_difference=(np.diag(actual-control)/np.diag(control)).tolist(),
        complete_nondimensional_difference=((actual-control)/denominator).tolist(),
        scope='Single same-box depth4 static comparison, not a converged reference. '
              'No interface P1 restriction; all15 polynomial trace modes share each triangle. '
              'Capytaine omitted because the existing instantaneous helper assumes infinite '
              'depth and no lateral walls. No nonlinear, outgoing-wave or performance qualification.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-whole-current-canoe.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        report = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
