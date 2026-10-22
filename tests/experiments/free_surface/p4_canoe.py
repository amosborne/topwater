"""Bounded current-canoe P4 inertia falsifier, not an unbounded A-infinity gate.

The exact current moving-canoe local boundary is meshed afresh. Its finite
side walls and one-metre bed differ from Capytaine's unbounded exterior.
MeshPy/TetGen is optional preprocessing only; no project dependency is added.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from .canoe import body_load
from .moving_canoe import make_domain, geometry
from .operator import assemble_energy
from .pocket import boundary
from . import p4_tetra


def triangle_keys(triangles):
    return {tuple(sorted(tuple(point) for point in triangle)) for triangle in triangles}


def tetrahedralize(points, faces, markers):
    from meshpy.tet import MeshInfo, Options, build
    info = MeshInfo()
    info.set_points(points)
    info.set_facets(faces, markers=markers.tolist())
    # Y forbids boundary facet subdivision. Reject changed triangles below;
    # quality or point-count failure must not silently change the canoe.
    result = build(info, options=Options('pq1.4YfQ'))
    coordinates = np.asarray(result.points)
    tets = np.asarray(result.elements)
    if len(tets) > 6000:
        raise RuntimeError(f'mesher exceeded the 6000-tetrahedron preflight cap: {len(tets)}')
    determinant = np.linalg.det(coordinates[tets[:, 1:]]-coordinates[tets[:, :1]])
    if not np.all(np.isfinite(determinant)) or np.any(determinant == 0):
        raise RuntimeError('invalid tetrahedral Jacobian')
    negative = determinant < 0
    tets[negative, 0], tets[negative, 1] = tets[negative, 1].copy(), tets[negative, 0].copy()
    actual_faces = boundary(coordinates, tets)
    marker_by_face = {tuple(sorted(face)): marker for face, marker in
                      zip(result.faces, result.face_markers) if marker}
    actual_markers = np.array([marker_by_face[tuple(sorted(face))] for face in actual_faces])
    for marker in (1, 2, 3):
        if triangle_keys(coordinates[actual_faces[actual_markers == marker]]) != triangle_keys(points[faces[markers == marker]]):
            raise RuntimeError(f'mesher changed the current boundary, marker {marker}')
    panels = TrianglePanels(points[faces], 1)
    boundary_volume = np.einsum('f,fi,fi->', panels.areas, panels.centroids, panels.normals)/3
    tetra_volume = np.abs(determinant).sum()/6
    if boundary_volume <= 0 or abs(tetra_volume-boundary_volume) > 1e-10*boundary_volume:
        raise RuntimeError('tetrahedral volume disagrees with the exact current boundary')
    return coordinates, tets, actual_faces, actual_markers, dict(
        tetra_volume=float(tetra_volume), boundary_volume=float(boundary_volume),
        relative_volume_difference=float(abs(tetra_volume/boundary_volume-1)),
        minimum_absolute_tet_determinant=float(np.min(abs(determinant))))


def harmonic_quartic(points, center):
    # Re[(u.x+i*v.x)^4] is harmonic for orthonormal u,v. Both vectors together
    # involve all world coordinates, so the patch exercises every gradient axis.
    u = np.array([.6, .8, 0.]); v = np.array([-.48, .36, .8])
    value = ((points-center)@u+1j*((points-center)@v))/3.5
    return value.real**4-6*value.real**2*value.imag**2+value.imag**4, (
        4*np.real(value[..., None]**3*(u+1j*v))/3.5)


def loads(mesh, faces, center):
    """Exact P4 rigid-mode loads and independent harmonic-quartic boundary load."""
    panels = TrianglePanels(mesh.points[faces], 1)
    samples, weights = triangle_samples(panels.triangles, 5)
    barycentric, _ = triangle_samples(np.eye(3)[None], 5)
    values, _ = p4_tetra.basis(barycentric[0])
    nodes = p4_tetra.face_nodes(mesh, faces)
    modes = np.concatenate((np.broadcast_to(panels.normals[:, None], samples.shape),
                            np.cross(samples-center, panels.normals[:, None])), axis=2)
    _, gradient = harmonic_quartic(samples, center)
    flux = np.einsum('fqi,fi->fq', gradient, panels.normals)
    load = np.zeros((len(mesh.points), 6)); patch = np.zeros(len(mesh.points))
    for i in range(15):
        np.add.at(load, nodes[:, i], np.einsum('f,q,q,fqk->fk', panels.areas, weights, values[:, i], modes))
        np.add.at(patch, nodes[:, i], np.einsum('f,q,q,fq->f', panels.areas, weights, values[:, i], flux))
    return load, patch


def diagnostics(matrix):
    return dict(reciprocity_relative=float(np.linalg.norm(matrix-matrix.T)/np.linalg.norm(matrix)),
                minimum_symmetric_eigenvalue=float(np.linalg.eigvalsh((matrix+matrix.T)/2).min()))


def run():
    domain = make_domain(12, False, 2)
    points, _, _ = geometry(domain, domain.rest.position, domain.rest.orientation,
                             np.zeros(len(domain.free_nodes)))
    current = points[domain.local_mapping]
    faces = domain.contact.faces
    markers = np.full(len(faces), 3)
    markers[domain.contact.body_faces] = 1
    markers[domain.contact.free_faces] = 2
    volume_points, tets, volume_faces, volume_markers, preflight = tetrahedralize(current, faces, markers)
    mesh = p4_tetra.make_mesh(volume_points, tets)
    if len(mesh.points) > 100000:
        raise RuntimeError(f'mesh exceeded the 100000-P4-node preflight cap: {len(mesh.points)}')
    operator = p4_tetra.stiffness(mesh)
    _, patch_load = loads(mesh, volume_faces, domain.rest.position)
    polynomial, _ = harmonic_quartic(mesh.points, domain.rest.position)
    patch_error = float(np.linalg.norm(operator@polynomial-patch_load)/np.linalg.norm(patch_load))
    if patch_error > 1e-9:
        raise RuntimeError(f'independent harmonic-quartic patch failed: {patch_error:.12g}')
    load, _ = loads(mesh, volume_faces[volume_markers == 1], domain.rest.position)
    fixed = np.unique(p4_tetra.face_nodes(mesh, volume_faces[volume_markers == 2]))
    unknown = np.setdiff1d(np.arange(len(mesh.points)), fixed)
    factor = splu(operator[unknown][:, unknown].tocsc())
    potential = np.zeros_like(load)
    potential[unknown] = factor.solve(load[unknown])
    actual = 1000*load.T@potential
    residual = float(np.linalg.norm((operator@potential-load)[unknown])/np.linalg.norm(load[unknown]))
    if residual > 1e-9:
        raise RuntimeError(f'current P4 system residual failed: {residual:.12g}')
    # Independently use the unchanged P1 BEM on exactly the original boundary.
    panels = TrianglePanels(current[faces], 1)
    bem_energy, bem_diagnostics = assemble_energy(faces, panels, 3)
    from types import SimpleNamespace
    bem_load = body_load(SimpleNamespace(vertices=current, faces=faces), panels,
                          domain.contact.body_faces, domain.rest.position)
    bem_fixed = np.unique(faces[domain.contact.free_faces])
    bem_unknown = np.setdiff1d(np.arange(len(current)), bem_fixed)
    bem_potential = np.zeros_like(bem_load)
    bem_potential[bem_unknown] = np.linalg.solve(bem_energy[np.ix_(bem_unknown, bem_unknown)], bem_load[bem_unknown])
    control = 1000*bem_load.T@bem_potential
    # Translation/rotation units are compared with the declared canoe scale.
    d = np.array([1., 1., 1., 3.5, 3.5, 3.5])
    normalized = (actual-control)/(1000*3.5**3*d[:, None]*d[None, :])
    report = dict(qualified=False, case='current canonical canoe, confined local-pocket inertia',
        wet_subdivisions=2, n=12, density=1000., rotation_center=domain.rest.position.tolist(),
        depth=1., tetrahedra=len(tets), p4_nodes=len(mesh.points), p4_unknowns=len(unknown),
        boundary_faces=len(faces), wetted_faces=int(np.sum(markers == 1)),
        current_wet_triangles_sha256=sha256(np.asarray(current[faces[markers == 1]][:, ::-1], dtype='<f8').tobytes()).hexdigest(),
        preflight=preflight, independent_harmonic_quartic_relative_patch_error=patch_error,
        relative_system_residual=residual, p4_inertia=actual.tolist(), p1_bem_inertia=control.tolist(),
        p4_diagnostics=diagnostics(actual), p1_bem_diagnostics=diagnostics(control),
        diagonal_relative_difference=(np.diag(actual-control)/np.diag(control)).tolist(),
        complete_nondimensional_difference=normalized.tolist(), bem_assembly_diagnostics=bem_diagnostics,
        scope='Single static confined-domain comparison, not a Capytaine qualification. '
              'Side walls and depth1 bed are the current local-pocket boundaries. '
              'No nonlinear moving-mesh, field, wave, domain-convergence or performance qualification.')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-current-canoe.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        report = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
