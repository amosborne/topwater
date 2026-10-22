"""Bounded same-mesh affine-rest K/B/M bridge for the actual canoe quarter.

Read an existing coarse physical mesh; no owner, remesher, fluid factor or
linear solve. Actual nonlinear curved helpers and frequency helpers act on
the SAME oriented mesh/CG. Zero stretching is essential: this does not bridge
the closed tank to a causal outgoing boundary or qualify nonlinear/Capytaine
accuracy. Run explicitly after source review.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from threadpoolctl import threadpool_limits

from .p4_tetra import make_mesh, face_nodes
from .pocket import boundary
from .p4_geometry import stiffness
from .p4_surface import integrals as surface_integrals
from .p4_canoe_field import body_load
from .p4_canoe import loads, harmonic_quartic
from .quarter_canoe import POCKET
from .quarter_canoe_control import scalar_operator
from .p4_first_variation import _hash, _progress, _revision


def _difference(actual, expected):
    delta = actual-expected
    values = delta.data if hasattr(delta, 'tocsr') else np.asarray(delta)
    reference = expected.data if hasattr(expected, 'tocsr') else np.asarray(expected)
    scale = float(np.linalg.norm(reference))
    return dict(maximum_absolute=float(np.max(abs(values))) if values.size else 0.,
        frobenius_error=float(np.linalg.norm(values)), reference_frobenius=scale,
        relative_frobenius=float(np.linalg.norm(values)/scale) if scale else None)


def audit(mesh_path):
    """Actual saved 69-tet/1049-node control only; never build a large case."""
    mesh_path = Path(mesh_path)
    metadata_path = mesh_path.with_suffix('.json')
    metadata = json.loads(metadata_path.read_text())
    with np.load(mesh_path) as saved:
        points, tets, body, center = [saved[name].copy() for name in
            ('points', 'tets', 'body_faces', 'center')]
    if (tets.ndim != 2 or tets.shape[1] != 4 or not 0 < len(tets) <= 69
            or points.ndim != 2 or points.shape[1] != 3 or len(points) > 1049
            or body.ndim != 2 or body.shape[1] != 3 or not len(body)
            or center.shape != (3,) or not all(np.all(np.isfinite(v)) for v in (points, center))):
        raise ValueError('bounded saved coarse actual canoe mesh required (at most69 tets)')
    mesh = make_mesh(points, tets)
    if (len(mesh.points) > 1049 or metadata['quarter_nodes'] != len(mesh.points)
            or metadata['quarter_tetrahedra'] != len(tets)
            or not np.array_equal(metadata['reference']['rotation_center'], center)):
        raise ValueError('coarse saved mesh differs from its physical control metadata')
    faces = boundary(points, tets)
    oriented = {tuple(np.roll(face, shift)) for face in faces for shift in range(3)}
    if any(tuple(face) not in oriented for face in body):
        raise ValueError('body faces must share the actual fluid-outward boundary winding')
    top = faces[np.all(points[faces, 2] == 0, axis=1)]
    if not len(top):
        raise ValueError('actual flat upward free-surface faces required')
    _progress('bounded saved canoe quarter reconstructed; no mesher/owner')

    # Identity vertical geometry coefficients at every local P4 height mode;
    # mass comes directly from the physical nonlinear surface primitive.
    top_nodes = face_nodes(mesh, top)
    local_mass = []
    for nodes in top_nodes:
        term = surface_integrals(mesh.points[nodes],
            (np.zeros((15, 15)), np.zeros((15, 15)), np.eye(15)), np.zeros((15, 3, 6)))
        local_mass.append(term['mass'])
    mass = coo_matrix((np.asarray(local_mass).ravel(),
        (np.repeat(top_nodes, 15, axis=1).ravel(), np.tile(top_nodes, (1, 15)).ravel())),
        shape=(len(mesh.points),)*2).tocsr()
    nonlinear_k = stiffness(mesh, mesh.points)
    nonlinear_b = body_load(mesh, body, mesh.points, center)
    frequency_b, _ = loads(mesh, body, center)
    _progress('actual nonlinear K/B/M and frequency six-column B assembled')
    frequencies = (2., 4.)
    operators = [scalar_operator(mesh, omega, POCKET+.5, [.5, .5], 0., order=12,
                                max_nodes=1049) for omega in frequencies]
    # A(omega)=K-omega²M/g. Two frequencies isolate K and M separately.
    frequency_mass = (operators[0]-operators[1])*9.81/(frequencies[1]**2-frequencies[0]**2)
    frequency_k = operators[0]+frequencies[0]**2/9.81*frequency_mass
    combined = [_difference(actual, nonlinear_k-omega**2/9.81*mass)
                for omega, actual in zip(frequencies, operators)]
    _progress('two zero-stretch operators assembled; no factor/solve')

    # This existing independent physical harmonic polynomial has nonzero
    # derivatives in all world axes. Full boundary flux includes body, planes,
    # surface and walls; body-only flux would not be the PDE patch identity.
    polynomial, _ = harmonic_quartic(mesh.points, center)
    _, patch = loads(mesh, faces, center)
    patch_nonlinear = _difference(nonlinear_k@polynomial, patch)
    patch_frequency = _difference(frequency_k@polynomial, patch)
    k_error, m_error, b_error = (_difference(frequency_k, nonlinear_k),
        _difference(frequency_mass, mass), _difference(frequency_b, nonlinear_b))
    b_columns = [_difference(frequency_b[:, i], nonlinear_b[:, i]) for i in range(6)]
    # Existing polynomial weak-form controls use1e-9; this is an algebraic
    # affine/quadrature identity threshold, not a hydrodynamic accuracy gate.
    tolerance = 1e-9
    checks = [k_error, m_error, b_error, patch_nonlinear, patch_frequency, *combined, *b_columns]
    passed = all(value['relative_frobenius'] is not None
                 and np.isfinite(value['relative_frobenius'])
                 and value['relative_frobenius'] <= tolerance for value in checks)
    sources = ['p4_rest_operator_bridge.py', 'p4_tetra.py', 'pocket.py', 'p4_geometry.py',
        'p4_surface.py', 'p4_canoe_field.py', 'p4_canoe.py', 'quarter_canoe.py',
        'quarter_canoe_control.py', 'pml_corner.py', 'pml_strip.py']
    report = dict(qualified=False, diagnostic_identity_passed=passed,
        case='actual coarse canoe same-mesh affine-rest nonlinear/frequency operator bridge',
        revision=_revision(), source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                            for name in sources},
        input_mesh=str(mesh_path), input_mesh_sha256=sha256(mesh_path.read_bytes()).hexdigest(),
        input_metadata_sha256=sha256(metadata_path.read_bytes()).hexdigest(),
        array_sha256=dict(reference_points=_hash(mesh.points, '<f8'), elements=_hash(mesh.elements, '<i8'),
            body_faces=_hash(body, '<i8'), free_faces=_hash(top, '<i8'), center=_hash(center, '<f8')),
        tetrahedra=len(tets), p4_nodes=len(mesh.points), physical_rotation_center=center.tolist(),
        frequencies=list(frequencies), stretch_strength=0., relative_identity_tolerance=tolerance,
        nonlinear_volume_quadrature_order=8, nonlinear_surface_quadrature_order=9,
        frequency_volume_surface_quadrature_order=12,
        stiffness_difference=k_error, mass_difference=m_error, body_load_difference=b_error,
        body_load_columns=b_columns, scalar_operator_differences=combined,
        nonlinear_harmonic_quartic_weak_flux_error=patch_nonlinear,
        frequency_harmonic_quartic_weak_flux_error=patch_frequency,
        scope='Same physical oriented affine P4 rest mesh and CG only. Calls the actual separate '
            'helpers without changing live operator code. Does not establish their identity on '
            'curved moving geometry, wet1/wet2 field convergence, nonlinear reaction/transients, '
            'causal/open-water boundary behavior, Capytaine or throughput qualification.')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, default=Path('artifacts/quarter-canoe-parity-control.npz'))
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-rest-operator-bridge.json'))
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('choose a new operator-bridge report; existing evidence is preserved')
    with threadpool_limits(1):
        report = audit(args.mesh)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
    if not report['diagnostic_identity_passed']:
        raise SystemExit(1)
