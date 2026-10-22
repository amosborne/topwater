"""Actual uncompressed PML quarter mesh-only preflight; no fluid operator.

The n16/n18 names identify previously phase-tested maximum spacing envelopes,
not a claim that this graded, canoe-containing mesh passes a propagation test.
Only the fixed geometry mesher and topological P4 enrichment run here. No
whole-domain mirror, stiffness assembly, fluid factorization or solve occurs.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np
from threadpoolctl import threadpool_limits

from .quarter_canoe import POCKET, make_plc, make_quarter


WIDTH = 7.877580409572782
BUFFER = .5
DEPTH_LEVELS = np.array([-4., -2.25, -1., -.25, 0.])


def envelope_grid(envelope):
    """Retain all physical/PML planes; never enlarge the cell-size envelope."""
    if envelope not in (16, 18):
        raise ValueError('use the predeclared n16 or n18 spacing envelope')
    maximum = np.array([4*np.pi/envelope, 2*np.pi/envelope])
    axes = []
    for pocket, spacing in zip(POCKET, maximum):
        planes = [0., pocket, pocket+BUFFER, pocket+BUFFER+WIDTH]
        levels = [np.array([0.])]
        for a, b in zip(planes[:-1], planes[1:]):
            count = int(np.ceil((b-a)/spacing))
            levels.append(np.linspace(a, b, count+1)[1:])
        coordinate = np.concatenate(levels)
        if not all(plane in coordinate for plane in planes) or np.max(np.diff(coordinate)) > spacing*(1+8*np.finfo(float).eps):
            raise RuntimeError('grid lost an exact plane or exceeded its spacing envelope')
        axes.append(coordinate)
    return axes[0], axes[1], DEPTH_LEVELS.copy()


def reflected_body(triangles, parents, parent_reflections):
    """Four BODY-outward copies with original physical-facet IDs, no volume."""
    copies = []; reflected_parents = []
    for row, (sx, sy) in enumerate(((1, 1), (1, -1), (-1, 1), (-1, -1))):
        reflected = triangles*np.array([sx, sy, 1])
        if sx*sy < 0:
            reflected = reflected[:, ::-1]
        target = parent_reflections[row, parents]
        if np.any(target < 0):
            raise ValueError('missing reflected original hull parent')
        copies.append(reflected); reflected_parents.append(target)
    return np.concatenate(copies), np.concatenate(reflected_parents)


def digest(array):
    """Hash shape, explicit little-endian dtype and all values."""
    array = np.asarray(array)
    normalized = np.asarray(array, dtype='<f8' if array.dtype.kind == 'f' else '<i8')
    value = sha256(str(normalized.shape).encode()+normalized.dtype.str.encode())
    value.update(normalized.tobytes())
    return value.hexdigest()


def run(envelope, reference_path):
    from .canoe_wet_reference import pinned_reference, verify_coverage

    x, y, z = envelope_grid(envelope)
    plc = make_plc(x, y, z, subdivisions=1)
    mesh, faces, body, parents, mesher_report = make_quarter(plc)
    # Body faces returned by make_quarter wind out of the fluid. The coverage
    # verifier takes body-outward triangles, including odd reflection reversal.
    candidate, candidate_parents = reflected_body(mesh.points[body][:, ::-1], parents,
                                                 plc.parent_reflections)
    domain, reference, _, provenance = pinned_reference(reference_path)
    coverage = verify_coverage(candidate, candidate_parents, reference,
                               domain.contact.body_source_faces)
    sx = mesh.points[:, 0] == 0; sy = mesh.points[:, 1] == 0
    # Inclusion-exclusion counts the conceptual four-copy P4 mesh exactly,
    # without constructing it. Shared axes use topological node identification.
    conceptual_full = 4*len(mesh.points)-2*int(sx.sum())-2*int(sy.sum())+int((sx&sy).sum())
    free = mesh.points[:, 2] == 0
    parity_unknown = {}
    for name, parity in [('EE', (1, 1)), ('OE', (-1, 1)), ('EO', (1, -1)), ('OO', (-1, -1))]:
        fixed = (sx if parity[0] == -1 else np.zeros(len(sx), bool)) | (sy if parity[1] == -1 else np.zeros(len(sy), bool))
        parity_unknown[name] = int((~fixed).sum())
    report = dict(qualified=False, mesh_only=True, spacing_envelope=envelope,
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('quarter_canoe_preflight.py', 'quarter_canoe.py', 'canoe_wet_reference.py',
                         'p4_canoe.py', 'p4_tetra.py')},
        x=x.tolist(), y=y.tolist(), z=z.tolist(),
        maximum_horizontal_spacing=[float(np.max(np.diff(x))), float(np.max(np.diff(y)))],
        physical_pocket=POCKET.tolist(), pml_onset=(POCKET+BUFFER).tolist(),
        pml_width=[WIDTH, WIDTH], pml_strength=8., real_terminal=1., wet_subdivisions=1,
        quarter_linear_vertices=mesh.vertex_count, quarter_p4_nodes=len(mesh.points),
        quarter_tetrahedra=len(mesh.tets), quarter_body_triangles=len(body),
        quarter_free_p4_nodes=int(free.sum()), quarter_parity_unknowns=parity_unknown,
        conceptual_full_p4_nodes=conceptual_full, conceptual_full_tetrahedra=4*len(mesh.tets),
        quarter_symmetry_plane_p4_nodes=[int(sx.sum()), int(sy.sum())],
        quarter_symmetry_axis_p4_nodes=int((sx&sy).sum()),
        array_sha256=dict(linear_points=digest(mesh.points[:mesh.vertex_count]),
            linear_tets=digest(mesh.tets), body_faces=digest(body), body_parent_ids=digest(parents),
            p4_master_points=digest(mesh.points), p4_elements=digest(mesh.elements)),
        mesher_preflight=mesher_report, reference=provenance, physical_coverage=coverage,
        scope='Actual uncompressed canonical-canoe quarter mesh-only preflight. '
              'Conceptual mirrored counts do not allocate a full mesh. No fluid matrix, '
              'LU-memory bound, propagation, nonlinear or Capytaine accuracy qualification.')
    saved = dict(points=mesh.points[:mesh.vertex_count], tets=mesh.tets,
                 body_faces=body, body_parents=parents, center=plc.center)
    return report, saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--envelope', type=int, choices=(16, 18), default=16)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('artifacts/quarter-canoe-actual-preflight.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        report, saved = run(args.envelope, args.reference)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output.with_suffix('.npz'), **saved)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(output=str(args.output), quarter_p4_nodes=report['quarter_p4_nodes'],
        conceptual_full_p4_nodes=report['conceptual_full_p4_nodes'],
        quarter_tetrahedra=report['quarter_tetrahedra'])), flush=True)
