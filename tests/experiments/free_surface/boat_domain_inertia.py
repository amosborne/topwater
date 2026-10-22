"""Separate finite-tank bias from body discretization at flat startup.

Zero free-surface acceleration potential is the instantaneous/infinite-frequency
limit, not a finite-frequency or nonlinear open-water correlation.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from topwater import HullMesh
from topwater.frames import rotation_matrix
from topwater.inertia import instantaneous_added_mass
from topwater.rankine import TrianglePanels
from tests.reference_boat import LENGTH, WATER_DENSITY, equilibrium, reference_body
from .boat_world_boundary import construct
from .boat_world_field import refine_boundary
from .canoe import body_load
from .p4_fluid_ordering import fresh


def mirrored_body(triangles):
    """Exact reflected closed body and odd-potential basis, no coordinate weld."""
    triangles = np.asarray(triangles, dtype=float)
    if np.any(triangles[:, :, 2] > 0):
        raise ValueError('body-outward wetted triangles must lie below z=0')
    full = np.concatenate((triangles, triangles[:, ::-1]*[1, 1, -1]))
    points, indices = np.unique(full.reshape(-1, 3), axis=0, return_inverse=True)
    mesh = HullMesh(points, indices.reshape(-1, 3))
    lower = np.flatnonzero(points[:, 2] < 0)
    lookup = {tuple(point): index for index, point in enumerate(points)}
    upper = np.array([lookup[tuple(points[index]*[1, 1, -1])] for index in lower])
    odd = np.zeros((len(points), len(lower)))
    odd[lower, np.arange(len(lower))] = 1
    odd[upper, np.arange(len(lower))] = -1
    return mesh, odd


def run(output, quadrature_order=3, capytaine=False):
    from .cuda_operator import prepare_energy, cp
    output = fresh(output)
    if quadrature_order not in (3, 6, 12):
        raise ValueError('this bounded diagnostic supports quadrature orders 3, 6 and 12')
    output.parent.mkdir(parents=True, exist_ok=True)
    report = dict(complete=False, qualified=False, levels=[],
        scope='Same canonical flat-startup wetted body: finite-tank P1 versus unbounded mirrored P1 and production halfspace P0. Not wave or trajectory qualification.',
        density_kg_per_m3=WATER_DENSITY, length_scale_m=LENGTH,
        quadrature_order=quadrature_order,
        capytaine_reference_requested=capytaine,
        source_components_sha256={str(path): sha256(path.read_bytes()).hexdigest() for path in
            (Path('src/topwater/inertia.py'), Path('src/topwater/bem.py'),
             Path('src/topwater/boundary.py'), Path('tests/capytaine_reference.py'))},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_domain_inertia.py', 'boat_world_boundary.py', 'boat_world_field.py',
                         'canoe.py', 'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh')})
    saved = {}
    started = perf_counter()

    def budget():
        if perf_counter()-started > 120:
            raise RuntimeError('120 s diagnostic budget exhausted between operations')

    def describe(matrix):
        scale = np.r_[np.ones(3), np.full(3, LENGTH)]
        normalized = matrix/scale[:, None]/scale[None, :]
        return dict(added_mass_world=matrix.tolist(),
            normalized_reciprocity_defect=float(np.linalg.norm(normalized-normalized.T)/np.linalg.norm(normalized)),
            normalized_energy_eigenvalues=np.linalg.eigvalsh((normalized+normalized.T)/2).tolist())

    try:
        cp.cuda.get_current_stream().synchronize()
        initial_path = output.with_name(output.stem+'-flat-geometry.json')
        _, arrays, metadata = construct(initial_path, parent_eta=np.zeros(15), parent_psi=np.zeros(15))
        report['flat_geometry'] = dict(path=str(initial_path),
            sha256=sha256(initial_path.with_suffix('.npz').read_bytes()).hexdigest(),
            fixture_sha256=metadata['fixture_sha256'])
        body, pose = reference_body(), equilibrium()[0]
        rotation = rotation_matrix(pose.orientation)
        rigid = np.zeros((6, 6))
        rigid[:3, :3] = body.mass*np.eye(3)
        rigid[3:, 3:] = rotation@body.inertia_cg@rotation.T
        report.update(world_cg=pose.position.tolist(), orientation_wxyz=pose.orientation.tolist(),
            mass_kg=body.mass, rigid_mass_world=rigid.tolist())
        saved['world_cg'] = pose.position
        for level in range(3):
            budget()
            if level:
                arrays = refine_boundary(arrays, region='body')
            row = dict(level=level, complete=False, finite_panels=len(arrays['faces']),
                       body_panels=len(arrays['body_faces']))
            report['levels'].append(row)
            points, faces, fixed = arrays['points'], arrays['faces'], arrays['free_nodes']
            mesh = HullMesh(points, faces)
            panels = TrianglePanels.from_mesh(mesh)
            before = perf_counter()
            finite_assembly = prepare_energy(faces, panels, order=quadrature_order)
            finite_energy = finite_assembly.energy
            load = body_load(mesh, panels, arrays['body_faces'], pose.position)
            unknown = np.setdiff1d(np.arange(len(points)), fixed)
            response = np.zeros(load.shape)
            response[unknown] = cho_solve(cho_factor(finite_energy[np.ix_(unknown, unknown)]), load[unknown])
            finite_mass = WATER_DENSITY*load.T@response
            row['finite_p1'] = dict(**describe(finite_mass), seconds=perf_counter()-before,
                                    assembly=finite_assembly.diagnostics)
            del finite_assembly
            budget()
            # Reversing the finite-fluid face gives BODY-outward geometry.
            wet = points[faces[arrays['body_faces']]][:, ::-1]
            saved[f'wetted_triangles_level{level}'] = wet
            mirrored, odd = mirrored_body(wet)
            row['wetted_triangles_sha256'] = sha256(np.ascontiguousarray(wet).tobytes()).hexdigest()
            row['mirrored_panels'] = len(mirrored.faces)
            before = perf_counter()
            exterior_panels = TrianglePanels.from_mesh(mirrored, normal_sign=-1)
            exterior_assembly = prepare_energy(mirrored.faces, exterior_panels, order=quadrature_order)
            exterior_energy = exterior_assembly.energy
            lower_load = body_load(mirrored, exterior_panels, np.arange(len(wet)), pose.position)
            # Full odd energy is twice the physical lower-half energy.
            half_energy = .5*(odd.T@exterior_energy@odd)
            lower_response = odd@cho_solve(cho_factor(half_energy), odd.T@lower_load)
            deep_mass = WATER_DENSITY*lower_load.T@lower_response
            row['deep_p1'] = dict(**describe(deep_mass), seconds=perf_counter()-before,
                                  assembly=exterior_assembly.diagnostics)
            del exterior_assembly
            budget()
            before = perf_counter()
            production = instantaneous_added_mass(TrianglePanels(wet, -1), pose.position, density=WATER_DENSITY)
            row['deep_p0'] = dict(**describe(production.added_mass), seconds=perf_counter()-before,
                                  diagnostics=production.diagnostics)
            for label, matrix in [('finite_p1', finite_mass), ('deep_p1', deep_mass),
                                  ('deep_p0', production.added_mass)]:
                row[label]['unit_world_load_acceleration'] = np.linalg.solve(rigid+matrix, np.eye(6)).tolist()
                saved[f'{label}_level{level}'] = matrix
            row['finite_vs_deep_p1_diagonal_relative_difference'] = (
                np.diag(finite_mass-deep_mass)/np.diag(deep_mass)).tolist()
            row['deep_p0_vs_p1_diagonal_relative_difference'] = (
                np.diag(production.added_mass-deep_mass)/np.diag(deep_mass)).tolist()
            if capytaine:
                from tests.capytaine_reference import instantaneous_inertia
                before = perf_counter()
                reference, damping = instantaneous_inertia(wet, pose.position, body.mass,
                                                           density=WATER_DENSITY)
                if not np.isfinite(reference).all() or np.any(damping):
                    raise ValueError('invalid infinite-frequency Capytaine reference')
                saved[f'capytaine_level{level}'] = reference
                row['capytaine'] = dict(**describe(reference), seconds=perf_counter()-before,
                    unit_world_load_acceleration=np.linalg.solve(rigid+reference, np.eye(6)).tolist())
                row['capytaine_diagonal_relative_errors'] = {
                    label: (np.diag(matrix-reference)/np.diag(reference)).tolist()
                    for label, matrix in [('finite_p1', finite_mass), ('deep_p1', deep_mass),
                                           ('deep_p0', production.added_mass)]}
            row['complete'] = True
            np.savez(output.with_suffix('.npz'), **saved)
            output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
            print(json.dumps(dict(level=level, complete=True,
                body_panels=len(wet), finite_vs_deep_diagonal=row['finite_vs_deep_p1_diagonal_relative_difference'])), flush=True)
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['seconds'] = perf_counter()-started
        np.savez(output.with_suffix('.npz'), **saved)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--quadrature-order', type=int, choices=(3, 6, 12), default=3)
    parser.add_argument('--capytaine', action='store_true')
    result = run(**vars(parser.parse_args()))
    print(json.dumps(dict(complete=result['complete'], seconds=result['seconds'])))
