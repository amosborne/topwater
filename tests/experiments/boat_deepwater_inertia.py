"""Same-panel canonical-boat A-infinity component comparison, not a trajectory.

python -m tests.experiments.boat_deepwater_inertia --output artifacts/boat-deepwater-inertia.json
Requires the optional Capytaine 3.0.0 validation dependency. Fixed three levels
retain the exact same equilibrium physical hull; there is no finite tank bed.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
from time import perf_counter

import numpy as np
import scipy
from threadpoolctl import threadpool_info, threadpool_limits

from topwater.inertia import instantaneous_added_mass
from topwater.rankine import TrianglePanels
from tests.capytaine_reference import MODES, instantaneous_inertia
from tests.reference_boat import LENGTH, WATER_DENSITY, GRAVITY, reference_body, equilibrium
from tests.experiments.capytaine_mesh_cost import subdivide


def array_hash(values):
    a = np.ascontiguousarray(values)
    return dict(shape=list(a.shape), dtype=a.dtype.str,
                sha256=sha256(a.tobytes()).hexdigest())


def matrix_diagnostics(matrix):
    scale = np.array([1., 1., 1., LENGTH, LENGTH, LENGTH])
    normalized = matrix/scale[:, None]/scale[None, :]/(WATER_DENSITY*LENGTH**3)
    eigenvalues = np.linalg.eigvalsh((normalized+normalized.T)/2)
    diagonal = np.diag(matrix)
    coupling = (matrix/np.sqrt(diagonal[:, None]*diagonal[None, :])
                if np.all(diagonal > 0) else None)
    return dict(nondimensional=normalized.tolist(),
        normalized_reciprocity_max=float(np.abs(normalized-normalized.T).max()),
        reciprocity_relative_frobenius=float(np.linalg.norm(normalized-normalized.T)/np.linalg.norm(normalized)),
        symmetric_part_eigenvalues=eigenvalues.tolist(), positive_energy=bool(eigenvalues.min() > 0),
        coupling_coefficients=None if coupling is None else coupling.tolist())


def run(output):
    import capytaine as cpt
    output = Path(output).with_suffix('.json')
    if output.exists() or output.with_suffix('.npz').exists():
        raise FileExistsError('preserve prior evidence; choose a new output stem')
    output.parent.mkdir(parents=True, exist_ok=True)
    source_paths = [Path(__file__), Path('tests/reference_boat.py'), Path('tests/capytaine_reference.py'),
        Path('tests/experiments/capytaine_mesh_cost.py'), Path('src/topwater/inertia.py'),
        Path('src/topwater/bem.py'), Path('src/topwater/boundary.py'), Path('src/topwater/rankine.py'),
        Path('src/topwater/_rankine_kernel.py'), Path('src/topwater/_rankine_extended.py')]
    report = dict(complete=False, qualified=False,
        scope='Same-panel deep-water zero-acceleration-potential-plane A-infinity component only; not finite-frequency, spatially converged, finite-tank, nonlinear pressure or full-model qualification',
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, capytaine=cpt.__version__,
        density_kg_per_m3=WATER_DENSITY, gravity_m_per_s2=GRAVITY,
        water_depth='infinite', acceleration_potential_free_surface_value=0.,
        normalization_length_m=LENGTH, mode_order=list(MODES),
        matrix_units='translation/translation kg; translation/rotation kg*m; rotation/rotation kg*m^2',
        source_sha256={str(p): sha256(p.read_bytes()).hexdigest() for p in source_paths},
        levels=[])
    arrays, started = {}, perf_counter()

    def checkpoint():
        report['elapsed_component_seconds'] = perf_counter()-started
        report['threadpools'] = threadpool_info()
        np.savez(output.with_suffix('.npz'), **arrays)
        report['array_archive_sha256'] = sha256(output.with_suffix('.npz').read_bytes()).hexdigest()
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    try:
        body, (pose, hydro) = reference_body(), equilibrium()
        triangles = hydro.wetted_triangles.copy()
        if len(triangles) != 17:
            raise ValueError('canonical equilibrium clip changed; expected17 original wet triangles')
        report.update(mass_kg=body.mass, material_cg_m=body.cg.tolist(),
            world_cg_m=pose.position.tolist(), orientation_wxyz=pose.orientation.tolist(),
            body_inertia_cg_kg_m2=body.inertia_cg.tolist(),
            fixture_vertices=array_hash(body.mesh.vertices), fixture_faces=array_hash(body.mesh.faces),
            physical_face_indices=array_hash(body.mesh.physical_faces),
            hydrostatic_force_N=hydro.force.tolist(), hydrostatic_moment_Nm=hydro.moment.tolist())
        arrays.update(world_cg=pose.position, orientation_wxyz=pose.orientation)
        for level in range(3):
            if level:
                triangles = subdivide(triangles)
            row = dict(level=level, complete=False, panels=len(triangles), triangle_identity=array_hash(triangles))
            report['levels'].append(row)
            arrays[f'triangles_{level}'] = triangles.copy()
            checkpoint()
            before = perf_counter()
            with threadpool_limits(1):
                actual = instantaneous_added_mass(TrianglePanels(triangles, -1), pose.position,
                                                  density=WATER_DENSITY)
            row['topwater_component_seconds'] = perf_counter()-before
            before = perf_counter()
            reference, damping = instantaneous_inertia(triangles, pose.position, body.mass,
                density=WATER_DENSITY, gravity=GRAVITY)
            row['capytaine_component_seconds'] = perf_counter()-before
            arrays[f'topwater_{level}'], arrays[f'capytaine_{level}'] = actual.added_mass, reference
            arrays[f'capytaine_damping_{level}'] = damping
            if not np.isfinite(actual.added_mass).all() or not np.isfinite(reference).all():
                raise FloatingPointError('nonfinite raw added mass')
            if not np.array_equal(damping, np.zeros((6, 6))):
                raise ValueError('A-infinity reference returned nonzero radiation damping')
            scale = np.array([1., 1., 1., LENGTH, LENGTH, LENGTH])
            delta = actual.added_mass-reference
            nd_delta = delta/scale[:, None]/scale[None, :]/(WATER_DENSITY*LENGTH**3)
            nd_reference = reference/scale[:, None]/scale[None, :]/(WATER_DENSITY*LENGTH**3)
            floor = max(1e-10, .005*float(np.abs(nd_reference).max()))
            ratio = np.abs(nd_delta)/(.05*np.abs(nd_reference)+floor)
            off = ~np.eye(6, dtype=bool)
            row.update(complete=True, added_mass=actual.added_mass.tolist(),
                capytaine_added_mass=reference.tolist(), raw_matrix_difference=delta.tolist(),
                diagonal_relative_errors=(np.diag(delta)/np.diag(reference)).tolist(),
                offdiagonal_raw_differences=delta[off].tolist(),
                offdiagonal_normalized_error_frobenius=float(np.linalg.norm(nd_delta[off])),
                normalized_error_frobenius=float(np.linalg.norm(nd_delta)),
                topwater=matrix_diagnostics(actual.added_mass), capytaine=matrix_diagnostics(reference),
                informational_existing_family_gate=dict(absolute_floor=floor,
                    tolerance_ratio_all36=ratio.tolist(), within_tolerance=bool(np.all(ratio <= 1)),
                    caveat='Global floor does not establish small-coupling or spatial convergence accuracy'),
                topwater_diagnostics=actual.diagnostics)
            if level:
                row['previous_level_diagonal_relative_change'] = {
                    name: ((np.diag(arrays[f'{name}_{level}'])-np.diag(arrays[f'{name}_{level-1}']))
                           /np.diag(arrays[f'{name}_{level}'])).tolist()
                    for name in ('topwater', 'capytaine')}
            checkpoint()
        report['complete'] = True
    except BaseException as error:
        report['error'] = dict(type=type(error).__name__, message=str(error))
        checkpoint()
        raise
    checkpoint()
    print(json.dumps(dict(output=str(output), complete=True, qualified=False)))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output)
