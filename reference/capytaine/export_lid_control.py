"""Isolated native68-hull internal-lid radiation control;24 calls, no repair."""
import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from importlib.metadata import distribution, version
from pathlib import Path

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'

import numpy as np
from scipy.io import loadmat, savemat
from scipy.spatial import ConvexHull
from threadpoolctl import threadpool_limits

MODES = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--controls', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    controls = read(args.controls)
    if controls['status'] != 'ready_for_coordinator_execution_allocation':
        raise ValueError('Reviewed READY controls and separate allocation required')
    frozen = controls['input_sha256']
    for path, digest in frozen.items():
        assert sha(path) == digest, path
    args.output.mkdir(parents=True, exist_ok=False)
    snapshot = args.output/'source'
    snapshot.mkdir()
    for i, path in enumerate(dict.fromkeys([str(Path(__file__)), str(args.controls), *frozen])):
        shutil.copyfile(path, snapshot/f'{i:02d}-{Path(path).name}')
    started = time.perf_counter()
    report = dict(complete=False, expected_calls=24, completed_calls=0, cpu_threads=1,
        command=[sys.executable, *sys.argv], input_sha256=frozen, source_sha256=sha(__file__),
        controls_sha256=sha(args.controls), scope='Isolated native internal-lid mechanism control; no replacement coefficients or physical qualification')

    def save():
        (args.output/'manifest.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    save()
    try:
        import capytaine as cpt
        from capytaine.green_functions.delhommeau import Delhommeau
        if version('capytaine') != '3.0.0':
            raise ValueError('Reviewed Capytaine3.0.0 API required')
        original = read(controls['original_manifest'])
        package = distribution('capytaine')
        actual = {str(p): sha(package.locate_file(p)) for p in package.files
                  if str(p).endswith(('.py', '.pyd', '.so', '.dll'))}
        assert actual == original['dependency_sha256']['capytaine']
        source = loadmat(controls['no_lid_mat'])
        native = loadmat(controls['native_mat'])
        triangles = source['triangles']
        assert triangles.shape == (68, 3, 3)
        assert all(np.array_equal(source[key], native[key]) for key in ('M', 'K', 'cg'))
        vertices, ids = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
        mesh = cpt.Mesh(vertices=vertices, faces=ids.reshape(-1, 3), auto_clean=False, auto_check=False)
        hull_vertices, hull_faces = mesh.vertices.copy(), mesh.faces.copy()
        hull_area = float(mesh.faces_areas.sum())
        # Native API generates interior horizontal panels, not a physical cap.
        lid = mesh.generate_lid(z=0.0, faces_max_radius=0.25)
        if not (0 < lid.nb_faces <= 512):
            raise ValueError('Native lid empty or exceeds bounded512-panel control')
        body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid,
            dofs=cpt.rigid_body_dofs(rotation_center=native['cg'].ravel()),
            center_of_mass=native['cg'].ravel(), mass=native['M'][0, 0], name='same68_hull_native_lid')
        assert list(body.dofs) == MODES
        assert np.array_equal(body.mesh.vertices, hull_vertices) and np.array_equal(body.mesh.faces, hull_faces)
        assert body.mesh.nb_faces == 68 and int(body.hull_mask.sum()) == 68
        assert body.mesh.faces_areas.sum() == hull_area
        assert np.allclose(body.lid_mesh.vertices[:, 2], 0, rtol=0, atol=1e-13)
        np.testing.assert_allclose(body.lid_mesh.faces_normals,
                                   np.tile([0., 0., -1.], (lid.nb_faces, 1)), rtol=0, atol=1e-13)
        rho, gravity = native['rho'].item(), native['g'].item()
        waterline = vertices[abs(vertices[:, 2]) < 1e-12, :2]
        polygon = ConvexHull(waterline)
        exact_waterplane = native['K'][2, 2]/rho/gravity
        # Convex test is a validity check on this boat, never an outside-hull fill.
        np.testing.assert_allclose(polygon.volume, exact_waterplane, rtol=2e-12, atol=2e-12)
        assert np.max(body.lid_mesh.vertices[:, :2]@polygon.equations[:, :2].T+polygon.equations[:, 2]) <= 2e-12
        lid_area = float(body.lid_mesh.faces_areas.sum())
        assert 0 < lid_area <= exact_waterplane+2e-12
        report['lid_invariants'] = dict(physical_hull_faces=68, hull_exact=True, physical_hull_area_m2=hull_area,
            waterplane_area_m2=float(exact_waterplane), lid_faces=lid.nb_faces, lid_area_m2=lid_area,
            waterplane_coverage_fraction=lid_area/exact_waterplane, horizontal=True, normals_down=True,
            vertices_inside_physical_waterplane=True, no_lid_dofs=True, metadata_M_K_CG_unchanged=True)
        np.savez(args.output/'lid-geometry.npz', physical_triangles=triangles,
            lid_vertices=body.lid_mesh.vertices, lid_faces=body.lid_mesh.faces, hull_mask=body.hull_mask)
        omega = np.array([1., 2., 12.])
        A, B, Ainf = np.empty((3, 6, 6)), np.empty((3, 6, 6)), np.empty((6, 6))
        with threadpool_limits(1):
            solver = cpt.BEMSolver(green_function=Delhommeau(tabulation_cache_dir=None))
            for i, w in enumerate([*omega, np.inf]):
                for column, mode in enumerate(MODES):
                    if time.perf_counter()-started > 300:
                        raise TimeoutError('Declared300s elapsed budget exceeded')
                    problem = cpt.RadiationProblem(body=body, omega=w, radiating_dof=mode,
                        free_surface=0.0, water_depth=np.inf, rho=rho, g=gravity)
                    # Native symbolic-infinity wrapper stores finite coefficient
                    # values; lid velocity stays exactly zero for every mode.
                    bc = problem.boundary_condition
                    bc_values = bc.value if hasattr(bc, 'symbol') else bc
                    assert np.all(bc_values[~body.hull_mask] == 0)
                    result = solver.solve(problem, keep_details=False)
                    report['completed_calls'] += 1
                    if np.isfinite(w):
                        A[i, :, column] = [result.added_mass[row] for row in MODES]
                        B[i, :, column] = [result.radiation_damping[row] for row in MODES]
                    else:
                        Ainf[:, column] = [result.added_mass[row] for row in MODES]
                save()
                print(f'Native lid control completed omega={w}', flush=True)
        assert report['completed_calls'] == 24 and all(np.isfinite(v).all() for v in (A, B, Ainf))
        output = args.output/'radiation-68-with-lid.mat'
        savemat(output, dict(omega_rad_s=omega[:, None], added_mass=A, radiation_damping=B, Ainf=Ainf,
            M=native['M'], K=native['K'], cg=native['cg'], cb=native['cb'], volume=native['volume'],
            rho=rho, g=gravity, dof_names=MODES))
        written = loadmat(output)
        for key, value in dict(omega_rad_s=omega[:, None], added_mass=A, radiation_damping=B, Ainf=Ainf,
                              M=native['M'], K=native['K'], cg=native['cg'], cb=native['cb'], volume=native['volume']).items():
            assert np.array_equal(written[key], value), key
        indices = [np.flatnonzero(source['omega_rad_s'].ravel() == w).item() for w in omega]
        diagnostic = dict(omega_rad_s=omega.tolist(),
            A_lid_minus_no_lid=(A-source['added_mass'][indices]).tolist(),
            B_lid_minus_no_lid=(B-source['radiation_damping'][indices]).tolist(),
            Ainf_lid_minus_no_lid=(Ainf-source['Ainf']).tolist(),
            interpretation='Raw isolated mechanism comparison only; no positivity repair, fit, cutofftail or kernel qualification')
        (args.output/'diagnostic.json').write_text(json.dumps(diagnostic, indent=2)+'\n')
        assert all(sha(path) == digest for path, digest in frozen.items())
        report.update(complete=True, elapsed_seconds=time.perf_counter()-started,
            output_sha256={p.name: sha(p) for p in args.output.iterdir() if p.is_file() and p.name != 'manifest.json'})
    except BaseException as error:
        report.update(error=repr(error), elapsed_seconds=time.perf_counter()-started)
    save()
    print(json.dumps({k: report.get(k) for k in ('complete', 'completed_calls', 'elapsed_seconds', 'error')}))
    return 0 if report['complete'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
