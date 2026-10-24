"""Excitation-only periodic table for the mild buoy incoming train.

SOURCE ONLY until separate review, READY controls and execution allocation.
No radiation solve, interpolation, trajectory, or mutation of original H5/CIC.
"""
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
        raise ValueError('Separate source/interface review and numerical allocation required')
    for path, digest in controls['input_sha256'].items():
        if sha(path) != digest:
            raise ValueError('Frozen input differs: '+path)
    args.output.mkdir(parents=True, exist_ok=False)
    snapshot = args.output/'source'
    snapshot.mkdir()
    for i, path in enumerate(dict.fromkeys([str(Path(__file__)), str(args.controls), *controls['input_sha256']])):
        shutil.copyfile(path, snapshot/f'{i:02d}-{Path(path).name}')
    started = time.perf_counter()
    report = dict(complete=False, completed_calls=0, expected_calls=816,
        new_radiation_calls=0, cpu_threads=1, command=[sys.executable, *sys.argv],
        source_sha256=sha(__file__), controls_sha256=sha(args.controls),
        input_sha256=controls['input_sha256'], versions={p: version(p) for p in
        ('capytaine', 'numpy', 'scipy', 'threadpoolctl')},
        scope='Excitation coverage only; interpolation and mesh accuracy unqualified')

    def save():
        (args.output/'manifest.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    save()
    try:
        import capytaine as cpt
        from capytaine.bem.airy_waves import froude_krylov_force
        from capytaine.green_functions.delhommeau import Delhommeau
        if version('capytaine') != '3.0.0':
            raise ValueError('Pinned Capytaine3.0.0 required')
        original = read(controls['original_manifest'])
        package = distribution('capytaine')
        dependencies = {str(p): sha(package.locate_file(p)) for p in package.files
                        if str(p).endswith(('.py', '.pyd', '.so', '.dll'))}
        if dependencies != original['dependency_sha256']['capytaine']:
            raise ValueError('Installed Capytaine differs from original producer')
        report['capytaine_dependency_sha256'] = dependencies
        case, native = read(controls['canonical_case']), loadmat(controls['native_mat'])
        with np.load(controls['native_npz']) as archive:
            triangles = archive['triangles'].copy()
        assert triangles.shape == (17, 3, 3)
        assert np.array_equal(native['M'], case['dry_mass_matrix_world'])
        assert np.array_equal(native['cg'].ravel(), case['world_cg_m'])
        vertices, indices = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
        mesh = cpt.Mesh(vertices=vertices, faces=indices.reshape(-1, 3), auto_clean=False, auto_check=False)
        body = cpt.FloatingBody(mesh=mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=native['cg'].ravel()),
            center_of_mass=native['cg'].ravel(), mass=case['body_mass_kg'], name='canonical_boat')
        assert list(body.dofs) == MODES and mesh.nb_faces == 17
        assert np.array_equal(mesh.vertices[mesh.faces[:, :3]], triangles)
        rho, gravity = native['rho'].item(), native['g'].item()
        omega = np.arange(1, 35, dtype=float)*.1
        theta = np.arange(24, dtype=float)*np.pi/12
        fk = np.empty((34, 24, 6), complex)
        diffraction = np.empty_like(fk)
        with threadpool_limits(1):
            solver = cpt.BEMSolver(green_function=Delhommeau(tabulation_cache_dir=None))
            for iw, w in enumerate(omega):
                for ih, heading in enumerate(theta):
                    if time.perf_counter()-started > controls['maximum_elapsed_seconds']:
                        raise TimeoutError('Excitation export wall budget reached')
                    problem = cpt.DiffractionProblem(body=body, omega=w, wave_direction=heading,
                        rho=rho, g=gravity, water_depth=np.inf)
                    result = solver.solve(problem, keep_details=False)
                    report['completed_calls'] += 1
                    native_fk = froude_krylov_force(problem)
                    fk[iw, ih] = [native_fk[mode] for mode in MODES]
                    diffraction[iw, ih] = [result.forces[mode] for mode in MODES]
                save()
                print(f"Excitation calls {report['completed_calls']}/816", flush=True)
        H = fk+diffraction
        assert all(np.isfinite(a).all() for a in (fk, diffraction, H))
        # Independent centroid incident-pressure integration at EVERY actual node.
        area_normal = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])/2
        centers = triangles.mean(axis=1)
        pressure = rho*gravity*np.exp((omega[:, None, None]**2/gravity)*(
            centers[None, None, :, 2]+1j*(np.cos(theta)[None, :, None]*centers[None, None, :, 0]
            +np.sin(theta)[None, :, None]*centers[None, None, :, 1])))
        panel_force = -pressure[:, :, :, None]*area_normal[None, None, :, :]
        expected_fk = np.concatenate((panel_force.sum(axis=2), np.cross(
            centers[None, None, :, :]-native['cg'].ravel(), panel_force).sum(axis=2)), axis=2)
        np.testing.assert_allclose(fk, expected_fk, rtol=3e-12, atol=1e-9)
        # Append an EXACT copy of heading0 at2pi; no extra BEM solve or extrapolation.
        seam = lambda a: np.concatenate((a, a[:, :1, :]), axis=1)
        bank = read(controls['bank_case'])
        table = dict(omega_rad_s=omega[:, None], heading_native_rad=np.r_[theta, 2*np.pi][:, None],
            fk_complex=seam(fk), diffraction_complex=seam(diffraction), excitation_complex=seam(H),
            dof_names=MODES, cg=native['cg'], M=native['M'], K=native['K'],
            rho=rho, g=gravity, triangles=triangles,
            yaw_global_rad=bank['mean_box_yaw_global_rad'],
            global_from_native_R=np.asarray(bank['global_from_boat_rotation']),
            global_from_native_t=np.asarray(bank['global_from_boat_translation_m']),
            bank_case_sha256=sha(controls['bank_case']),
            bank_components_sha256=sha(controls['bank_components']),
            convention='exp(-i omega t); dimensional per unit elevation; no conjugation/CG phase',
            layout='omega x periodic_heading x influenced_dof',
            interpolation='Consumer real/imag bilinear, modulo2pi, omega bounded [.1,3.4]; accuracy unqualified')
        savemat(args.output/'excitation.mat', table)
        np.savez(args.output/'excitation.npz', **table,
            heading_rad=table['heading_native_rad'], H_complex=table['excitation_complex'],
            FK_complex=table['fk_complex'])
        actual = loadmat(args.output/'excitation.mat')
        for key in ('omega_rad_s', 'heading_native_rad', 'fk_complex', 'diffraction_complex',
                    'excitation_complex', 'cg', 'M', 'K', 'triangles', 'global_from_native_R', 'global_from_native_t'):
            np.testing.assert_array_equal(actual[key].ravel(), np.asarray(table[key]).ravel())
        assert all(np.array_equal(actual[k][:, 0], actual[k][:, -1]) for k in
                   ('fk_complex', 'diffraction_complex', 'excitation_complex'))
        assert np.array_equal(actual['excitation_complex'], actual['fk_complex']+actual['diffraction_complex'])
        with np.load(args.output/'excitation.npz') as archive:
            for key, expected in [('omega_rad_s', table['omega_rad_s']),
                                  ('heading_rad', table['heading_native_rad']),
                                  ('H_complex', table['excitation_complex']),
                                  ('FK_complex', table['fk_complex']),
                                  ('diffraction_complex', table['diffraction_complex'])]:
                np.testing.assert_array_equal(archive[key], expected)
        assert abs(np.real(1j*np.exp(-1j*np.pi/2))-1) < 1e-15
        for path, digest in controls['input_sha256'].items():
            assert sha(path) == digest, path
        report.update(complete=True, elapsed_s=time.perf_counter()-started,
            independent_fk_max_abs_error=float(np.max(abs(fk-expected_fk))),
            checks=dict(all_816_nodes_finite=True, full_six_dof=True, exact_metadata_readback=True,
                periodic_seam_exact=True, FK_plus_diff=True, quarter_cycle_native_sign=True,
                input_and_original_radiation_H5_unchanged=True),
            output_sha256={p.name: sha(p) for p in (args.output/'excitation.mat', args.output/'excitation.npz')})
        save()
    except BaseException as error:
        report.update(error=repr(error), elapsed_s=time.perf_counter()-started)
        save()
        raise


if __name__ == '__main__':
    main()
