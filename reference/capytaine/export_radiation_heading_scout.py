"""One bounded diagnostic:66 radiation solves plus384 Jan05 offset90 pairs.

No automatic refinement cascade, fitting, wave edits or physical qualification.
Run only after independent review, frozen manifest and coordinator allocation.
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

from export_case import triangles_at

MODES = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def moments(triangles):
    a, b, c = triangles.transpose(1, 0, 2)
    cross = np.cross(b-a, c-a)
    area = np.linalg.norm(cross, axis=1)/2
    first = area[:, None]*(a+b+c)/3
    total = a+b+c
    second = area[:, None, None]*(np.einsum('ni,nj->nij', total, total)
        +sum(np.einsum('ni,nj->nij', v, v) for v in (a, b, c)))/12
    return cross/2, area, first, second


def surface_check(parent, children):
    # triangles_at concatenates each of the four child classes, each in parent order.
    groups = children.reshape(4, len(parent), 3, 3).transpose(1, 0, 2, 3)
    parent_moments = moments(parent)
    child_moments = moments(groups.reshape(-1, 3, 3))
    for expected, values in zip(parent_moments, child_moments):
        summed = values.reshape((len(parent), 4)+values.shape[1:]).sum(axis=1)
        np.testing.assert_allclose(summed, expected, rtol=2e-13, atol=2e-12)
    normal = parent_moments[0]/parent_moments[1][:, None]
    for p, group, n in zip(parent, groups, normal):
        np.testing.assert_allclose((group-p[0])@n, 0, atol=2e-13)
        # Child vertices must remain inside the same physical parent triangle.
        basis = np.column_stack((p[1]-p[0], p[2]-p[0]))
        barycentric = np.linalg.lstsq(basis, (group-p[0]).reshape(-1, 3).T, rcond=None)[0]
        assert np.min(barycentric) >= -2e-13 and np.max(barycentric.sum(axis=0)) <= 1+2e-13
        assert np.all(moments(group)[0]@n > 0)
    return dict(parent_faces=17, child_faces=68,
                checked='Per-parent oriented area, area, first/second moments, plane, winding and barycentric containment')


def check_pairs(path, heading, global_bank, triangles, mat, expected_hashes):
    values = loadmat(path)
    bank = heading['boat_frame_waves']
    expected = dict(component_index=np.arange(1, 385)[:, None],
        omega_rad_s=np.asarray(bank['angular_frequencies'])[:, None],
        heading_native_rad=np.asarray(bank['directions'])[:, None],
        amplitude_m=np.asarray(bank['amplitudes'])[:, None],
        phase_local_rad=np.asarray(bank['phases'])[:, None], cg=mat['cg'],
        rho=mat['rho'], g=mat['g'], yaw_global_rad=np.array([[heading['mean_box_yaw_global_rad']]]),
        global_from_native_R=np.asarray(heading['global_from_boat_rotation']),
        global_from_native_t=np.atleast_2d(heading['global_from_boat_translation_m']))
    for key, array in expected.items():
        assert np.array_equal(values[key], array), key
    assert np.char.strip(values['dof_names'].astype(str)).ravel().tolist() == MODES
    for key, digest in expected_hashes.items():
        assert values[key].item().strip() == digest, key
    for key in ('fk_complex', 'diffraction_complex', 'excitation_complex'):
        assert values[key].shape == (384, 6) and np.isfinite(values[key]).all()
    assert np.array_equal(values['excitation_complex'], values['fk_complex']+values['diffraction_complex'])
    assert abs(np.real(1j*np.exp(-1j*np.pi/2))-1) < 1e-15
    index = int(np.flatnonzero(values['amplitude_m'].ravel() > 0)[0])
    w, theta = values['omega_rad_s'][index, 0], values['heading_native_rad'][index, 0]
    gravity, rho = mat['g'].item(), mat['rho'].item()
    center = triangles.mean(axis=1)
    area_normal = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])/2
    k = w*w/gravity
    pressure = rho*gravity*np.exp(k*center[:, 2]+1j*k*(center[:, 0]*np.cos(theta)+center[:, 1]*np.sin(theta)))
    force = -pressure[:, None]*area_normal
    fk = np.r_[force.sum(axis=0), np.cross(center-mat['cg'].ravel(), force).sum(axis=0)]
    np.testing.assert_allclose(fk, values['fk_complex'][index], rtol=2e-12, atol=1e-9)
    times = np.array([0., .37, 3.1])
    w, phase, amplitude = (values[key].ravel() for key in ('omega_rad_s', 'phase_local_rad', 'amplitude_m'))
    psi = phase[:, None]-w[:, None]*times
    H = values['excitation_complex']
    complex_sum = np.real((amplitude[:, None]*np.exp(1j*psi)).T@H)
    scalar_sum = np.array([[sum(amplitude[j]*(H[j, d].real*np.cos(psi[j, i])-H[j, d].imag*np.sin(psi[j, i]))
                               for j in range(384)) for d in range(6)] for i in range(3)])
    np.testing.assert_allclose(complex_sum, scalar_sum, rtol=2e-12, atol=1e-9)
    local_points = np.array([[0., 0., 0.], [2., -1., 0.], [-3., 2., 0.]])
    global_points = local_points@values['global_from_native_R'].T+values['global_from_native_t'].ravel()

    def eta(waves, points):
        w, direction, phase, amplitude = (np.asarray(waves[key]) for key in
            ('angular_frequencies', 'directions', 'phases', 'amplitudes'))
        vector = w[:, None]**2/gravity*np.column_stack((np.cos(direction), np.sin(direction)))
        return np.einsum('j,jpt->pt', amplitude, np.cos((vector@points[:, :2].T)[:, :, None]
            +phase[:, None, None]-w[:, None, None]*times))

    closure = eta(bank, local_points)-eta(global_bank['waves'], global_points)
    np.testing.assert_allclose(closure, 0, atol=2e-12)
    return dict(schema_exact=True, harmonic_signs=True,
        independent_FK_max_abs_error=float(np.max(abs(fk-values['fk_complex'][index]))),
        scalar_sum_max_abs_error=float(np.max(abs(complex_sum-scalar_sum))),
        frame_eta_max_abs_error=float(np.max(abs(closure))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--controls', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    controls = read(args.controls)
    if controls['status'] != 'ready_for_coordinator_execution_allocation':
        raise ValueError('Independent review and READY controls required')
    frozen = controls['input_sha256']
    for path, digest in frozen.items():
        assert sha(path) == digest, path
    started = time.perf_counter()
    args.output.mkdir(parents=True, exist_ok=False)
    snapshot = args.output/'source'
    snapshot.mkdir()
    for i, path in enumerate(dict.fromkeys([str(Path(__file__)), str(args.controls), *frozen])):
        shutil.copyfile(path, snapshot/f'{i:02d}-{Path(path).name}')
    report = dict(complete=False, expected_new_calls=450, completed_calls=0,
        command=[sys.executable, *sys.argv], source_sha256=sha(__file__), controls_sha256=sha(args.controls),
        input_sha256=frozen, cpu_threads=1, versions={name: version(name) for name in
            ('capytaine', 'numpy', 'scipy', 'threadpoolctl')}, scope='Diagnostic only; no physical qualification or canonical heading adoption')

    def save():
        (args.output/'manifest.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    save()
    try:
        import capytaine as cpt
        from capytaine.bem.airy_waves import froude_krylov_force
        from capytaine.green_functions.delhommeau import Delhommeau
        if version('capytaine') != '3.0.0':
            raise ValueError('Reviewed Capytaine3.0.0 required')
        original = read(controls['original_manifest'])
        package = distribution('capytaine')
        actual_dependency = {str(p): sha(package.locate_file(p)) for p in package.files
                             if str(p).endswith(('.py', '.pyd', '.so', '.dll'))}
        assert actual_dependency == original['dependency_sha256']['capytaine']
        mat = loadmat(controls['native_mat'])
        stored = loadmat(controls['radiation_mat'])
        case = read(controls['case'])
        assert np.array_equal(mat['M'], case['dry_mass_matrix_world'])
        with np.load(controls['native_npz']) as raw:
            base = raw['triangles'].copy()
        assert base.shape == (17, 3, 3)
        refined = triangles_at(base, 1)
        report['surface_invariants'] = surface_check(base, refined)
        rho, gravity = mat['rho'].item(), mat['g'].item()

        def body_for(triangles):
            vertices, ids = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
            mesh = cpt.Mesh(vertices=vertices, faces=ids.reshape(-1, 3), auto_clean=False, auto_check=False)
            return cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=mat['cg'].ravel()),
                center_of_mass=mat['cg'].ravel(), mass=case['body_mass_kg'], name=f'scout_{len(triangles)}')

        bodies = {17: body_for(base), 68: body_for(refined)}
        with threadpool_limits(1):
            solver = cpt.BEMSolver(green_function=Delhommeau(tabulation_cache_dir=None))

            def solve(problem):
                if time.perf_counter()-started > 300:
                    raise TimeoutError('Declared300s budget exceeded')
                result = solver.solve(problem, keep_details=False)
                report['completed_calls'] += 1
                if report['completed_calls'] % 6 == 0:
                    save()
                return result

            for panels, new in ((17, [16., 24., 32., 48., 64.]), (68, [1., 2., 12., 24., 48., np.inf])):
                matrices = {}
                if panels == 17:
                    for w in (1., 2., 12.):
                        index = np.flatnonzero(stored['omega_rad_s'].ravel() == w)
                        assert len(index) == 1
                        matrices[w] = (stored['added_mass'][index[0]].copy(), stored['radiation_damping'][index[0]].copy())
                    Ainf = mat['Ainf'].copy()
                for w in new:
                    A, B = np.empty((6, 6)), np.empty((6, 6))
                    for column, mode in enumerate(MODES):
                        result = solve(cpt.RadiationProblem(body=bodies[panels], omega=w,
                            radiating_dof=mode, rho=rho, g=gravity, water_depth=np.inf))
                        A[:, column] = [result.added_mass[row] for row in MODES]
                        if np.isfinite(w):
                            B[:, column] = [result.radiation_damping[row] for row in MODES]
                    if np.isfinite(w):
                        matrices[w] = (A, B)
                    else:
                        Ainf = A
                omega = np.array(sorted(matrices))
                A, B = np.stack([matrices[w][0] for w in omega]), np.stack([matrices[w][1] for w in omega])
                assert all(np.isfinite(v).all() for v in (A, B, Ainf))
                savemat(args.output/f'radiation-{panels}.mat', dict(omega_rad_s=omega[:, None],
                    added_mass=A, radiation_damping=B, Ainf=Ainf, dof_names=MODES,
                    M=mat['M'], K=mat['K'], cg=mat['cg'], rho=rho, g=gravity,
                    reused_frequency_mask=np.isin(omega, [1., 2., 12.]) if panels == 17 else np.zeros(len(omega), bool),
                    reused_Ainf=panels == 17, triangles=base if panels == 17 else refined))
                readback = loadmat(args.output/f'radiation-{panels}.mat')
                for key, value in dict(omega_rad_s=omega[:, None], added_mass=A,
                        radiation_damping=B, Ainf=Ainf, M=mat['M'], K=mat['K'], cg=mat['cg'],
                        triangles=base if panels == 17 else refined).items():
                    assert np.array_equal(readback[key], value), key
                print(f'Completed radiation scout{panels} panels', flush=True)
            assert report['completed_calls'] == 66
            heading = read(controls['heading_file'])
            assert heading['relative_heading_offset_deg'] == 90 and heading['observation_id'] == 'ndbc-46086-20250105T1140Z'
            bank = heading['boat_frame_waves']
            assert all(len(bank[key]) == 384 for key in ('amplitudes', 'angular_frequencies', 'directions', 'phases'))
            fk, diffraction = np.empty((384, 6), complex), np.empty((384, 6), complex)
            for i, (w, theta) in enumerate(zip(bank['angular_frequencies'], bank['directions'])):
                problem = cpt.DiffractionProblem(body=bodies[17], omega=w, wave_direction=theta,
                                                rho=rho, g=gravity, water_depth=np.inf)
                result = solve(problem)
                native_fk = froude_krylov_force(problem)
                fk[i], diffraction[i] = [native_fk[m] for m in MODES], [result.forces[m] for m in MODES]
            paired = args.output/'jan05-offset90-pairs.mat'
            savemat(paired, dict(component_index=np.arange(1, 385)[:, None],
                omega_rad_s=np.asarray(bank['angular_frequencies'])[:, None],
                heading_native_rad=np.asarray(bank['directions'])[:, None],
                amplitude_m=np.asarray(bank['amplitudes'])[:, None], phase_local_rad=np.asarray(bank['phases'])[:, None],
                fk_complex=fk, diffraction_complex=diffraction, excitation_complex=fk+diffraction,
                dof_names=MODES, cg=mat['cg'], rho=rho, g=gravity,
                yaw_global_rad=heading['mean_box_yaw_global_rad'],
                global_from_native_R=heading['global_from_boat_rotation'],
                global_from_native_t=heading['global_from_boat_translation_m'],
                globalbank_sha256=sha(controls['global_bank']), headingfile_sha256=sha(controls['heading_file']),
                convention='Native exp(-iomega*t); internal rejection hypothesis only, no canonical heading adoption'))
        report['pair_checks'] = check_pairs(paired, heading, read(controls['global_bank']), base, mat,
            dict(globalbank_sha256=sha(controls['global_bank']), headingfile_sha256=sha(controls['heading_file'])))
        assert report['completed_calls'] == 450
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
