"""Exactly3 prescribed native radiation witnesses, no motion/repair/band choice."""
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
from export_radiation_heading_scout import surface_check

MODES = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def work_fixture():
    # exp(-iwt): force=i*w*B*q for pure damping. This is an algebra fixture,
    # not another native solve or a synthetic body trajectory.
    omega, B = 2., 3.
    displacement = np.array([1., -.5])
    velocity = -1j*omega*displacement
    force = 1j*omega*B*displacement
    outward = -.5*np.real(np.vdot(velocity, force))
    np.testing.assert_allclose(outward, .5*omega**2*B*np.dot(displacement, displacement), rtol=0, atol=1e-14)
    A = np.array([[2., 3.], [1., 4.]])
    np.testing.assert_allclose(-.5*np.real(np.vdot(velocity, omega**2*A@displacement)), 0., rtol=0, atol=1e-14)
    assert outward > 0
    return dict(positive_damping_work_sign=True, real_motion_reactive_work_zero=True,
                convention='Ffluid=omega²Aq+iomegaBq; v=-iomegaq; Pout=-.5Re(v*F)')


def independent_fields(triangles, cg, coefficients):
    centers = triangles.mean(axis=1)
    area_normals = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])/2
    areas = np.linalg.norm(area_normals, axis=1)
    normals = area_normals/areas[:, None]
    displacement = coefficients[:3]+np.cross(coefficients[3:], centers-cg)
    return centers, areas, normals, displacement


def examine(solver, body, triangles, omega, coefficients, radiating_name, stored, out, on_solved):
    import capytaine as cpt
    from capytaine.post_pro.kochin import compute_kochin
    cg, rho, gravity = stored['cg'].ravel(), stored['rho'].item(), stored['g'].item()
    centers, areas, normals, displacement = independent_fields(triangles, cg, coefficients)
    np.testing.assert_allclose(body.mesh.faces_centers, centers, rtol=0, atol=2e-12)
    np.testing.assert_allclose(body.mesh.faces_normals, normals, rtol=0, atol=2e-12)
    np.testing.assert_allclose(body.mesh.faces_areas, areas, rtol=2e-13, atol=2e-12)
    problem = cpt.RadiationProblem(body=body, omega=omega, radiating_dof=radiating_name,
        rho=rho, g=gravity, free_surface=0., water_depth=np.inf, forward_speed=0.)
    normal_velocity = -1j*omega*np.sum(displacement*normals, axis=1)
    np.testing.assert_allclose(problem.boundary_condition, normal_velocity, rtol=2e-13, atol=2e-12)
    result = solver.solve(problem, method='indirect', keep_details=True)
    on_solved()
    pressure, potential, sources = result.pressure, result.potential, result.sources
    assert sources is not None and all(np.isfinite(v).all() for v in (pressure, potential, sources))
    np.testing.assert_allclose(pressure, 1j*omega*rho*potential, rtol=2e-13, atol=2e-12)
    # Independent physical force/moment integration, without integrate_pressure.
    face_force = -pressure[:, None]*areas[:, None]*normals
    physical_force = np.r_[face_force.sum(axis=0), np.cross(centers-cg, face_force).sum(axis=0)]
    native_force = np.array([result.forces[name] for name in MODES])
    np.testing.assert_allclose(physical_force, native_force, rtol=3e-12, atol=3e-9)
    generalized_force = np.dot(coefficients, physical_force)
    if radiating_name == 'WorkWitness':
        np.testing.assert_allclose(result.forces[radiating_name], generalized_force, rtol=3e-12, atol=3e-9)
    velocity = -1j*omega*coefficients
    force_work = float(-.5*np.real(np.vdot(velocity, physical_force)))
    panel_work = float(.5*np.real(np.sum(pressure*np.conj(normal_velocity)*areas)))
    native_work = float(.5*omega*np.imag(np.dot(coefficients, native_force)))
    np.testing.assert_allclose([force_work, native_work], panel_work, rtol=3e-12, atol=3e-9)
    details = dict(omega_rad_s=omega, panels=len(triangles), radiating_name=radiating_name,
        physical_coefficients_m_rad_per_dimensionless_amplitude=coefficients.tolist(),
        outward_panel_work_W_per_amplitude_squared=panel_work,
        outward_generalized_force_work_W_per_amplitude_squared=force_work,
        native_projected_work_W_per_amplitude_squared=native_work,
        raw_generalized_B_kg_m2_per_s=float(np.imag(generalized_force)/omega),
        physical_force_units=['N', 'N', 'N', 'Nm', 'Nm', 'Nm'],
        normalized_motion_amplitude=1., illustrative_small_amplitude=1e-3,
        outward_work_for_illustrative_small_amplitude_W=panel_work*1e-6,
        absolute_Kochin_power_claim=False)
    arrays = dict(omega_rad_s=omega, coefficients=coefficients, physical_triangles=triangles,
        centers=centers, normals=normals, areas=areas, displacement_field=displacement,
        normal_velocity=normal_velocity, boundary_condition=problem.boundary_condition,
        potential=potential, pressure=pressure, sources=sources, force_physical=native_force,
        force_independent=physical_force, M=stored['M'], K=stored['K'], cg=stored['cg'],
        cb=stored['cb'], volume=stored['volume'], rho=rho, g=gravity, dof_names=MODES)
    norms = []
    for count in (360, 720):
        theta = np.arange(count)*2*np.pi/count
        H = compute_kochin(result, theta, ref_point=tuple(cg[:2]))
        wavevector_projection = (centers[:, :2]-cg[:2])@np.vstack((np.cos(theta), np.sin(theta)))
        manual = (np.exp(result.wavenumber*centers[:, 2, None]-1j*result.wavenumber*wavevector_projection)
                  *areas[:, None]).T@sources/(4*np.pi)
        np.testing.assert_allclose(H, manual, rtol=3e-12, atol=3e-10)
        norm = float((2*np.pi/count)*np.sum(abs(H)**2))
        assert norm >= 0 and np.isfinite(norm)
        norms.append(norm)
        arrays[f'kochin_theta_{count}'], arrays[f'kochin_complex_{count}'] = theta, H
    details.update(Kochin_integral_abs_H_squared_360_720=norms,
        Kochin_angular_norm_difference=norms[1]-norms[0],
        Kochin_norm_units='m6/s2 per dimensionlessamplitude²; no power prefactor fitted or claimed')
    savemat(out/'details.mat', arrays)
    written = loadmat(out/'details.mat')
    for key in ('pressure', 'potential', 'sources', 'normal_velocity', 'force_physical', 'coefficients'):
        np.testing.assert_array_equal(written[key].ravel(), np.asarray(arrays[key]).ravel())
    (out/'work.json').write_text(json.dumps(details, indent=2, allow_nan=False)+'\n')
    return details, native_force


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
    report = dict(complete=False, expected_calls=3, completed_calls=0, cpu_threads=1,
        command=[sys.executable, *sys.argv], source_sha256=sha(__file__), controls_sha256=sha(args.controls),
        input_sha256=frozen, scope='Prescribed native harmonic energy witness only; no time trajectory/cutoff/repair/bug attribution', cases=[])

    def save():
        (args.output/'manifest.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    def on_solved():
        report['completed_calls'] += 1
        save()

    save()
    try:
        import capytaine as cpt
        from capytaine.green_functions.delhommeau import Delhommeau
        assert version('capytaine') == '3.0.0'
        report['version'] = version('capytaine')
        package = distribution('capytaine')
        hashes = {str(p): sha(package.locate_file(p)) for p in package.files
                  if str(p).endswith(('.py', '.pyd', '.so', '.dll'))}
        assert hashes == read(controls['original_manifest'])['dependency_sha256']['capytaine']
        report['algebra_fixture'] = work_fixture()
        stored = loadmat(controls['native_mat'])
        raw = loadmat(controls['radiation_mat'])
        with np.load(controls['native_npz']) as archive:
            base = archive['triangles'].copy()
        assert base.shape == (17, 3, 3)
        refined = triangles_at(base, 1)
        np.testing.assert_array_equal(refined, loadmat(controls['scout68_mat'])['triangles'])
        report['surface_invariants'] = surface_check(base, refined)
        case = read(controls['case'])
        assert np.array_equal(stored['M'], case['dry_mass_matrix_world'])
        cg = stored['cg'].ravel()
        coefficients = np.array(controls['fixed_physical_coefficients_m_rad'])
        scaled = np.array(read(controls['stored_witness'])['incident_band_native_B_symmetric_energy']['eigenvector_scaled_coordinates'])
        np.testing.assert_array_equal(coefficients, np.diag([2.7432]*3+[1.]*3)@scaled)
        with threadpool_limits(1):
            solver = cpt.BEMSolver(green_function=Delhommeau(tabulation_cache_dir=None), method='indirect')
            definitions = [('fixed17_omega2_2', base, 2.2, coefficients, 'WorkWitness'),
                           ('heave17_omega1_positive_control', base, 1., np.array([0., 0., 1., 0., 0., 0.]), 'Heave'),
                           ('same_fixed68_omega2_2', refined, 2.2, coefficients, 'WorkWitness')]
            for label, triangles, omega, prescribed, mode in definitions:
                if time.perf_counter()-started > 180:
                    raise TimeoutError('Declared180s elapsed budget exceeded')
                vertices, ids = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
                mesh = cpt.Mesh(vertices=vertices, faces=ids.reshape(-1, 3), auto_clean=False, auto_check=False)
                body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=cg),
                    center_of_mass=cg, mass=case['body_mass_kg'], name=label)
                if mode == 'WorkWitness':
                    body.dofs['WorkWitness'] = coefficients[:3]+np.cross(coefficients[3:], mesh.faces_centers-cg)
                assert body.lid_mesh is None and body.mesh.nb_faces == len(triangles)
                out = args.output/label
                out.mkdir()
                details, force = examine(solver, body, triangles, omega, prescribed, mode, stored, out, on_solved)
                if len(triangles) == 17:
                    index = np.flatnonzero(raw['omega_rad_s'].ravel() == omega).item()
                    A, B = raw['added_mass'][index], raw['radiation_damping'][index]
                    expected_force = (omega**2*A+1j*omega*B)@prescribed
                    np.testing.assert_allclose(force, expected_force, rtol=2e-9, atol=2e-8)
                    expected_work = float(.5*omega**2*(prescribed@B@prescribed))
                    np.testing.assert_allclose(details['outward_panel_work_W_per_amplitude_squared'], expected_work,
                                               rtol=2e-9, atol=2e-8)
                    details['exact_stored_omega_used'] = float(raw['omega_rad_s'][index, 0])
                    details['stored_raw_B_work_W_per_amplitude_squared'] = expected_work
                (out/'work.json').write_text(json.dumps(details, indent=2, allow_nan=False)+'\n')
                report['cases'].append(details)
                save()
                print(f'Completed prescribed native witness {label}', flush=True)
        assert report['completed_calls'] == 3 and all(sha(path) == digest for path, digest in frozen.items())
        report.update(complete=True, elapsed_seconds=time.perf_counter()-started,
            output_sha256={p.relative_to(args.output).as_posix(): sha(p) for p in args.output.rglob('*')
                           if p.is_file() and 'source' not in p.relative_to(args.output).parts and p.name != 'manifest.json'})
    except BaseException as error:
        report.update(error=repr(error), elapsed_seconds=time.perf_counter()-started)
    save()
    print(json.dumps({k: report.get(k) for k in ('complete', 'completed_calls', 'elapsed_seconds', 'error')}))
    return 0 if report['complete'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
