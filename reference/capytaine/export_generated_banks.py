"""Native radiation and exact paired excitation; requires a separate execution grant.

No Cartesian assembly of the component pairs, wave modification, response solver,
or coefficient symmetry repair. Existing baseline exports remain immutable.
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
import xarray as xr
from scipy.io import loadmat, savemat
from threadpoolctl import threadpool_limits

MODES = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    contract = read(args.contract)
    if contract['review_status'] != 'ready_for_coordinator_execution_allocation':
        raise ValueError('Joint source/interface review not ready')
    frozen = contract['input_sha256']
    for path, digest in frozen.items():
        if sha(path) != digest:
            raise ValueError('Frozen input differs: '+path)
    case_path = Path(contract['case'])
    case = read(case_path)
    native = Path(contract['accepted_native'])
    accepted = read(native/'manifest.json')
    if not accepted['complete']:
        raise ValueError('Accepted native source incomplete')
    for path, digest in accepted['output_sha256'].items():
        if sha(native/path) != digest:
            raise ValueError('Accepted native source changed: '+path)
    mat = {key: value for key, value in loadmat(native/'native.mat').items()
           if not key.startswith('__')}
    if not np.array_equal(mat['M'], case['dry_mass_matrix_world']):
        raise ValueError('Full dry matrix identity differs')
    started = time.perf_counter()
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    source = out/'source'
    source.mkdir()
    paths = [Path(__file__), args.contract, *map(Path, frozen)]
    for index, path in enumerate(dict.fromkeys(paths)):
        shutil.copyfile(path, source/f'{index:02d}-{path.name}')
    report = dict(complete=False, command=[sys.executable, *sys.argv],
                  input_sha256=frozen, source_sha256=sha(__file__),
                  contract_sha256=sha(args.contract), completed_calls=0,
                  expected_calls=2838, panels=17, cpu_threads=1,
                  scope='Native coefficients only; no motion or mesh-accuracy claim',
                  versions={name: version(name) for name in
                            ('capytaine', 'numpy', 'scipy', 'xarray', 'threadpoolctl')})

    def save():
        (out/'manifest.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    save()
    try:
        import capytaine as cpt
        from capytaine.bem.airy_waves import froude_krylov_force
        from capytaine.green_functions.delhommeau import Delhommeau
        if version('capytaine') != '3.0.0':
            raise ValueError('Reviewed API is Capytaine3.0.0')
        # Exact baseline hydrodynamic geometry and centroid quadrature: GL2 is
        # used only for previously accepted hydrostatic metadata, copied below.
        original = read(native/'source'/'original-native-manifest.json')
        report['dependency_sha256'] = {}
        for name in report['versions']:
            package = distribution(name)
            report['dependency_sha256'][name] = {
                str(path): sha(package.locate_file(path)) for path in package.files
                if str(path).endswith(('.py', '.pyd', '.so', '.dll'))}
        if report['dependency_sha256']['capytaine'] != original['dependency_sha256']['capytaine']:
            raise ValueError('Capytaine source/binary differs from baseline')
        with np.load(native/'native.npz') as archive:
            triangles = archive['triangles'].copy()
        if triangles.shape != (17, 3, 3):
            raise ValueError('Expected unchanged17 triangles')
        vertices, indices = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
        mesh = cpt.Mesh(vertices=vertices, faces=indices.reshape(-1, 3),
                        auto_clean=False, auto_check=False)
        body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=case['world_cg_m']),
                               center_of_mass=case['world_cg_m'], mass=case['body_mass_kg'], name='canonical_boat')
        coordinates = dict(influenced_dof=MODES, radiating_dof=MODES)
        body.inertia_matrix = xr.DataArray(mat['M'], dims=list(coordinates), coords=coordinates)
        body.hydrostatic_stiffness = xr.DataArray(mat['K'], dims=list(coordinates), coords=coordinates)
        rho, gravity = float(mat['rho'].item()), float(mat['g'].item())
        common = dict(body=body, rho=rho, g=gravity, water_depth=np.inf)
        frequencies = np.arange(1, 241, dtype=float)*0.05
        results = []
        with threadpool_limits(1):
            solver = cpt.BEMSolver(green_function=Delhommeau(tabulation_cache_dir=None))

            def solve(problem):
                if time.perf_counter()-started > contract['maximum_elapsed_seconds']:
                    raise TimeoutError('Bounded producer elapsed budget exceeded')
                result = solver.solve(problem, keep_details=False)
                report['completed_calls'] += 1
                if report['completed_calls'] % 60 == 0:
                    save()
                    print(f"Native calls {report['completed_calls']}/2838", flush=True)
                return result

            Ainf = np.empty((6, 6))
            for column, mode in enumerate(MODES):
                result = solve(cpt.RadiationProblem(omega=np.inf, radiating_dof=mode, **common))
                Ainf[:, column] = [result.added_mass[row] for row in MODES]
            for omega in frequencies:
                for mode in MODES:
                    results.append(solve(cpt.RadiationProblem(omega=omega, radiating_dof=mode, **common)))
                # Real native heading0 data solely for official BEMIO schema.
                results.append(solve(cpt.DiffractionProblem(omega=omega, wave_direction=0., **common)))
            data = cpt.assemble_dataset(results, hydrostatics=True).sel(**coordinates)
            data['center_of_buoyancy'].data = mat['cb'].ravel()
            data['disp_mass'].data = np.asarray(rho*mat['volume'].item())
            data['body_name'] = data['body'].copy()
            data.attrs['excitation_role'] = 'Actual heading0 only; generated forcing uses separate paired tables'
            for name in ('added_mass', 'radiation_damping', 'diffraction_force',
                         'Froude_Krylov_force', 'excitation_force', 'inertia_matrix', 'hydrostatic_stiffness'):
                if not np.isfinite(data[name]).all():
                    raise ValueError('Nonfinite native field '+name)
            with xr.set_options(netcdf_engine_order=['scipy']):
                cpt.export_dataset(out/'canonical.nc', data)
            radiation = dict(omega_rad_s=frequencies[:, None],
                             added_mass=data.added_mass.transpose('omega', 'influenced_dof', 'radiating_dof').values,
                             radiation_damping=data.radiation_damping.transpose('omega', 'influenced_dof', 'radiating_dof').values)
            savemat(out/'radiation.mat', dict(**radiation, Ainf=Ainf, M=mat['M'], K=mat['K'], dof_names=MODES))
            for record in contract['records']:
                heading_path = Path(record['heading_file'])
                heading = read(heading_path)
                bank = heading['boat_frame_waves']
                if heading['relative_heading_offset_deg'] != 0 or not heading['default_heading_case']:
                    raise ValueError('Only unchanged default offset0 bank is authorized')
                omega = np.asarray(bank['angular_frequencies'])
                theta = np.asarray(bank['directions'])
                amplitude = np.asarray(bank['amplitudes'])
                phase = np.asarray(bank['phases'])
                if not all(array.shape == (384,) for array in (omega, theta, amplitude, phase)):
                    raise ValueError('Expected every one of384 bank components')
                fk, diffraction = np.empty((384, 6), complex), np.empty((384, 6), complex)
                for index, (w, direction) in enumerate(zip(omega, theta)):
                    problem = cpt.DiffractionProblem(omega=w, wave_direction=direction, **common)
                    result = solve(problem)
                    native_fk = froude_krylov_force(problem)
                    fk[index] = [native_fk[mode] for mode in MODES]
                    diffraction[index] = [result.forces[mode] for mode in MODES]
                excitation = fk+diffraction
                if not all(np.isfinite(array).all() for array in (fk, diffraction, excitation)):
                    raise ValueError('Nonfinite paired coefficients')
                savemat(out/(record['observation_id']+'-pairs.mat'), dict(
                    component_index=np.arange(1, 385)[:, None], omega_rad_s=omega[:, None],
                    heading_native_rad=theta[:, None], amplitude_m=amplitude[:, None],
                    phase_local_rad=phase[:, None], fk_complex=fk, diffraction_complex=diffraction,
                    excitation_complex=excitation, dof_names=MODES, cg=mat['cg'], rho=rho, g=gravity,
                    globalbank_sha256=sha(heading_path.parent/'realization-global.json'),
                    headingfile_sha256=sha(heading_path), yaw_global_rad=heading['mean_box_yaw_global_rad'],
                    global_from_native_R=heading['global_from_boat_rotation'],
                    global_from_native_t=heading['global_from_boat_translation_m'],
                    convention='Dimensional perunit elevation; exp(-i omega t); no conjugation'))
                # Serialized table identity is part of the producer interface.
                written = loadmat(out/(record['observation_id']+'-pairs.mat'))
                for key, expected in dict(omega_rad_s=omega[:, None],
                                         component_index=np.arange(1, 385)[:, None],
                                         heading_native_rad=theta[:, None], amplitude_m=amplitude[:, None],
                                         phase_local_rad=phase[:, None], fk_complex=fk,
                                         diffraction_complex=diffraction, excitation_complex=excitation,
                                         cg=mat['cg'], rho=np.array([[rho]]), g=np.array([[gravity]]),
                                         yaw_global_rad=np.array([[heading['mean_box_yaw_global_rad']]]),
                                         global_from_native_R=np.asarray(heading['global_from_boat_rotation']),
                                         global_from_native_t=np.atleast_2d(heading['global_from_boat_translation_m'])).items():
                    if not np.array_equal(written[key], expected):
                        raise ValueError('MAT paired readback differs: '+key)
                if np.char.strip(written['dof_names'].astype(str)).ravel().tolist() != MODES:
                    raise ValueError('MAT paired DOF order differs')
                for key, expected in dict(globalbank_sha256=sha(heading_path.parent/'realization-global.json'),
                                          headingfile_sha256=sha(heading_path)).items():
                    if written[key].item().strip() != expected:
                        raise ValueError('MAT paired source identity differs: '+key)
        if report['completed_calls'] != report['expected_calls'] or not np.isfinite(Ainf).all():
            raise ValueError('Incomplete call count or Ainf')
        mat['Ainf'] = Ainf
        savemat(out/'native.mat', mat)
        for name in ('boat.stl', 'Hydrostatics.dat', 'KH.dat'):
            shutil.copyfile(native/name, out/name)
        for path, digest in frozen.items():
            if sha(path) != digest:
                raise ValueError('Frozen input changed during solve: '+path)
        report.update(complete=True, elapsed_seconds=time.perf_counter()-started,
                      output_sha256={path.name: sha(path) for path in out.iterdir()
                                     if path.is_file() and path.name != 'manifest.json'})
    except BaseException as error:
        report.update(error=repr(error), elapsed_seconds=time.perf_counter()-started)
    save()
    print(json.dumps({key: report.get(key) for key in ('complete', 'completed_calls', 'elapsed_seconds', 'error')}))
    return 0 if report['complete'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
