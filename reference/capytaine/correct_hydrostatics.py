"""Recompute native hydrostatics on unchanged triangles; reuse every BEM solve.

Capytaine 3 native Gauss-Legendre 2 integrates these polynomial hydrostatic
moments exactly on the padded triangular faces. This does not refine the mesh
or replace native hydrodynamics. Run only under an allocated numerical lease.
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

for thread_variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[thread_variable] = '1'

import numpy as np
import xarray as xr
from scipy.io import loadmat, savemat
from threadpoolctl import threadpool_limits

MODES = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']
HYDRODYNAMIC_FIELDS = ('added_mass', 'radiation_damping', 'diffraction_force',
                      'Froude_Krylov_force', 'excitation_force')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_sha(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, default=Path('reference/capytaine/canonical-boat.json'))
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source_manifest = json.loads((args.source/'manifest.json').read_text())
    if not source_manifest['complete'] or source_manifest['completed_calls'] != 1086:
        raise ValueError('Source is not the completed minimum producer export')
    if sha(args.case) != source_manifest['case_sha256']:
        raise ValueError('Physical input differs from original hydrodynamic solve')
    for name, digest in source_manifest['output_sha256'].items():
        if sha(args.source/name) != digest:
            raise ValueError('Source output identity changed: '+name)
    args.output.mkdir(parents=True, exist_ok=False)
    snapshot = args.output/'source'
    snapshot.mkdir()
    for path in (Path(__file__), args.case):
        shutil.copyfile(path, snapshot/path.name)
    shutil.copyfile(args.source/'manifest.json', snapshot/'original-native-manifest.json')
    started = time.perf_counter()
    report = dict(complete=False, command=[sys.executable, *sys.argv],
                  scope='Native hydrostatics-only correction and version-interface alias; no new BEM solve',
                  physical_input_sha256=sha(args.case), correction_sha256=sha(Path(__file__)),
                  source_directory=str(args.source), source_manifest_sha256=sha(args.source/'manifest.json'),
                  reused_output_sha256=source_manifest['output_sha256'],
                  versions={name: version(name) for name in
                            ('capytaine', 'numpy', 'scipy', 'xarray', 'threadpoolctl', 'pandas')},
                  hydrostatic_quadrature='Capytaine native Gauss-Legendre 2 on unchanged padded triangular faces',
                  hydrodynamic_quadrature='Original centroid/collocation results, unchanged',
                  native_thread_count=1, new_BEM_solves=0)

    def save():
        (args.output/'manifest.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    save()
    try:
        import capytaine as cpt
        from capytaine.io.legacy import export_hydrostatics
        if version('capytaine') != source_manifest['versions']['capytaine']:
            raise ValueError('Capytaine version differs from original solve')
        package = distribution('capytaine')
        # Match every original Capytaine implementation/binary, not only version.
        hashes = {str(path): sha(package.locate_file(path)) for path in package.files
                  if str(path).endswith(('.py', '.pyd', '.so', '.dll'))}
        if hashes != source_manifest['dependency_sha256']['capytaine']:
            raise ValueError('Capytaine implementation differs from original solve')
        case = json.loads(args.case.read_text())
        with np.load(args.source/'native.npz') as archive:
            raw = {key: archive[key] for key in archive.files}
        mat = {key: value for key, value in loadmat(args.source/'native.mat').items()
               if not key.startswith('__')}
        if not np.array_equal(mat['M'], case['dry_mass_matrix_world']):
            raise ValueError('Source full dry inertia differs from canonical physical input')
        triangles = raw['triangles']
        if triangles.shape != (17, 3, 3):
            raise ValueError('Hydrostatic correction requires the unchanged minimum17 triangles')
        vertices, indices = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
        mesh = cpt.Mesh(vertices=vertices, faces=indices.reshape(-1, 3),
                        quadrature_method='Gauss-Legendre 2', auto_clean=False, auto_check=False)
        body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=case['world_cg_m']),
                               center_of_mass=case['world_cg_m'], mass=case['body_mass_kg'], name='canonical_boat')
        rho = case['water']['density_kg_m3']
        gravity = case['water']['gravity_m_s2']
        body.inertia_matrix = xr.DataArray(case['dry_mass_matrix_world'],
            dims=['influenced_dof', 'radiating_dof'],
            coords=dict(influenced_dof=MODES, radiating_dof=MODES))
        with threadpool_limits(1):
            body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness(rho=rho, g=gravity)
            stiffness = body.hydrostatic_stiffness.sel(influenced_dof=MODES, radiating_dof=MODES).values
            center = body.center_of_buoyancy
            volume = body.disp_volume
            export_hydrostatics(str(args.output), body)
        with xr.open_dataset(args.source/'canonical.nc', engine='scipy') as original:
            data = original.load().copy(deep=True)
        # Change only physical hydrostatic outputs; retain every frequency-domain array.
        data['hydrostatic_stiffness'].data = stiffness
        data['center_of_buoyancy'].data = center
        data['disp_mass'].data = np.asarray(rho*volume)
        # Official WEC-Sim v7.1 readCAPYTAINE requires this old coordinate name.
        # Capytaine3 exports `body`; preserve it and add an explicit compatibility alias.
        if 'body_name' not in data:
            data['body_name'] = data['body'].copy()
        data.attrs['hydrostatic_quadrature'] = report['hydrostatic_quadrature']
        data.attrs['hydrodynamic_source_manifest_sha256'] = report['source_manifest_sha256']
        with xr.set_options(netcdf_engine_order=['scipy']):
            cpt.export_dataset(args.output/'canonical.nc', data)
        raw.update(hydrostatic_stiffness=stiffness, center_of_buoyancy=center,
                   displaced_volume=np.array(volume))
        mat.update(K=stiffness, cb=center, volume=volume)
        np.savez(args.output/'native.npz', **raw)
        savemat(args.output/'native.mat', mat)
        shutil.copyfile(args.source/'boat.stl', args.output/'boat.stl')
        with xr.open_dataset(args.source/'canonical.nc', engine='scipy') as original, \
             xr.open_dataset(args.output/'canonical.nc', engine='scipy') as corrected:
            reused = {key: bool(np.array_equal(original[key], corrected[key])) for key in original.variables
                      if key not in ('hydrostatic_stiffness', 'center_of_buoyancy', 'disp_mass')}
            if not all(reused.values()):
                raise ValueError('Correction changed a reused native field')
            coefficient_hashes = {key: array_sha(corrected[key].values) for key in HYDRODYNAMIC_FIELDS}
        if not all(np.isfinite(value).all() for value in (stiffness, center, volume)):
            raise ValueError('Nonfinite corrected hydrostatics')
        report.update(complete=True, elapsed_seconds=time.perf_counter()-started,
                      reused_variables_exact=reused, hydrodynamic_array_sha256=coefficient_hashes,
                      original_center_of_buoyancy_m=loadmat(args.source/'native.mat')['cb'].ravel().tolist(),
                      corrected_center_of_buoyancy_m=np.asarray(center).tolist(),
                      corrected_volume_m3=float(volume), corrected_stiffness=stiffness.tolist(),
                      output_sha256={path.name: sha(path) for path in args.output.iterdir()
                                     if path.is_file() and path.name != 'manifest.json'})
    except BaseException as error:
        report.update(error=dict(type=type(error).__name__, message=str(error)),
                      elapsed_seconds=time.perf_counter()-started)
    save()
    print(json.dumps({key: report.get(key) for key in
                     ('complete', 'new_BEM_solves', 'elapsed_seconds', 'error')}))
    return 0 if report['complete'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
