"""Post-export identity, phase and independent incident-pressure checks; no BEM."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
import xarray as xr


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--proposal', type=Path, required=True)
    args = parser.parse_args()
    out = args.output
    manifest = read(out/'manifest.json')
    proposal = read(args.proposal)
    contract = read('reference/capytaine/generated-bank-contract.json')
    checks = {}
    assert manifest['complete'] and manifest['completed_calls'] == 2838
    checks['terminal_manifest_call_count'] = True
    assert all(sha(path) == digest for path, digest in proposal['input_sha256'].items())
    assert all(sha(out/name) == digest for name, digest in manifest['output_sha256'].items())
    checks['frozen_input_and_output_hashes'] = True
    native = Path(contract['accepted_native'])
    original, actual = loadmat(native/'native.mat'), loadmat(out/'native.mat')
    for key in ('M', 'K', 'cg', 'cb', 'volume', 'rho', 'g'):
        assert np.array_equal(original[key], actual[key]), key
    assert actual['Ainf'].shape == (6, 6) and np.isfinite(actual['Ainf']).all()
    checks['exact_full_mean_metadata_and_finite_Ainf'] = True
    with xr.open_dataset(out/'canonical.nc', engine='scipy') as dataset:
        assert dataset.sizes['omega'] == 240 and dataset.sizes['wave_direction'] == 1
        for key in ('added_mass', 'radiation_damping', 'diffraction_force',
                    'Froude_Krylov_force', 'excitation_force', 'inertia_matrix', 'hydrostatic_stiffness'):
            assert np.isfinite(dataset[key]).all(), key
        assert 'body_name' in dataset
        assert list(dataset.influenced_dof.values) == ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']
        assert list(dataset.radiating_dof.values) == list(dataset.influenced_dof.values)
    checks['official_schema_full_finite_arrays_and_dof_order'] = True
    # Pure analytic signs for exp(-iwt); H=i has +1 response at t=T/4.
    assert abs(np.real(1*np.exp(0j))-1) < 1e-15
    assert abs(np.real(1j*np.exp(-1j*np.pi/2))-1) < 1e-15
    checks['analytic_quarter_cycle_signs'] = True
    with np.load(native/'native.npz') as archive:
        triangles = archive['triangles']
    assert triangles.shape == (17, 3, 3)
    area_normals = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])/2
    centers = triangles.mean(axis=1)
    cg = actual['cg'].ravel()
    rho, gravity = actual['rho'].item(), actual['g'].item()
    records = []
    for record in contract['records']:
        heading = read(record['heading_file'])
        bank = heading['boat_frame_waves']
        values = loadmat(out/(record['observation_id']+'-pairs.mat'))
        for field, source in (('omega_rad_s', 'angular_frequencies'), ('heading_native_rad', 'directions'),
                              ('amplitude_m', 'amplitudes'), ('phase_local_rad', 'phases')):
            assert np.array_equal(values[field].ravel(), bank[source]), field
        assert np.array_equal(values['component_index'].ravel(), np.arange(1, 385))
        assert np.array_equal(values['excitation_complex'], values['fk_complex']+values['diffraction_complex'])
        assert all(values[key].shape == (384, 6) and np.isfinite(values[key]).all()
                   for key in ('excitation_complex', 'fk_complex', 'diffraction_complex'))
        index = int(np.flatnonzero(values['amplitude_m'].ravel() > 0)[0])
        omega, theta = values['omega_rad_s'][index, 0], values['heading_native_rad'][index, 0]
        k = omega**2/gravity
        pressure = rho*gravity*np.exp(k*centers[:, 2]+1j*k*(centers[:, 0]*np.cos(theta)+centers[:, 1]*np.sin(theta)))
        forces = -pressure[:, None]*area_normals
        expected_fk = np.r_[forces.sum(axis=0), np.cross(centers-cg, forces).sum(axis=0)]
        error = float(np.max(abs(expected_fk-values['fk_complex'][index])))
        np.testing.assert_allclose(expected_fk, values['fk_complex'][index], atol=1e-9, rtol=2e-12)
        times = np.array([0., .37, 3.1])
        amplitude = values['amplitude_m'].ravel()
        phase = values['phase_local_rad'].ravel()
        frequencies = values['omega_rad_s'].ravel()
        H = values['excitation_complex']
        matrix = np.real((amplitude[:, None]*np.exp(1j*(phase[:, None]-frequencies[:, None]*times))) .T@H)
        scalar = np.array([[sum(amplitude[j]*(H[j, dof].real*np.cos(phase[j]-frequencies[j]*t)
                           -H[j, dof].imag*np.sin(phase[j]-frequencies[j]*t)) for j in range(384))
                           for dof in range(6)] for t in times])
        np.testing.assert_allclose(matrix, scalar, rtol=2e-12, atol=1e-9)
        global_waves = read(Path(record['heading_file']).parent/'realization-global.json')['waves']
        points = np.array([[0., 0., 0.], [2., -1., 0.], [-3., 2., 0.]])
        R, translation = values['global_from_native_R'], values['global_from_native_t'].ravel()
        global_points = points@R.T+translation

        def eta(waves, positions):
            w, direction, p, a = (np.asarray(waves[key]) for key in
                                  ('angular_frequencies', 'directions', 'phases', 'amplitudes'))
            wavevector = w[:, None]**2/gravity*np.column_stack((np.cos(direction), np.sin(direction)))
            return np.einsum('j,jpt->pt', a, np.cos((wavevector@positions[:, :2].T)[:, :, None]
                +p[:, None, None]-w[:, None, None]*times))

        # Explicit broadcasting avoids hiding a spatial-origin correction.
        local_eta, global_eta = eta(bank, points), eta(global_waves, global_points)
        np.testing.assert_allclose(local_eta, global_eta, rtol=0, atol=2e-12)
        records.append(dict(id=record['observation_id'], independent_fk_max_abs_error=error,
                            component=index+1, frame_eta_max_abs_error=float(np.max(abs(local_eta-global_eta))),
                            scalar_matrix_max_abs_error=float(np.max(abs(matrix-scalar)))))
    checks['all_exact_pair_arrays_indices_finite_and_sum'] = True
    checks['independent_actual_incident_pressure_FK'] = True
    checks['multidirectional_scalar_matrix_phase_reconstruction'] = True
    checks['global_native_spatial_eta_frame_closure'] = True
    report = dict(complete=True, checks=checks, records=records, check_source_sha256=sha(__file__),
                  manifest_sha256=sha(out/'manifest.json'), proposal_sha256=sha(args.proposal),
                  limitation='Coefficient/interface checks only; no mesh/kernel/motion/relative-geometry accuracy claim')
    (out/'postcheck.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
