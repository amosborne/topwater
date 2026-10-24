"""Analytic incident-field diagnostics only; no generation or boat response."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIELDS = ('eta_m', 'eta_x', 'eta_y', 'eta_xx_per_m', 'eta_xy_per_m',
          'eta_yy_per_m', 'eta_t_m_per_s')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def footprint(physical):
    hull = physical['minimum_hull']
    w, x, y, z = np.asarray(hull['equilibrium_quaternion_wxyz'], dtype=float)
    np.testing.assert_allclose(w*w+x*x+y*y+z*z, 1., atol=2e-15, rtol=0)
    rotation = np.array([[1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
                         [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
                         [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)]])
    cg = np.asarray(physical['world_cg_m'])
    vertices = (np.asarray(hull['vertices_material_m'])-hull['cg_material_m'])@rotation.T+cg
    triangles = vertices[np.asarray(hull['triangular_faces'])]
    samples = [vertices[:, :2], triangles.mean(axis=1)[:, :2]]
    samples += [((triangles[:, i]+triangles[:, (i+1) % 3])/2)[:, :2] for i in range(3)]
    return cg[:2], np.unique(np.concatenate(samples), axis=0), vertices[:, :2]


def fields(amplitude, omega, direction, phase, points, times, gravity):
    """Re sum a exp(i(k.x-wt+phase)); exact derivative multipliers."""
    k = omega**2/gravity
    kx, ky = k*np.cos(direction), k*np.sin(direction)
    spatial = amplitude[None, :]*np.exp(1j*(points[:, 0, None]*kx+points[:, 1, None]*ky+phase))
    multipliers = np.array([np.ones_like(k), 1j*kx, 1j*ky,
                            -kx*kx, -kx*ky, -ky*ky, -1j*omega])
    # Equal frequencies are grouped algebraically; every original component remains.
    frequencies, inverse = np.unique(omega, return_inverse=True)
    coefficients = np.zeros((len(FIELDS), len(points), len(frequencies)), dtype=complex)
    for band in range(len(frequencies)):
        selected = inverse == band
        coefficients[:, :, band] = multipliers[:, selected]@spatial[:, selected].T
    temporal = np.exp(-1j*frequencies[:, None]*times)
    return (coefficients.reshape(-1, len(frequencies))@temporal).real.reshape(len(FIELDS), len(points), -1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    config = read(args.config)
    bank = ROOT/config['bank_directory']
    bank_manifest = read(bank/'manifest.json')
    assert bank_manifest['complete'] is True
    physical_path = ROOT/config['physical_input']
    physical = read(physical_path)
    assert sha(physical_path) == bank_manifest['input_sha256'][config['physical_input']]
    cg, sample_xy, vertex_xy = footprint(physical)
    points = np.vstack((cg, sample_xy))
    offsets = points-cg
    gravity = physical['water']['gravity_m_s2']
    times = np.arange(config['samples'], dtype=float)*config['dt_s']
    assert times[0] == 0 and times[-1] == config['duration_s']
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError('Use a fresh output directory; preserve prior evidence')
    inputs = [Path(__file__).resolve(), args.config.resolve(), physical_path, bank/'manifest.json']
    for record in config['observation_ids']:
        inputs += [bank/record/'summary.json', bank/record/'realization-global.json']
        inputs += [bank/record/f'heading-offset-{offset:g}deg.json'
                   for offset in config['relative_heading_offsets_deg']]
    binaries = list((Path(np.__file__).parent/'_core').glob('*.pyd'))
    frozen = {str(path.resolve()): sha(path) for path in inputs+binaries}
    # Verify accepted immutable exported inputs before evaluating any field.
    for path in inputs:
        if path.is_relative_to(bank) and path.name != 'manifest.json':
            assert sha(path) == bank_manifest['output_sha256'][path.relative_to(bank).as_posix()]
    output.mkdir(parents=True)
    for path in (Path(__file__).resolve(), args.config.resolve(), physical_path):
        archived = output/'source'/path.relative_to(ROOT)
        archived.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, archived)
    manifest = dict(complete=False, scope='Incident field ONLY: slope angles are input proxies, never boat response',
        started_at=datetime.now(timezone.utc).isoformat(), command=sys.argv, python=sys.version,
        python_executable_sha256=sha(sys.executable), numpy_version=np.__version__, frozen_sha256=frozen,
        fields=list(FIELDS), field_convention='eta=sum a cos(kx*x+ky*y-omega*t+phase); Hessian is graph second derivative, not principal surface curvature',
        time_samples=len(times), dt_s=config['dt_s'], duration_s=config['duration_s'],
        footprint_sampling='CG then unique projected full equilibrium-hull vertices, triangle centroids and edge midpoints; extrema are sampled, not continuous extrema',
        initialization='Original stationary banks at t=0; no rephase, scale, force ramp, pose or radiation-state prescription',
        selection='Internal candidates only; no canonical choice or motion qualification from slope proxies', cases=[])
    write(output/'manifest.json', manifest)
    try:
        for record in config['observation_ids']:
            source_summary = read(bank/record/'summary.json')
            for offset in config['relative_heading_offsets_deg']:
                case = read(bank/record/f'heading-offset-{offset:g}deg.json')
                waves = case['boat_frame_waves']
                a, omega, direction, phase = [np.asarray(waves[name]) for name in
                    ('amplitudes', 'angular_frequencies', 'directions', 'phases')]
                values = fields(a, omega, direction, phase, points, times, gravity)
                assert np.isfinite(values).all()
                # Independent direct ungrouped t=0 identities for every field/point.
                k = omega**2/gravity
                kx, ky = k*np.cos(direction), k*np.sin(direction)
                arguments = points@np.vstack((kx, ky))+phase
                cosine, sine = np.cos(arguments), np.sin(arguments)
                direct = np.array([cosine@a, -sine@(a*kx), -sine@(a*ky),
                    -cosine@(a*kx*kx), -cosine@(a*kx*ky), -cosine@(a*ky*ky), sine@(a*omega)])
                tolerance = 4096*np.finfo(float).eps*max(1., float(np.sum(a*(1+k*k+omega))))
                np.testing.assert_allclose(values[:, :, 0], direct, rtol=0, atol=tolerance)
                plane = values[0, :1]+offsets[:, 0, None]*values[1, :1]+offsets[:, 1, None]*values[2, :1]
                departure = values[0]-plane
                proxy = np.rad2deg(np.arctan(values[[2, 1], 0]))
                proxy[1] *= -1
                # Global-in-time continuous Taylor bound applies over the entire hull footprint.
                max_radius_squared = float(np.max(np.sum((vertex_xy-cg)**2, axis=1)))
                hessian_bound = float(np.sum(a*k*k))
                bound = .5*max_radius_squared*hessian_bound
                assert float(np.max(abs(departure))) <= bound+tolerance
                weights = a*a/2
                vector = np.sum(weights*np.exp(1j*direction))
                prefix = f'{record}-offset-{offset:g}deg'
                np.savez(output/f'{prefix}.npz', time_s=times, point_xy_body_m=points,
                    fields=values, cg_tangent_plane_eta_m=plane, departure_from_cg_plane_m=departure,
                    cg_input_slope_proxy_roll_pitch_deg=proxy)
                summary = dict(observation_id=record, offset_deg=offset,
                    source_global_realization_sha256=sha(bank/record/'realization-global.json'),
                    hm0_m=float(4*np.sqrt(np.sum(weights))), tm02_s=source_summary['tm02_s'],
                    body_mean_propagation_deg=float(np.rad2deg(np.angle(vector)) % 360),
                    body_mean_arrival_deg=float((np.rad2deg(np.angle(vector))+180) % 360),
                    first_resultant=float(abs(vector)/np.sum(weights)), yaw_global_rad=case['mean_box_yaw_global_rad'],
                    point_count=len(points), t0_fields_by_point=values[:, :, 0].tolist(),
                    t0_cg_slope_proxy_roll_pitch_deg=proxy[:, 0].tolist(),
                    sampled_time_max_abs_cg_slope_proxy_roll_pitch_deg=np.max(abs(proxy), axis=1).tolist(),
                    sampled_field_min=np.min(values, axis=(1, 2)).tolist(),
                    sampled_field_max=np.max(values, axis=(1, 2)).tolist(),
                    t0_plane_departure_by_point_m=departure[:, 0].tolist(),
                    t0_sampled_max_abs_plane_departure_m=float(np.max(abs(departure[:, 0]))),
                    sampled_time_footprint_max_abs_plane_departure_m=float(np.max(abs(departure))),
                    continuous_all_time_footprint_taylor_departure_bound_m=bound,
                    global_all_time_hessian_operator_norm_bound_per_m=hessian_bound,
                    direct_t0_field_check_max_error=float(np.max(abs(values[:, :, 0]-direct))),
                    interpretation='Prescribed incident surface only. Local tangent-plane slope angles do not predict hull rotations, restoring, forces, or initialized fluid/body state.')
                write(output/f'{prefix}.json', summary)
                manifest['cases'].append(summary)
        assert frozen == {str(path.resolve()): sha(path) for path in inputs+binaries}
        manifest.update(complete=True, elapsed_s=time.perf_counter()-started,
            output_sha256={path.relative_to(output).as_posix(): sha(path)
                           for path in output.rglob('*') if path.is_file() and path.name != 'manifest.json'})
        write(output/'manifest.json', manifest)
        print(json.dumps(dict(complete=True, cases=len(manifest['cases']), elapsed_s=manifest['elapsed_s'])))
    except BaseException as error:
        manifest.update(error=repr(error), elapsed_s=time.perf_counter()-started)
        write(output/'manifest.json', manifest)
        raise


if __name__ == '__main__':
    main()
