"""Export unchanged observed-sea generator candidates, never boat dynamics."""
import os
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
from topwater import sea_states


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, allow_nan=False, indent=2)+'\n', encoding='utf-8')


def mean_direction(weights, cosine, sine):
    total = float(np.sum(weights))
    vector = np.array([np.dot(weights, cosine), np.dot(weights, sine)])
    resultant = float(np.linalg.norm(vector)/total)
    if resultant <= 64*np.finfo(float).eps:
        raise ValueError('Mean directional resultant is numerically zero; heading is undefined')
    return float(np.arctan2(vector[1], vector[0]) % (2*np.pi)), resultant


def elevation(amplitude, omega, direction, phase, points, times, gravity):
    vectors = omega[:, None]**2/gravity*np.column_stack((np.cos(direction), np.sin(direction)))
    arguments = vectors@points.T-omega[:, None]*times+phase[:, None]
    return np.sum(amplitude[:, None]*np.cos(arguments), axis=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text(encoding='utf-8-sig'))
    physical_path = ROOT/config['physical_input']
    physical = json.loads(physical_path.read_text(encoding='utf-8-sig'))
    vertices = np.asarray(physical['minimum_hull']['vertices_material_m'])
    beam = float(np.ptp(vertices[:, 1]))
    gravity = float(physical['water']['gravity_m_s2'])
    cg = np.asarray(physical['world_cg_m'], dtype=float)
    assert config['rescale_generated_amplitudes'] is False
    assert config['default_relative_heading_offset_deg'] == 0
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError('Use a fresh output directory; accepted data must not be overwritten')
    output.mkdir(parents=True)
    source_paths = [Path(__file__).resolve(), config_path, physical_path]
    source_paths += [ROOT/'src/topwater'/name for name in
                     ('__init__.py', '_validation.py', 'inputs.py', 'sea_states.py', 'data/socal_sea_states.json')]
    frozen = {path.relative_to(ROOT).as_posix(): sha(path) for path in source_paths}
    for path in source_paths:
        archived = output/'source'/path.relative_to(ROOT)
        archived.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, archived)
    dependency_paths = list(Path(np.random.__file__).parent.glob('*.pyd'))
    dependency_paths += list((Path(np.__file__).parent/'_core').glob('*.pyd'))
    dependency_hashes = {str(path): sha(path) for path in dependency_paths}
    manifest = dict(scope='Generator data export only: no BEM, dynamics, motion fitting or nonlinear qualification',
        started_at=datetime.now(timezone.utc).isoformat(), complete=False, input_sha256=frozen,
        command=sys.argv, python=sys.version, python_executable_sha256=sha(sys.executable),
        dependency_versions={name: importlib.metadata.version(name) for name in ('numpy', 'scipy')},
        numpy_binary_sha256=dependency_hashes,
        generation=dict(seed=config['seed'], components_per_band=config['components_per_band']),
        beam_m=beam, gravity_m_s2=gravity,
        frequency_limit_hz=min(.6, float(np.sqrt(gravity/(4*np.pi*beam)))), records=[])
    write(output/'manifest.json', manifest)
    started = time.perf_counter()
    try:
        for record_id in config['observation_ids']:
            state = sea_states.observed(record_id)
            waves = sea_states.realize(state, beam_m=beam, gravity=gravity,
                seed=config['seed'], components_per_band=config['components_per_band'])
            folder = output/record_id
            folder.mkdir()
            serialized = sea_states.dumps_realization(waves)
            (folder/'realization-global.json').write_text(serialized, encoding='utf-8')
            replay = sea_states.loads_realization(serialized, beam_m=beam)
            for name in ('amplitudes', 'angular_frequencies', 'directions', 'phases'):
                np.testing.assert_array_equal(getattr(waves, name), getattr(replay, name))
            count = config['components_per_band']
            realized_band_variance = np.sum(waves.amplitudes.reshape(-1, count)**2/2, axis=1)
            observed_band_variance = state.band_amplitudes_m**2/2
            np.testing.assert_allclose(realized_band_variance, observed_band_variance, rtol=8e-15, atol=0)
            sea_states.validate_envelope(waves, beam_m=beam, gravity=gravity)
            component_energy = waves.amplitudes**2/2
            observed_theta, observed_r = mean_direction(observed_band_variance,
                -state.direction_moments[:, 1], -state.direction_moments[:, 0])
            theta, resultant = mean_direction(component_energy, np.cos(waves.directions), np.sin(waves.directions))
            summary = dict(observation_id=record_id, observation_time_utc=state.time_utc,
                hm0_m=float(state.hm0_m), tm02_s=state.tm02_s,
                peak_density_band_period_s=state.peak_band_period_s, components=len(waves.amplitudes),
                realized_variance_m2=float(np.sum(component_energy)), observed_variance_m2=state.variance_m2,
                band_variance_max_absolute_error_m2=float(np.max(abs(realized_band_variance-observed_band_variance))),
                observed_mean_propagation_rad=observed_theta, observed_first_resultant=observed_r,
                realized_mean_propagation_rad=theta, realized_first_resultant=resultant,
                heading_statistic='Energy-weighted first circular moment of the actual generated components; not the angular PDF mode',
                mean_direction_warning='A small resultant makes this mean ambiguous; report it, never replace the directional spread',
                frequency_range_hz=[float(waves.angular_frequencies.min()/(2*np.pi)), float(waves.angular_frequencies.max()/(2*np.pi))],
                generator_json_sha256=sha(folder/'realization-global.json'),
                startup_status='Stationary harmonic realization only; no calm-start packet or nonlinear fluid initialization supplied',
                selection_status=config['selection_status'], heading_candidates=[])
            global_points = np.array([[0., 0.], [beam, -beam/2], [-beam, beam]])
            times = np.array([0., .37, 3.1])
            global_eta = elevation(waves.amplitudes, waves.angular_frequencies,
                waves.directions, waves.phases, global_points, times, gravity)
            for offset in config['relative_heading_offsets_deg']:
                yaw = float((theta+np.pi+np.deg2rad(offset)) % (2*np.pi))
                cosine, sine = np.cos(yaw), np.sin(yaw)
                rotation = np.array([[cosine, -sine, 0.], [sine, cosine, 0.], [0., 0., 1.]])
                translation = cg-rotation@cg
                directions = (waves.directions-yaw) % (2*np.pi)
                k_global = waves.angular_frequencies[:, None]**2/gravity*np.column_stack(
                    (np.cos(waves.directions), np.sin(waves.directions)))
                phases = waves.phases+k_global@translation[:2]
                boat_points = (global_points-translation[:2])@rotation[:2, :2]
                boat_eta = elevation(waves.amplitudes, waves.angular_frequencies,
                    directions, phases, boat_points, times, gravity)
                tolerance = 2048*np.finfo(float).eps*max(1., float(np.sum(waves.amplitudes)))
                np.testing.assert_allclose(boat_eta, global_eta, rtol=0, atol=tolerance)
                case = dict(format_version=1, purpose=config['purpose'], observation_id=record_id,
                    global_realization='realization-global.json', relative_heading_offset_deg=offset,
                    default_heading_case=offset == config['default_relative_heading_offset_deg'],
                    forward_convention=config['forward_convention'], mean_box_yaw_global_rad=yaw,
                    global_axes='east/north/up; propagation angles CCW from east',
                    boat_aligned_axes='initial forward/port/up; propagation angles CCW from forward',
                    global_from_boat_rotation=rotation.tolist(), global_from_boat_translation_m=translation.tolist(),
                    global_phase_origin_xy_m=[0., 0.],
                    boat_frame_waves=dict(amplitudes=waves.amplitudes.tolist(),
                        angular_frequencies=waves.angular_frequencies.tolist(), directions=directions.tolist(), phases=phases.tolist()),
                    frame_elevation_closure_max_error_m=float(np.max(abs(boat_eta-global_eta))),
                    physical_input=config['physical_input'], physical_input_sha256=frozen[physical_path.relative_to(ROOT).as_posix()],
                    startup_status=summary['startup_status'], selection_status=config['selection_status'])
                file_name = f'heading-offset-{offset:g}deg.json'
                write(folder/file_name, case)
                summary['heading_candidates'].append(dict(file=file_name, offset_deg=offset, yaw_global_rad=yaw))
            write(folder/'summary.json', summary)
            manifest['records'].append(summary)
        assert frozen == {path.relative_to(ROOT).as_posix(): sha(path) for path in source_paths}
        assert dependency_hashes == {str(path): sha(path) for path in dependency_paths}
        manifest.update(complete=True, elapsed_s=time.perf_counter()-started,
            checks='Every measured band variance preserved; exact JSON replay; full incident support valid; global/body phase-coordinate closure',
            output_sha256={path.relative_to(output).as_posix(): sha(path)
                           for path in output.rglob('*.json') if path.name != 'manifest.json'})
        write(output/'manifest.json', manifest)
        print(json.dumps(dict(complete=True, records=len(manifest['records']),
                              elapsed_s=manifest['elapsed_s'], output=str(output))))
    except BaseException as error:
        manifest.update(error=repr(error), elapsed_s=time.perf_counter()-started)
        write(output/'manifest.json', manifest)
        raise


if __name__ == '__main__':
    main()
