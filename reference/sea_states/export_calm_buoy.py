"""Export one frozen mild observation through the existing generator, no dynamics."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'

import argparse
import base64
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from topwater import sea_states


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    config = json.loads(args.config.read_text(encoding='utf-8-sig'))
    source = Path(config['source_library'])
    assert sha(source) == config['source_library_sha256']
    assert config['rescale_generated_amplitudes'] is False
    assert config['relative_heading_offset_deg'] == 0
    physical = json.loads((ROOT/config['physical_input']).read_text())
    catalog = json.loads((ROOT/config['catalog_provenance']).read_text())
    data = json.loads(source.read_text())
    library = data['sampling_library']
    assert library['version'] == 'marine-sea-library/1'
    region = library['regions'][config['region']]
    edges = np.asarray(region['frequency_edges_Hz'])
    bands, count = len(edges)-1, region['count']
    assert bands == 48 and region['record_bytes'] == 8+8*bands
    raw = gzip.decompress(base64.b64decode(region['rows_gzip_base64'], validate=True))
    assert len(raw) == count*region['record_bytes']
    rows = np.frombuffer(raw, dtype=np.dtype([('time', '<f8'), ('amplitude', '<f8', (bands,))]))
    index = config['source_row']
    stamp = datetime.fromtimestamp(float(rows['time'][index]), timezone.utc).isoformat().replace('+00:00', 'Z')
    assert stamp == config['time_utc']
    raw_moments = gzip.decompress(base64.b64decode(region['direction_moments_int16_gzip_base64'], validate=True))
    assert len(raw_moments) == count*bands*4*2
    packed = np.frombuffer(raw_moments, dtype='<i2').reshape(count, bands, 4)[index]
    assert not np.any(packed == region['direction_moment_missing'])
    source_files = [item for item in catalog['sources'] if '2023.' in item['file']]
    assert len(source_files) == 5
    provenance = dict(station_id='46086', source_library_sha256=sha(source),
        source_library_row=str(index), source_region=config['region'],
        source_files=json.dumps(source_files, separators=(',', ':')),
        direction_measurement='paired four Fourier FROM moments, source int16/32767',
        selection='actual mild observation; fixed seed42; no amplitude, phase or response fitting',
        coverage='all 48 observed usable bands .030-.495Hz; unmeasured tails unknown')
    state = sea_states.ObservedSeaState(config['observation_id'], stamp, edges,
        rows['amplitude'][index], packed/32767., provenance)
    assert abs(state.hm0_m-config['expected_hm0_m']) <= 1e-14
    assert abs(state.tm02_s-config['expected_tm02_s']) <= 1e-12
    vertices = np.asarray(physical['minimum_hull']['vertices_material_m'])
    beam = float(np.ptp(vertices[:, 1]))
    gravity = float(physical['water']['gravity_m_s2'])
    waves = sea_states.realize(state, beam_m=beam, seed=config['seed'],
        components_per_band=config['components_per_band'], gravity=gravity)
    text = sea_states.dumps_realization(waves)
    replay = sea_states.loads_realization(text, beam_m=beam)
    for name in ('amplitudes', 'angular_frequencies', 'directions', 'phases'):
        np.testing.assert_array_equal(getattr(waves, name), getattr(replay, name))
    per_band = config['components_per_band']
    variance = waves.amplitudes**2/2
    np.testing.assert_allclose(variance.reshape(bands, per_band).sum(axis=1),
        state.band_amplitudes_m**2/2, rtol=8e-15, atol=0)
    sea_states.validate_envelope(waves, beam_m=beam, gravity=gravity)
    vector = np.array([variance@np.cos(waves.directions), variance@np.sin(waves.directions)])
    resultant = float(np.linalg.norm(vector)/variance.sum())
    assert resultant > 64*np.finfo(float).eps
    propagation = float(np.arctan2(vector[1], vector[0]) % (2*np.pi))
    yaw = float((propagation+np.pi) % (2*np.pi))
    rotation = np.array([[np.cos(yaw), -np.sin(yaw), 0.],
                         [np.sin(yaw), np.cos(yaw), 0.], [0., 0., 1.]])
    cg = np.asarray(physical['world_cg_m'])
    translation = cg-rotation@cg
    k_global = waves.angular_frequencies[:, None]**2/gravity*np.column_stack((np.cos(waves.directions), np.sin(waves.directions)))
    direction_body = (waves.directions-yaw) % (2*np.pi)
    phase_body = waves.phases+k_global@translation[:2]
    xy = np.array([[0., 0.], [beam, -beam/2], [-beam, beam]])
    times = np.array([0., .37, 3.1])
    eta_global = np.sum(waves.amplitudes[:, None]*np.cos(k_global@xy.T-waves.angular_frequencies[:, None]*times+waves.phases[:, None]), axis=0)
    xy_body = (xy-translation[:2])@rotation[:2, :2]
    k_body = waves.angular_frequencies[:, None]**2/gravity*np.column_stack((np.cos(direction_body), np.sin(direction_body)))
    eta_body = np.sum(waves.amplitudes[:, None]*np.cos(k_body@xy_body.T-waves.angular_frequencies[:, None]*times+phase_body[:, None]), axis=0)
    tolerance = 2048*np.finfo(float).eps*max(1., float(waves.amplitudes.sum()))
    np.testing.assert_allclose(eta_body, eta_global, rtol=0, atol=tolerance)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/'realization-global.json').write_text(text, encoding='utf-8')
    write(args.output/'observation.json', dict(id=state.id, time_utc=stamp,
        frequency_edges_hz=edges.tolist(), band_amplitudes_m=state.band_amplitudes_m.tolist(),
        direction_moments_int16=packed.tolist(), provenance=provenance))
    np.savez(args.output/'components.npz', amplitudes=waves.amplitudes,
        angular_frequencies=waves.angular_frequencies, directions=waves.directions, phases=waves.phases,
        band_index=np.repeat(np.arange(bands), per_band), component_id=np.arange(bands*per_band),
        directions_body=direction_body, phases_body=phase_body, frequency_edges_hz=edges)
    write(args.output/'case.json', dict(observation_id=state.id, seed=waves.seed,
        mean_box_yaw_global_rad=yaw, relative_heading_offset_deg=0,
        heading_statistic='energy-weighted first circular moment of actual generated propagation, plus pi',
        global_from_boat_rotation=rotation.tolist(), global_from_boat_translation_m=translation.tolist(),
        boat_frame_waves=dict(amplitudes=waves.amplitudes.tolist(), angular_frequencies=waves.angular_frequencies.tolist(),
            directions=direction_body.tolist(), phases=phase_body.tolist()),
        physical_input=config['physical_input'], physical_input_sha256=sha(ROOT/config['physical_input']),
        startup_status=config['arrival_status']))
    summary = dict(complete=True, observation_id=state.id, time_utc=stamp, source_row=index,
        seed=waves.seed, components=bands*per_band, hm0_m=float(state.hm0_m), tm02_s=state.tm02_s,
        peak_density_band_period_s=state.peak_band_period_s, variance_m2=state.variance_m2,
        realized_variance_m2=float(variance.sum()), mean_propagation_rad=propagation,
        directional_first_resultant=resultant, bow_yaw_global_rad=yaw,
        frame_elevation_error_m=float(np.max(abs(eta_global-eta_body))),
        source_library_sha256=sha(source), startup_status=config['arrival_status'],
        elapsed_s=time.monotonic()-started,
        output_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()})
    write(args.output/'summary.json', summary)
    assert summary['elapsed_s'] <= config['maximum_elapsed_seconds']
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
