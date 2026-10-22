"""Import nine joint observed SoCal spectra from a frozen marine-gyro export.

Offline data preparation only. No downloader or marine-gyro runtime dependency.
Verify all frozen NOAA source hashes and selected native energy/direction rows.
"""
import argparse
import base64
from datetime import datetime, timezone
import gzip
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def source_rows(path, epochs):
    lines = gzip.decompress(path.read_bytes()).decode('ascii').splitlines()
    frequencies = np.asarray(list(map(float, lines[0].split()[5:])))
    rows = {}
    for line in lines[1:]:
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.split()
        epoch = int(datetime(*map(int, parts[:5]), tzinfo=timezone.utc).timestamp())
        if epoch in epochs:
            if epoch in rows:
                raise ValueError('duplicate native observation')
            rows[epoch] = np.asarray(list(map(float, parts[5:])))
    return frequencies, rows


def build(library_path, source_manifest, output, native_control):
    exported = json.loads(library_path.read_text())['sampling_library']
    config = json.loads(source_manifest.read_text())
    metadata = exported['regions']['california']
    if exported['version'] != 'marine-sea-library/1' or metadata.get('directional_unknown_hours') != 0:
        raise ValueError('complete current directional marine sea library required')
    edges = np.asarray(metadata['frequency_edges_Hz'])
    count, bands = int(metadata['count']), len(edges)-1
    if metadata['record_bytes'] != 8+8*bands:
        raise ValueError('exported row size differs')
    dtype = np.dtype([('time', '<f8'), ('amplitudes', '<f8', (bands,))])
    raw = gzip.decompress(base64.b64decode(metadata['rows_gzip_base64'], validate=True))
    if len(raw) != count*dtype.itemsize:
        raise ValueError('exported row count differs')
    rows = np.frombuffer(raw, dtype=dtype)
    times, amplitudes = rows['time'], rows['amplitudes']
    packed = gzip.decompress(base64.b64decode(metadata['direction_moments_int16_gzip_base64'], validate=True))
    moments = np.frombuffer(packed, dtype='<i2').reshape(count, bands, 4)
    if (not np.all(np.isfinite(amplitudes)) or np.any(amplitudes < 0)
            or np.any(moments == -32768) or len(np.unique(times)) != count):
        raise ValueError('finite unique observations with known paired directions required')
    variance = amplitudes**2/2
    m0 = variance.sum(axis=1)
    mean_f2 = (edges[1:]**3-edges[:-1]**3)/(3*np.diff(edges))
    features = np.column_stack((4*np.sqrt(m0), np.sqrt(m0/(variance@mean_f2))))
    quantiles = np.quantile(features, [.1, .25, .5, .75, .9], axis=0)
    scale = quantiles[3]-quantiles[1]
    if not np.all(scale > 0) or not np.all(np.isfinite(features)):
        raise ValueError('nondegenerate observed Hm0/Tm02 population required')
    selection, used = [], set()
    # Median first; remaining target order makes duplicate avoidance reproducible.
    for ih in (2, 0, 4):
        for it in (2, 0, 4):
            target = quantiles[[ih, it], [0, 1]]
            distance = np.sum(((features-target)/scale)**2, axis=1)
            order = np.lexsort((times, distance))
            index = next(int(i) for i in order if int(i) not in used)
            used.add(index)
            selection.append(dict(index=index, hm0_quantile={0: .1, 2: .5, 4: .9}[ih],
                tm02_quantile={0: .1, 2: .5, 4: .9}[it], distance_squared=float(distance[index])))
    station = next(s for s in config['stations'] if s['id'] == '46086')
    source_root = source_manifest.parent
    sources = {}
    directions = {source['energy_file']: source['variables'] for source in station['population_directional']}
    for energy in station['population_spectra']:
        for source in [energy, *directions[energy['file']].values()]:
            path = source_root/source['file']
            if not path.resolve().is_relative_to(source_root.resolve()) or digest(path) != source['sha256']:
                raise ValueError('frozen NOAA source checksum/path failed: '+source['file'])
            sources[source['file']] = {key: source[key] for key in ('file', 'url', 'sha256')}
    selected_epochs = {int(times[item['index']]) for item in selection}
    native, native_frequency = {}, None
    for source in sources.values():
        frequency, records = source_rows(source_root/source['file'], selected_epochs)
        if native_frequency is not None and not np.array_equal(native_frequency, frequency):
            raise ValueError('native NOAA frequency grids differ')
        native_frequency = frequency
        native[source['file']] = records
    expected = np.r_[.02, np.arange(13)*.005+.0325, np.arange(26)*.01+.1, np.arange(7)*.02+.365]
    if not np.allclose(native_frequency, expected, rtol=0, atol=1e-6):
        raise ValueError('unexpected NOAA spectral grid')
    widths = np.r_[.005, np.full(13, .005), np.full(26, .01), np.full(7, .02)]
    frequency, widths = native_frequency[1:], widths[1:]  # documented .020-Hz noise band excluded
    overlap = np.maximum(0, np.minimum((frequency+widths/2)[:, None], edges[None, 1:])
                         -np.maximum((frequency-widths/2)[:, None], edges[None, :-1]))
    if not np.allclose(overlap.sum(axis=0), np.diff(edges), rtol=0, atol=1e-6):
        raise ValueError('incomplete measured band coverage')
    records, control = [], None
    for item in selection:
        index = item.pop('index')
        epoch = int(times[index])
        stamp = datetime.fromtimestamp(epoch, timezone.utc)
        energy_source = next(source for source in station['population_spectra'] if epoch in native[source['file']])
        direction_sources = directions[energy_source['file']]
        density = native[energy_source['file']][epoch][1:]
        if np.any((density < 0) | (density >= 999)):
            raise ValueError('invalid selected native energy')
        original = np.sqrt(2*density@overlap)
        if not np.allclose(original, amplitudes[index], rtol=1e-13, atol=1e-14):
            raise ValueError('exported observation differs from frozen native energy')
        alpha1, alpha2, r1, r2 = [native[direction_sources[key]['file']][epoch][1:]
                                 for key in ('alpha1', 'alpha2', 'r1', 'r2')]
        if not (np.all((alpha1 >= 0)&(alpha1 <= 360)) and np.all((alpha2 >= 0)&(alpha2 <= 360))
                and np.all((r1 >= 0)&(r1 <= 100)) and np.all((r2 >= 0)&(r2 <= 100))):
            raise ValueError('invalid selected native directional observation')
        alpha1, alpha2 = np.deg2rad(alpha1), np.deg2rad(alpha2)
        native_moments = np.column_stack((.01*r1*np.cos(alpha1), .01*r1*np.sin(alpha1),
                                         .01*r2*np.cos(2*alpha2), .01*r2*np.sin(2*alpha2)))
        numerator = np.einsum('n,nk,nb->bk', density, native_moments, overlap, optimize=True)
        denominator = density@overlap
        rebinned = np.divide(numerator, denominator[:, None], out=np.zeros_like(numerator),
                             where=denominator[:, None] > 0)
        discrepancy = int(np.max(abs(np.rint(np.clip(rebinned, -1, 1)*32767)-moments[index])))
        if discrepancy > 1:
            raise ValueError('exported paired directions differ from frozen native observation')
        case_id = 'ndbc-46086-'+stamp.strftime('%Y%m%dT%H%MZ')
        files = [energy_source['file'], *[direction_sources[k]['file'] for k in ('alpha1', 'alpha2', 'r1', 'r2')]]
        record = dict(id=case_id, time_utc=stamp.isoformat(timespec='seconds').replace('+00:00', 'Z'),
            band_amplitudes_m=amplitudes[index].tolist(), direction_moments_int16=moments[index].tolist(),
            source_files=files, native_direction_verification_max_int16_units=discrepancy)
        records.append(record)
        item['record_id'] = case_id
        if control is None:
            control = dict(id=case_id, time_utc=record['time_utc'], energy_source=sources[energy_source['file']],
                frequency_hz=frequency.tolist(), bandwidth_hz=widths.tolist(), density_m2_per_hz=density.tolist(),
                direction_moments=native_moments.tolist())
    result = dict(version='topwater-observed-seas/1',
        attribution='NOAA National Data Buoy Center, station 46086, San Clemente Basin; frozen marine-gyro export.',
        station={key: station[key] for key in ('id', 'name', 'latitude_deg_north', 'longitude_deg_east', 'water_depth_m', 'station_url', 'scope')},
        source_library_sha256=digest(library_path), source_manifest_sha256=digest(source_manifest),
        builder_sha256=digest(Path(__file__)), source_retrieved_date_utc=config['retrieved_date_utc'],
        source_energy_conversion=dict(density_kg_m3=exported['density_kg_m3'],
                                      gravity_m_s2=exported['gravity_m_s2']),
        sources=list(sources.values()),
        coverage='0.03–0.495 Hz measured; unmeasured tails unknown; .020-Hz noise band excluded.',
        units=dict(amplitudes='cosine crest metres; variance=a^2/2', frequency='Hz',
                   direction_moments='int16/32767; (a1,b1,a2,b2), FROM clockwise from true north'),
        frequency_edges_hz=edges.tolist(), baseline_id=records[0]['id'],
        selection=dict(population_hours=count, population_years=config['years'],
            method='Nearest joint Hm0/Tm02 target in observed IQR units; earliest UTC tie; median first, '
                   'then targets in Hm0/Tm02 quantile order .5,.1,.9; previously selected hours excluded.',
            descriptors='Hm0=4sqrt(m0); Tm02=sqrt(m0/m2); uniform density within each retained interval.',
            population_quantiles=[dict(probability=p, hm0_m=float(row[0]), tm02_s=float(row[1]))
                                  for p, row in zip((.1, .25, .5, .75, .9), quantiles)],
            targets=selection, note='Nine coverage examples, not population weights or independent marginal draws.'), records=records)
    for path, value in ((output, result), (native_control, control)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('w', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(value, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(baseline_id=result['baseline_id'], records=len(records), population_hours=count)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('library', 'source-manifest', 'output', 'native-control'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    build(args.library, args.source_manifest, args.output, args.native_control)
