"""Two-phase, same-state extension sensitivity of the linear radiation probe.

The alternative fill is verification-only. It minimizes squared discrete
Laplacians on hidden nodes AND their exterior neighbors. Physical exterior
nodal values are held fixed; finite Fourier interpolation can still leak the
change between nodes. No alternate fill enters the evolving solver.
"""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import shutil
import sys
from time import perf_counter

import numpy as np
from scipy.spatial import ConvexHull

from topwater import SimulationInputs
from topwater.constants import GRAVITY


def biharmonic_fill(size, hidden):
    hidden = np.asarray(hidden)
    lookup = {int(index): col for col, index in enumerate(hidden)}

    def neighbors(index):
        row, col = divmod(int(index), size)
        return [row*size+(col-1) % size, row*size+(col+1) % size,
                ((row-1) % size)*size+col, ((row+1) % size)*size+col]

    rows = sorted(set(hidden) | {j for index in hidden for j in neighbors(index)})
    support = set(rows) | {j for index in rows for j in neighbors(index)}
    boundary = np.array(sorted(support-set(hidden)))
    external = {int(index): col for col, index in enumerate(boundary)}
    a, b = np.zeros((len(rows), len(hidden))), np.zeros((len(rows), len(boundary)))
    for row, index in enumerate(rows):
        for column, coefficient in [(index, 4), *[(j, -1) for j in neighbors(index)]]:
            if column in lookup:
                a[row, lookup[column]] += coefficient
            else:
                b[row, external[column]] += coefficient
    weights, _, rank, singular = np.linalg.lstsq(a, -b, rcond=None)
    if rank != len(hidden):
        raise ValueError('hidden least-squares extension must have full column rank')
    return boundary, weights, dict(rows=len(rows), hidden_count=len(hidden),
        exterior_support_count=len(boundary), rank=int(rank),
        singular_value_ratio=float(singular[-1]/singular[0]),
        stationarity_relative_residual=float(np.linalg.norm(a.T@(a@weights+b))/np.linalg.norm(a.T@b)))


def exterior_probes(system):
    points = system.panels.triangles.reshape(-1, 3)
    footprint = ConvexHull(points[points[:, 2] == 0, :2])
    low = np.floor((footprint.points.min(axis=0)-system.origin)/system.spacing).astype(int)-2
    high = np.ceil((footprint.points.max(axis=0)-system.origin)/system.spacing).astype(int)+2
    x, y = np.meshgrid(system.origin[0]+(np.arange(low[0], high[0]+1)+.5)*system.spacing,
                        system.origin[1]+(np.arange(low[1], high[1]+1)+.5)*system.spacing)
    xy = np.column_stack((x.ravel(), y.ravel()))
    plane = footprint.equations
    return xy[np.any(xy@plane[:, :2].T+plane[:, 2] > 0, axis=1)]


def sample_spectrum(system, points, spectra):
    cp = system.cp
    kx, ky = cp.meshgrid(2*np.pi*cp.fft.fftfreq(system.size, d=system.spacing),
                         2*np.pi*cp.fft.fftfreq(system.size, d=system.spacing))
    values = []
    for first in range(0, len(points), 8):
        xy = cp.asarray(points[first:first+8]-system.origin)
        phase = xy[:, 0, None]*kx.ravel()+xy[:, 1, None]*ky.ravel()
        weight = cp.exp(1j*phase)/(system.size**2)
        values.append(system.apply_columns(weight, spectra).real)
    return cp.concatenate(values)


def snapshot(system, state, time, omega, ramp, boundary, weights, points, reference, drive):
    cp = system.cp
    _, force, _, _ = system.rhs(time, state, omega, ramp)
    rate = -GRAVITY*state[0]-system.sponge[:, None]*state[1]
    original = system.extend(rate).reshape(-1, 2)
    alternate = original.copy()
    alternate[system.device_hidden] = cp.asarray(weights)@alternate[cp.asarray(boundary)]
    delta = alternate-original
    spectrum = cp.fft.fft2(delta.reshape(system.size, system.size, 2), axes=(0, 1)).reshape(-1, 2)
    force_change = -system.apply_columns(system.impulse_trace, spectrum).real
    original_spectrum = cp.fft.fft2(original.reshape(system.size, system.size, 2), axes=(0, 1)).reshape(-1, 2)
    leakage = sample_spectrum(system, points, spectrum)
    original_trace = sample_spectrum(system, points, original_spectrum)
    reactive = omega*np.asarray(reference['added_mass'])*system.amplitudes
    speed, acceleration = drive(time, omega, ramp)
    return dict(time_seconds=time, velocity_multiplier=speed, acceleration_multiplier=acceleration,
        current_force=cp.asnumpy(force).tolist(), alternative_force=cp.asnumpy(force+force_change).tolist(),
        force_change=cp.asnumpy(force_change).tolist(), capytaine_reactive_force_amplitude=reactive.tolist(),
        hidden_rate_change_max=cp.asnumpy(cp.max(cp.abs(delta), axis=0)).tolist(),
        exterior_nodal_change_max=cp.asnumpy(cp.max(cp.abs(delta[system.device_exterior]), axis=0)).tolist(),
        offgrid_exterior_trace_change_max=cp.asnumpy(cp.max(cp.abs(leakage), axis=0)).tolist(),
        offgrid_exterior_trace_change_rms=cp.asnumpy(cp.sqrt(cp.mean(leakage**2, axis=0))).tolist(),
        original_offgrid_trace_max=cp.asnumpy(cp.max(cp.abs(original_trace), axis=0)).tolist())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--refinement', type=int, default=1)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('new audit output required')
    report = json.loads(args.report.read_text(encoding='utf-8'))
    inputs = SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8'))
    if report['input_sha256'] != sha256(inputs.to_json().encode()).hexdigest() or not report['complete']:
        raise ValueError('completed report with identical physical inputs required')
    arrays = args.report.with_suffix('.npz')
    if report['arrays_sha256'] != sha256(arrays.read_bytes()).hexdigest():
        raise ValueError('owned state arrays differ from report')
    reference_path = args.report.with_suffix('.source')/'reference.json'
    if sha256(reference_path.read_bytes()).hexdigest() != report['reference_sha256']:
        raise ValueError('owned independent reference differs')
    references = json.loads(reference_path.read_text(encoding='utf-8'))
    # Continue the exact owned driver. All its imported producers must still
    # match, rather than treating one reconstructed force as equation parity.
    driver_path = Path('tests/experiments/free_surface/boat_radiation_time.py')
    parent_archive = args.report.with_suffix('.source')
    for name, digest in report['producer_sha256'].items():
        source = parent_archive/name if Path(name) == driver_path else Path(name)
        if sha256(source.read_bytes()).hexdigest() != digest:
            raise ValueError(f'owned producer mismatch: {name}')
    spec = importlib.util.spec_from_file_location(
        __package__+'._owned_radiation', parent_archive/driver_path)
    owned = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = owned
    spec.loader.exec_module(owned)
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    sources = [Path(__file__), *[Path(name) for name in report['producer_sha256']]]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        if relative == driver_path:
            source = parent_archive/driver_path
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    start = perf_counter()
    system = owned.Radiation(inputs, report['grid_size'], report['box_length_m'], args.refinement,
                       report['sponge_strength_per_second'])
    if (report['driven_modes'] != ['Heave', 'Roll']
            or report['force_rows'] != ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']
            or report['integrator'] != 'fixed-step RK4'
            or not np.array_equal(report['velocity_amplitudes'], system.amplitudes)):
        raise ValueError('owned mode ordering, amplitudes and integrator must match')
    if sha256(system.panels.triangles.tobytes()).hexdigest() != report['panel_sha256']:
        raise ValueError('reconstructed physical panels differ')
    data = np.load(arrays, allow_pickle=False)
    np.testing.assert_array_equal(system.exterior, data['exterior'])
    np.testing.assert_array_equal(system.hidden, data['hidden'])
    cp = system.cp
    state = cp.asarray(data['final_state'])
    time, step, ramp = report['duration_seconds'], report['time_step_seconds'], report['ramp_seconds']
    omega = 2*np.pi*report['frequency_hz']
    initial_force = cp.asnumpy(system.rhs(time, state, omega, ramp)[1])
    np.testing.assert_allclose(initial_force, data['forces'][-1], atol=1e-10, rtol=1e-11)
    reference = next(row['reference'] for row in references['rows']
        if row['panel_sha256'] == report['panel_sha256']
        and row['reference']['angular_frequency'] != 'infinite'
        and np.isclose(row['reference']['angular_frequency'], omega, rtol=0, atol=1e-14))
    boundary, weights, diagnostics = biharmonic_fill(system.size, system.hidden)
    probes = exterior_probes(system)
    rows = [snapshot(system, state, time, omega, ramp, boundary, weights, probes, reference, owned.drive)]
    count = round((.25/report['frequency_hz'])/step)
    if not np.isclose(count*step, .25/report['frequency_hz']):
        raise ValueError('quarter period must lie on the owned time grid')
    for index in range(count):
        t = time+index*step
        k1 = system.rhs(t, state, omega, ramp)[0]
        k2 = system.rhs(t+step/2, state+step*k1/2, omega, ramp)[0]
        k3 = system.rhs(t+step/2, state+step*k2/2, omega, ramp)[0]
        k4 = system.rhs(t+step, state+step*k3, omega, ramp)[0]
        state += step*(k1+2*k2+2*k3+k4)/6
    rows.append(snapshot(system, state, time+count*step, omega, ramp, boundary, weights, probes, reference, owned.drive))
    cp.cuda.get_current_stream().synchronize()
    result = dict(complete=True, qualified=False,
        scope='Same-state extension sensitivity including off-grid Fourier leakage; no alternative runtime fill or trajectory qualification',
        parent_report_sha256=sha256(args.report.read_bytes()).hexdigest(),
        parent_arrays_sha256=report['arrays_sha256'], producer_sha256=hashes,
        numerical_producers='Archived parent driver; all imported producers must match parent hashes exactly',
        units=dict(force_rows=['N']*3+['Nm']*3, potential_rate_trace='m^2/s^2'),
        initial_force_reconstruction_max_error=float(np.max(np.abs(initial_force-data['forces'][-1]))),
        continuation='Same harmonic-fill discrete equations, owned absolute time/step/drive; alternate fill is used only to audit instantaneous pressure',
        biharmonic=diagnostics, offgrid_probe_count=len(probes), elapsed_seconds=perf_counter()-start,
        rows=rows, limits=['Exterior nodal equality is not continuum trace equality; finite Fourier changes leak between nodes',
                          'Two phases screen instantaneous force sensitivity, not history feedback or total trajectory error'])
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
