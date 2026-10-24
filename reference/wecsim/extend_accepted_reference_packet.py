"""Append only future samples of the frozen accepted REF001 field; no new packet design."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.io import loadmat, savemat
from threadpoolctl import threadpool_limits
from incident_packet import IncidentPacket


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    control_path = Path('reference/wecsim/reference-extension-controls.json')
    control = json.loads(control_path.read_text())
    original = Path(control['original_packet'])
    prior = json.loads((original/'diagnostic.json').read_text())
    region_path = Path('artifacts/reference/minimum-20261024/packet-full-planform/diagnostic.json')
    region = json.loads(region_path.read_text())
    assert prior['accepted'] and region['accepted']
    assert sha(original/'packet.npz') == region['frozen_packet_sha256']
    physical_path = Path('reference/capytaine/canonical-boat.json')
    old_controls = Path('reference/wecsim/packet-controls.json')
    prior_hashes = {key.replace('\\','/'):value for key,value in prior['input_sha256'].items()}
    assert sha(physical_path) == prior_hashes[physical_path.as_posix()]
    assert sha(old_controls) == prior_hashes[old_controls.as_posix()]
    field = IncidentPacket(json.loads(physical_path.read_text()),json.loads(old_controls.read_text()))
    packet = np.load(original/'packet.npz')
    for name, actual in [('w',field.w),('k',field.k),('c',field.c),('weights',field.weights),
                         ('direction',field.theta),('carrier_omega',field.w0),
                         ('carrier_amplitude',field.a),('focus_time',field.tc)]:
        assert np.array_equal(actual,packet[name]), name
    dt = control['dt_s']
    assert dt == .02 and control['forcing_time_s'] == [-20.,80.]
    assert control['output_time_s'] == [0.,60.] and control['radiation_memory_s'] == 10.
    args.output.mkdir(parents=True,exist_ok=False)
    future_time = np.arange(1501,4001)*dt
    future_eta = np.array([field.elevation([[0.,0.]],t)[0] for t in future_time])
    old_time = packet['times']
    old_eta = packet['elevation']
    assert old_time.shape == (2501,) and old_eta.shape == (2501,3)
    assert abs(old_time[0]+20) < 1e-12 and abs(old_time[-1]-30) < 1e-12
    tf = np.r_[old_time,future_time]
    eta = np.vstack([old_eta,future_eta])
    assert np.max(abs(np.diff(tf)-dt)) < 1e-12
    for j in range(3):
        old_path = original/f'wave_{j+1}.txt'
        old_text = old_path.read_bytes()
        assert old_text.endswith(b'\n')
        parsed = np.loadtxt(old_path)
        assert np.array_equal(parsed[:,0],old_time) and np.array_equal(parsed[:,1],old_eta[:,j])
        dest = args.output/f'wave_{j+1}.txt'
        dest.write_bytes(old_text)
        with dest.open('ab') as stream:
            np.savetxt(stream,np.c_[future_time,future_eta[:,j]],fmt='%.16e')
        assert dest.read_bytes().startswith(old_text)
    np.savez(args.output/'packet.npz',w=field.w,k=field.k,c=field.c,weights=field.weights,
        direction=field.theta,times=tf,elevation=eta,carrier_omega=field.w0,
        carrier_amplitude=field.a,focus_time=field.tc)
    bemio_path = Path('artifacts/reference/minimum-20261024/consumer/bemio.mat')
    h = loadmat(bemio_path,simplify_cells=True)['h']
    physical = json.loads(physical_path.read_text())
    kernel = h['ex_K']*physical['water']['density_kg_m3']*physical['water']['gravity_m_s2']
    assert np.max(abs(np.asarray(h['ex_t']).reshape(-1)-np.arange(-1000,1001)*dt)) < 1e-12
    force = sum(np.stack([np.convolve(eta[:,j],kernel[q,j],mode='valid')*dt
        for q in range(6)],axis=1) for j in range(3))
    assert force.shape == (3001,6)
    old_force = np.load(original/'startup-forces.npz')['force_convolution']
    prefix_error = np.max(abs(force[:501]-old_force),axis=0)
    assert np.array_equal(force[:501],old_force),'Original retained force prefix must be bit-identical'
    expected_path = Path('artifacts/reference/minimum-20261024/expected-retained-force.mat')
    expected = loadmat(expected_path,simplify_cells=True)
    tolerance = 1e-9+1e-10*np.max(abs(expected['force']),axis=0)
    native_prefix_error = np.max(abs(force[:501]-expected['force']),axis=0)
    assert np.all(native_prefix_error <= tolerance)
    limits = np.asarray(prior['force_limits'])
    assert np.all(abs(force[0]) <= limits)
    assert np.all(np.asarray(prior['initial_spectral_force']) <= limits)
    assert np.max(abs(abs(force[0])-prior['initial_actual_force'])) < 1e-12
    elapsed = time.monotonic()-start
    assert elapsed <= control['maximum_elapsed_seconds']
    to = np.arange(3001)*dt
    savemat(args.output/'body-force.mat',dict(time=to[:,None],force=force))
    report = dict(complete=True,accepted=True,unchanged_physical_field=True,
        inherited_fullplanform_report=region,initial_retained_force=force[0].tolist(),
        inherited_initial_direct_force_abs=prior['initial_spectral_force'],force_limits=limits.tolist(),
        old_wave_prefix_bytes_identical=True,old_spectral_arrays_identical=True,
        first501_retained_force_bit_identical=True,retained_prefix_error=prefix_error.tolist(),
        expected_native_prefix_error=native_prefix_error.tolist(),
        inherited_arrival_fraction=prior['arrival_fraction_by_component'],elapsed_s=elapsed,
        operator='ORIGINAL REF001 H5/CIC10, unqualified accuracy; duration only changed',
        input_sha256={p.as_posix():sha(p) for p in [control_path,physical_path,old_controls,
            original/'packet.npz',original/'diagnostic.json',region_path,bemio_path,expected_path]})
    (args.output/'diagnostic.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    return 0


if __name__ == '__main__':
    with threadpool_limits(1):
        raise SystemExit(main())
