"""Coordinator-requested reductions of completed stored arrays; no solver call."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extrema(values, omega):
    low, high = int(np.argmin(values)), int(np.argmax(values))
    return dict(min=float(values[low]), max=float(values[high]),
                minimum_component_index_1based=low+1, maximum_component_index_1based=high+1,
                minimum_omega_rad_s=float(omega[low]), maximum_omega_rad_s=float(omega[high]),
                minimum_frequency_hz=float(omega[low]/(2*np.pi)),
                maximum_frequency_hz=float(omega[high]/(2*np.pi)),
                negative_count=int(np.count_nonzero(values < 0)),
                zero_count=int(np.count_nonzero(values == 0)), positive_count=int(np.count_nonzero(values > 0)))


base = Path('artifacts/reference/generated-bank-consumer-20261024')
out = base/'signed-band-report.json'
assert not out.exists()
rad_path = Path('artifacts/reference/generated-bank-native-20261024/attempt1/radiation.mat')
rad = loadmat(rad_path)
records = []
for path in sorted((base/'compare-attempt1').glob('*-comparison.mat')):
    data = loadmat(path)
    omega = data['omega_rad_s'].ravel()
    difference = data['A_native_interpolated'][:, 2, 2]-data['A_infinite_time_finiteband'][:, 2, 2]
    records.append(dict(id=path.name.removesuffix('-comparison.mat'), source_sha256=sha(path),
        signed_A33_native_minus_band_kg=extrema(difference, omega),
        interpolated_native_B33_N_s_m=extrema(data['B_native_interpolated'][:, 2, 2], omega)))
report = dict(scope='Only reductions of completed immutable stored arrays; no BEM, interpolation change, ODE or body run',
    records=records,
    raw_native_B33_N_s_m=extrema(rad['radiation_damping'][:, 2, 2], rad['omega_rad_s'].ravel()),
    radiation_source_sha256=sha(rad_path), source_sha256=sha(__file__),
    sign_rule='For passive B33, omitted upper nu>12 tail is positive for incident omega<12; omitted lower nu<.05 tail is negative',
    limitation='Stored native B33 sign does not prove passivity of uncomputed tails or all coupled modes; attribution unresolved')
out.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report))
