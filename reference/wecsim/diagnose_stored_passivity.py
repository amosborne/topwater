"""Allocated small stored-matrix reductions only; no ODE, pole or solver study."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.spatial.transform import Rotation


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def symmetric_energy(matrices, coordinate):
    values, vectors = np.linalg.eigh((matrices+np.swapaxes(matrices, -1, -2))/2)
    i = int(np.argmin(values[:, 0]))
    return dict(minimum_eigenvalue=float(values[i, 0]), coordinate_at_minimum=float(coordinate[i]),
                lowest_eigenvalue_range=[float(values[:, 0].min()), float(values[:, 0].max())],
                negative_minimum_count=int(np.count_nonzero(values[:, 0] < 0)),
                eigenvector_scaled_coordinates=vectors[i, :, 0].tolist(),
                maximum_absolute_eigenvalue=float(np.max(abs(values))))


base = Path('artifacts/reference/generated-bank-consumer-20261024')
case = base/'case-jan90-attempt1'
out = case/'stored-passivity-diagnosis.json'
assert not out.exists()
npath = Path('artifacts/reference/generated-bank-native-20261024/attempt1/native.mat')
rpath = npath.with_name('radiation.mat')
comparison_path = base/'jan90-attempt1/ndbc-46086-20250105T1140Z-comparison.mat'
n, rad, comparison = loadmat(npath), loadmat(rpath), loadmat(comparison_path)
motion = loadmat(case/'motion.mat', simplify_cells=True)
preflight = json.loads((case/'preflight.json').read_text())
D = np.diag([2.7432]*3+[1.]*3)
Di = np.linalg.inv(D)
B = rad['radiation_damping']
BK = comparison['B_kernel40']
w = rad['omega_rad_s'].ravel()
omega = comparison['omega_rad_s'].ravel()
M = n['M']
Ainf = n['Ainf']
solid = np.array(preflight['solid_local_matrix'])
Ar = np.array(preflight['residual_added_mass_world_matrix'])
np.testing.assert_allclose(solid+Ar, M+Ainf, atol=1e-8, rtol=0)
time = np.asarray(motion['time']).ravel()
angles = np.asarray(motion['position_world'])[:, 3:6]
rotation = Rotation.from_euler('xyz', angles).as_matrix()
radius = np.zeros(len(time))
norm = np.zeros(len(time))
real_min = np.zeros(len(time))
real_max = np.zeros(len(time))
selected = []
for i, R in enumerate(rotation):
    S = solid.copy()
    S[3:6, 3:6] = R@solid[3:6, 3:6]@R.T
    C = Di@np.linalg.solve(S, Ar)@D
    eigen = np.linalg.eigvals(C)
    radius[i] = np.max(abs(eigen))
    norm[i] = np.linalg.norm(C, 2)
    real_min[i], real_max[i] = eigen.real.min(), eigen.real.max()
    if i in [0, 1, 50, 500, 1000, 1500, 1700, 1732, 1750, 2000, 3000]:
        selected.append(dict(time_s=float(time[i]), eigenvalue_real=eigen.real.tolist(),
                             eigenvalue_imag=eigen.imag.tolist(), spectral_radius=float(radius[i]),
                             scaled_singular_norm=float(norm[i])))
minimum_mass = symmetric_energy(np.array([D@(M+Ainf)@D]), np.array([0.]))
report = dict(scope='Only stored matrix eigen/SVD reductions; no ODE/pole/simulation/control executed',
    coordinate_scaling='q=D*qhat, D=diag(beam,beam,beam,1,1,1); energy matrices D.T*M*D and D.T*B*D',
    raw_native_B_symmetric_energy=symmetric_energy(D@B@D, w),
    incident_band_native_B_symmetric_energy=symmetric_energy(D@comparison['B_native_interpolated']@D, omega),
    incident_band_kernel40_B_symmetric_energy=symmetric_energy(D@BK@D, omega),
    parameter_M_plus_Ainf_symmetric_energy=minimum_mass,
    dry_M_symmetric_energy=symmetric_energy(np.array([D@M@D]), np.array([0.])),
    Ainf_antisymmetric_max_abs=float(np.max(abs(Ainf-Ainf.T))),
    feedback=dict(selected_poses=selected,
        maximum_spectral_radius=float(radius.max()), time_at_max_s=float(time[int(np.argmax(radius))]),
        maximum_scaled_singular_norm=float(norm.max()),
        real_eigenvalue_min=float(real_min.min()), real_eigenvalue_max=float(real_max.max()),
        delay_s=1e-7, dt_s=.02, extrapolation_r=1-1e-7/.02,
        two_tap_formula='aDelayed[n]=(1+r)*a[n-1]-r*a[n-2]; aFeedback=-S_world^-1*Ar*aDelayed',
        caveat='Instantaneous residual gain only; neither spectral_radius(C)<1 nor >1 proves closed-loop stability; no characteristic roots/poles calculated'),
    passivity_meaning='Negative eigenvalue of symmetric real radiation impedance gives a harmonic velocity direction with negative mean radiation power at that sampled frequency; not proof that it causes observed runaway',
    limitations=['FiniteK test is only at actual384incident frequencies, not all frequencies',
                 'Negative instantaneous time-kernel eigenvalues alone would not establish nonpassivity',
                 'Coarsemesh/nativefrequency/Ainf inconsistency and nonlinear rigidbody effects remain unresolved'],
    source_sha256=sha(__file__), input_sha256={str(p):sha(p) for p in
        [npath, rpath, comparison_path, case/'motion.mat', case/'preflight.json']})
out.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
print(json.dumps(report))
