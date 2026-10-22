"""Identical-rest closed-tank linear control for saved P4 prescribed heave.

One fixed rest interior LU applies S*psi on demand. No dense Schur matrix or
interior-by-surface extension is built. DOP853 retains the zero-start natural
response and samples exactly the nonlinear run's accepted times.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from topwater.hydrostatics import hydrostatics
from .canoe_radiation import ramped_motion
from .linear_canoe_transient import comparison_norms
from .p4_whole_geometry import P4WholeCanoeGeometry
from .p4_canoe_field import prepare
from .p4_canoe_stage import verify_master_identity
from .p4_first_variation import _hash, _progress, _revision


def compare(report_path, *, rtol=1e-10, atol=1e-12):
    """Run only a rest linear control; never repeat the nonlinear simulation."""
    report_path = Path(report_path)
    report = json.loads(report_path.read_text())
    if (report['case'] != 'closed-tank nonlinear P4 prescribed heave'
            or report['integrated_fluid_variables'] != 'eta,Pi'
            or report['boundary'] != 'closed finite tank, no sponge or PML'
            or report['mode'] != 'heave' or report['density'] != 1000. or report['gravity'] != 9.81):
        raise ValueError('requires the closed-tank nonlinear P4 prescribed-heave case')
    if not all(np.isfinite(v) and v > 0 for v in (rtol, atol)):
        raise ValueError('positive finite linear integration tolerances required')
    npz_path = report_path.with_suffix('.npz')
    with np.load(npz_path) as saved:
        times, physical, loads, prescribed = [saved[name].copy() for name in
            ('times', 'states', 'hydrodynamic_load', 'prescribed_motion')]
    owner = P4WholeCanoeGeometry.make()
    identity = verify_master_identity(owner, npz_path)
    count = len(owner.free_nodes)
    hashes = dict(reference_points_sha256=_hash(owner.mesh.points, '<f8'),
        elements_sha256=_hash(owner.mesh.elements, '<i8'), free_nodes_sha256=_hash(owner.free_nodes, '<i8'))
    if (hashes != report['master_hashes'] or report['eta_modes'] != count
            or not np.array_equal(report['rest_position'], owner.domain.rest.position)
            or not np.array_equal(report['rest_orientation'], owner.domain.rest.orientation)):
        raise ValueError('linear geometry/master IDs or equilibrium differ from nonlinear case')
    # The comparison must retain the exact field, geometry and motion equations
    # used by the nonlinear run, not merely nearby physical coordinates.
    shared = ['p4_canoe_field.py', 'p4_whole_geometry.py', 'p4_volume_map.py',
        'p4_canoe_contact.py', 'p4_geometry.py', 'p4_surface.py', 'p4_tetra.py', 'canoe_radiation.py',
        'p4_prescribed_canoe.py', 'p4_variational_canoe.py', 'p4_boundary_action.py',
        'canonical_canoe.py', 'canoe_contact.py', 'p4_canoe_release.py']
    for name in shared:
        if sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() != report['source_sha256'][name]:
            raise ValueError('matched-rest source differs from nonlinear case: '+name)
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0)
            or not np.all(np.isfinite(times)) or times[-1] != report['duration']
            or physical.shape != (len(times), 13+2*count) or loads.shape != (len(times), 6)
            or prescribed.shape != (len(times), 3)
            or not all(np.all(np.isfinite(v)) for v in (physical, loads, prescribed))):
        raise ValueError('invalid saved accepted physical fields/times/loads/motion')
    if np.any(physical[0, 7:] != 0):
        raise ValueError('matched transient requires exactly zero initial body velocity and eta/psi')

    def motion(time):
        value = np.zeros((3, 6))
        value[:, 2] = ramped_motion(time, report['amplitude'], report['omega'], report['ramp_periods'])
        return value

    motions = np.asarray([motion(t) for t in times])
    if not np.array_equal(motions[:, :, 2], prescribed):
        raise ValueError('saved heave displacement/velocity/acceleration differs from metadata')
    expected_body = np.tile(np.r_[owner.domain.rest.position, owner.domain.rest.orientation,
                                  np.zeros(6)], (len(times), 1))
    expected_body[:, 2] += prescribed[:, 0]
    expected_body[:, 9] = prescribed[:, 1]
    if not np.allclose(physical[:, :13], expected_body, rtol=1e-13, atol=1e-15):
        raise ValueError('saved physical body states differ from the prescribed heave case')
    q = np.r_[owner.domain.rest.position, np.zeros(3+count)]
    data = prepare(owner, q)
    _progress('one identical real P4 rest operator prepared')
    mass = data['mass']; mass_factor = splu(mass.tocsc())
    if float(abs(data['eta_map']-mass).max()) > 1e-13 or np.max(abs(data['body_transport'])) > 1e-13:
        raise ValueError('rest Ceta=M and Cv=0 identities failed')
    free, inside = data['free'], data['inside']
    k_if = data['energy'][inside][:, free]
    k_fi = data['energy'][free][:, inside]
    k_ff = data['energy'][free][:, free]
    coupling = data['coupling']-data['body_transport']
    hydro = hydrostatics(owner.domain.body, owner.domain.rest)
    buoyancy = np.r_[hydro.force, hydro.moment]
    calls = 0

    def rate(time, state):
        nonlocal calls
        calls += 1
        eta, psi = state[:count], state[count:]
        interior = data['factor'].solve(-k_if@psi)
        schur_psi = k_ff@psi+k_fi@interior
        return np.r_[mass_factor.solve(schur_psi-coupling@motion(time)[1]), -9.81*eta]

    solved = solve_ivp(rate, (0., times[-1]), np.zeros(2*count), t_eval=times,
                       method='DOP853', rtol=rtol, atol=atol)
    if not solved.success:
        raise RuntimeError(solved.message)
    linear_states = solved.y.T
    psi_rate = -9.81*linear_states[:, :count]
    linear_loads = (buoyancy-motions[:, 0]@hydro.stiffness.T
        -motions[:, 2]@data['added_mass'].T-1000*psi_rate@coupling)
    linear = dict(times=times, eta=linear_states[:, :count], psi=linear_states[:, count:],
        hydrodynamic_load=linear_loads, prescribed_motion=motions)
    summary = dict(qualified=False, case='matched-rest closed-tank P4 heave transient comparison',
        source=str(report_path), source_report_sha256=sha256(report_path.read_bytes()).hexdigest(),
        source_npz_sha256=sha256(npz_path.read_bytes()).hexdigest(), revision=_revision(),
        comparison_script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        comparison_norms_source_sha256=sha256(Path(comparison_norms.__code__.co_filename).read_bytes()).hexdigest(),
        independent_hydrostatics_source_sha256=sha256(Path(hydrostatics.__code__.co_filename).read_bytes()).hexdigest(),
        master_hashes=hashes, master_identity=identity, integrator='DOP853',
        linear_rtol=rtol, linear_atol=atol, linear_rhs_evaluations=calls,
        rest_fluid_factorizations=1, dense_surface_schur_or_extension_built=False,
        expected_equations='eta_dot=M^-1(S*psi-C*V), psi_dot=-g*eta; '
            'hydro=buoyancy-K*xi-A*acceleration-rho*C.T*psi_dot; no mean removal or damping.',
        scope='Identical finite closed-tank linear zero-start control, including startup and all '
            'six load rows at exact nonlinear accepted times. No amplitude/time/spatial '
            'qualification from a single comparison; no outgoing/Capytaine or throughput claim.',
        **comparison_norms(times, physical, loads, linear, mass, buoyancy))
    return summary, linear


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--atol', type=float, default=1e-12)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    with threadpool_limits(1):
        summary, arrays = compare(args.report, rtol=args.rtol, atol=args.atol)
    output = args.output or args.report.with_name(args.report.stem+'-linear-control.json')
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output.with_suffix('.npz'), **arrays)
    output.write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)
