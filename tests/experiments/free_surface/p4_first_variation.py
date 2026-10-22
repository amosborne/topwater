"""One mixed-direction rest linearization of the actual nonlinear P4 canoe.

One real rest factor supplies the independent linear block prediction; four
fresh CURRENT nonlinear stages supply centered differences. This is not six
independent columns, a radiation transient, or Capytaine qualification.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.sparse.linalg import splu
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from topwater.hydrostatics import hydrostatics
from .p4_whole_geometry import P4WholeCanoeGeometry
from .p4_canoe_field import prepare
from .p4_variational_canoe import encode, evaluate
from .p4_canoe_stage import verify_master_identity


def _hash(array, dtype):
    return sha256(np.asarray(array, dtype=dtype).tobytes()).hexdigest()


def _progress(label):
    print(json.dumps(dict(progress=True, stage=label)), flush=True)


def _revision():
    return subprocess.check_output(['git', 'rev-parse', 'HEAD'],
        cwd=Path(__file__).resolve().parents[3], text=True).strip()


def _summary(actual, expected, groups):
    result = {}
    for name, (rows, units) in groups.items():
        error = actual[rows]-expected[rows]
        scale = float(np.linalg.norm(expected[rows]))
        result[name] = dict(units=units,
            maximum_absolute_error=float(np.max(abs(error))),
            reference_euclidean_norm=scale,
            error_euclidean_norm=float(np.linalg.norm(error)),
            relative_euclidean_error=float(np.linalg.norm(error)/scale) if scale else None)
    return result


def run(stage_path, rest_checks_path, physical_audit_path, *, derivative_step=1e-4):
    """Five fresh factorizations; no frozen operators or dense Schur matrix."""
    stage_path = Path(stage_path)
    rest_checks_path, physical_audit_path = Path(rest_checks_path), Path(physical_audit_path)
    rest_checks = json.loads(rest_checks_path.read_text())
    physical_audit = json.loads(physical_audit_path.read_text())
    hydrostatics_hash = sha256(Path(hydrostatics.__code__.co_filename).read_bytes()).hexdigest()
    if hydrostatics_hash != rest_checks['independent_hydrostatics_source_sha256']:
        raise ValueError('independent hydrostatics differs from the rest geometry check')
    source_hash = sha256(stage_path.read_bytes()).hexdigest()
    if physical_audit['source_npz_sha256'] != source_hash:
        raise ValueError('physical audit belongs to a different source NPZ')
    if any(value != 0 for value in physical_audit['fresh_physical_base_maximum_absolute_errors'].values()):
        raise ValueError('source ordering requires the documented exact fresh physical base match')
    with np.load(stage_path) as saved:
        source_q = np.array(saved['q'])
        physical = np.array(saved['physical'])
        canonical = np.array(saved['canonical'])
    if (source_q.ndim != 1 or len(source_q) < 7
            or canonical.shape != (2*len(source_q),)
            or physical.shape != (13+2*(len(source_q)-6),)
            or not np.array_equal(canonical[:len(source_q)], source_q)
            or not all(np.all(np.isfinite(a)) for a in (source_q, physical, canonical))):
        raise ValueError('incompatible finite saved physical/canonical stage')

    owner = P4WholeCanoeGeometry.make()
    master_identity = verify_master_identity(owner, stage_path)
    _progress('whole geometry built')
    hashes = dict(reference_points_sha256=_hash(owner.mesh.points, '<f8'),
        elements_sha256=_hash(owner.mesh.elements, '<i8'),
        free_nodes_sha256=_hash(owner.free_nodes, '<i8'))
    for name, value in hashes.items():
        if value != rest_checks[name]:
            raise ValueError('regenerated geometry differs from rest check: '+name)
    for name in ('elements_sha256', 'free_nodes_sha256'):
        if hashes[name] != physical_audit['regenerated_'+name]:
            raise ValueError('regenerated ordering differs from physical audit: '+name)
    count = len(owner.free_nodes)
    size = 6+count
    if len(source_q) != size:
        raise ValueError('source eta dimension differs from regenerated master space')
    orientation = (Rotation.from_rotvec(source_q[3:6])*Rotation.from_quat(
        owner.domain.rest.orientation, scalar_first=True)).as_quat(scalar_first=True)
    if (not np.array_equal(physical[:3], source_q[:3])
            or not np.array_equal(physical[13:13+count], source_q[6:])
            or min(np.max(abs(physical[3:7]-orientation)),
                   np.max(abs(physical[3:7]+orientation))) > 1e-13):
        raise ValueError('source physical configuration does not match its canonical chart')
    qrest = np.r_[owner.domain.rest.position, np.zeros(3+count)]
    direction = source_q-qrest
    velocity = physical[7:13].copy()
    psi = physical[13+count:].copy()
    if not all(np.any(abs(value) > 0) for value in (direction[:3], direction[3:6], direction[6:], psi)):
        raise ValueError('source must contain pose, surface height and potential disturbances')
    if np.any(velocity == 0):
        raise ValueError('source must exercise all six world velocity components')

    rest = prepare(owner, qrest)
    _progress('real rest operator prepared')
    mass = rest['mass']
    mass_factor = splu(mass.tocsc())
    free, inside = rest['free'], rest['inside']
    energy = rest['energy']
    interior = rest['factor'].solve(-energy[inside][:, free]@psi)
    schur_psi = energy[free][:, free]@psi+energy[free][:, inside]@interior
    # C is physical body flux after elimination, whereas stored coupling is D=C+Cv.
    coupling = rest['coupling']-rest['body_transport']
    eta_rate = mass_factor.solve(schur_psi-coupling@velocity)
    stiffness = hydrostatics(owner.domain.body, owner.domain.rest).stiffness
    expected_rate = np.r_[velocity, eta_rate, -stiffness@direction[:6],
                          -1000*9.81*(mass@direction[6:])]
    expected_momentum = np.r_[rest['effective_mass']@velocity+1000*coupling.T@psi,
                              1000*(mass@psi)]
    rest_identity = dict(eta_map_minus_mass_maximum=float(abs(rest['eta_map']-mass).max()),
                        body_transport_maximum=float(np.max(abs(rest['body_transport']))))
    # Do not retain the rest fluid LU while preparing any perturbed stage.
    del rest, energy, interior

    rate_groups = dict(position_rate=(slice(0, 3), 'm/s'),
        local_rotation_rate=(slice(3, 6), 'rad/s'),
        surface_height_rate=(slice(6, size), 'm/s'),
        body_linear_momentum_rate=(slice(size, size+3), 'N'),
        body_angular_momentum_rate=(slice(size+3, size+6), 'N m'),
        surface_momentum_rate=(slice(size+6, 2*size), 'N'))
    momentum_groups = dict(body_linear_momentum=(slice(0, 3), 'kg m/s'),
        body_angular_momentum=(slice(3, 6), 'kg m^2/s'),
        surface_momentum=(slice(6, size), 'kg m/s'))
    arrays = dict(rest_configuration=qrest, source_configuration=source_q,
        configuration_direction=direction, world_velocity_direction=velocity,
        surface_potential_direction=psi, hydrostatic_stiffness=stiffness,
        expected_rate=expected_rate, expected_momentum=expected_momentum,
        rest_surface_schur_potential_flux=schur_psi, rest_body_surface_coupling=coupling)
    levels = []
    for level, amplitude in enumerate((.02, .01)):
        samples = []
        diagnostics = []
        for sign in (1, -1):
            label = f'amplitude {sign*amplitude:+g}'
            q = qrest+sign*amplitude*direction
            data = prepare(owner, q)
            _progress(label+' current operator prepared')
            state = encode(data, sign*amplitude*psi, sign*amplitude*velocity)
            stage = evaluate(owner, data, state[size:], derivative_step=derivative_step)
            intended = np.r_[q[:3], data['orientation'], sign*amplitude*velocity,
                              q[6:], sign*amplitude*psi]
            mismatch = float(np.max(abs(stage['physical_state']-intended)))
            if mismatch > 1e-10:
                raise ValueError('fresh canonical decoding differs from intended physical state')
            suffix = 'plus' if sign > 0 else 'minus'
            arrays[f'level{level}_{suffix}_canonical'] = state.copy()
            arrays[f'level{level}_{suffix}_rate'] = stage['rate'].copy()
            arrays[f'level{level}_{suffix}_physical'] = stage['physical_state'].copy()
            samples.append((stage['rate'].copy(), state[size:].copy()))
            diagnostics.append(dict(signed_amplitude=sign*amplitude,
                physical_decode_maximum_error=mismatch,
                weak_normal_relative_residual=stage['full_weak_normal_velocity_relative_residual'],
                fluid_work_identity_absolute_residual_joules=stage['fluid_work_identity_absolute_residual']))
            _progress(label+' nonlinear stage complete')
            del data, stage, state
        rate = (samples[0][0]-samples[1][0])/(2*amplitude)
        momentum = (samples[0][1]-samples[1][1])/(2*amplitude)
        for name, value in dict(centered_rate=rate, centered_momentum=momentum,
                rate_error=rate-expected_rate, momentum_error=momentum-expected_momentum).items():
            arrays[f'level{level}_{name}'] = value
        error_eta = rate[6:size]-eta_rate
        eta_norm = float(np.sqrt(max(0., eta_rate@(mass@eta_rate))))
        eta_error = float(np.sqrt(max(0., error_eta@(mass@error_eta))))
        levels.append(dict(amplitude=amplitude, stages=diagnostics,
            rate=_summary(rate, expected_rate, rate_groups),
            momentum=_summary(momentum, expected_momentum, momentum_groups),
            surface_height_rate_rest_mass_L2_error_m2_per_second=eta_error,
            surface_height_rate_rest_mass_L2_relative_error=eta_error/eta_norm if eta_norm else None))
    ratios = {}
    for category in ('rate', 'momentum'):
        ratios[category] = {}
        for name, coarse in levels[0][category].items():
            denominator = coarse['error_euclidean_norm']
            ratios[category][name] = (levels[1][category][name]['error_euclidean_norm']/denominator
                                     if denominator else None)
    source_files = ['p4_first_variation.py', 'p4_canoe_field.py', 'p4_variational_canoe.py',
                    'p4_whole_geometry.py', 'p4_volume_map.py', 'p4_canoe_contact.py',
                    'p4_geometry.py', 'p4_surface.py', 'p4_boundary_action.py', 'p4_tetra.py']
    report = dict(qualified=False, case='actual P4 canoe one mixed-direction rest first variation',
        revision=_revision(), source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                            for name in source_files},
        source_npz=str(stage_path), source_npz_sha256=source_hash,
        rest_checks=str(rest_checks_path), rest_checks_sha256=sha256(rest_checks_path.read_bytes()).hexdigest(),
        physical_audit=str(physical_audit_path),
        physical_audit_sha256=sha256(physical_audit_path.read_bytes()).hexdigest(),
        independent_hydrostatics_source_sha256=hydrostatics_hash,
        master_hashes=hashes, active_nodes=int(len(np.unique(owner.mesh.elements))), eta_modes=count,
        master_ordering_assumption=master_identity,
        rest_identities=rest_identity, derivative_step=derivative_step, amplitudes=levels,
        fine_over_coarse_error_ratios=ratios,
        surface_norm='Constant flat-rest projected consistent mass L^2; units m^2/s for height-rate norm.',
        scope='One real mixed direction, not six independent matrix columns. Undamped Q=0 canonical stages only. No radiation transient, prescribed-motion reaction, spatial convergence, Capytaine or throughput qualification. Ratios are raw diagnostics; truncation and numerical floors are not assigned an order automatically.')
    return report, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', type=Path, default=Path('artifacts/p4-nonlinear-canoe-stage.npz'))
    parser.add_argument('--rest-checks', type=Path, default=Path('artifacts/p4-rest-geometry-checks.json'))
    parser.add_argument('--physical-audit', type=Path, default=Path('artifacts/p4-stage-physical-work.json'))
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-first-variation.json'))
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        report, arrays = run(args.stage, args.rest_checks, args.physical_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output.with_suffix('.npz'), **arrays)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
