"""Offline SAMEcase loose/tight nonlinear P4 release diagnostics.

All endpoint differences use the exact shared final time. Unequal accepted
grids are not interpolated and misreported as temporal error. No trajectory,
field operator or fluid solve is rerun; raw histories remain untouched.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation


def read(path):
    report = json.loads(path.read_text())
    with np.load(path.with_suffix('.npz')) as saved:
        arrays = {name: saved[name].copy() for name in ('times', 'states', 'canonical_states',
            'hamiltonian', 'fluid_volume', 'rim_clearance_lower_bound',
            'weak_normal_relative_residual', 'fluid_work_identity_absolute_residual')}
    times = arrays['times']; count = report['eta_modes']; states = arrays['states']
    if (report['case'] != 'fully nonlinear P4 canonical canoe release' or len(times) < 2
            or times.ndim != 1 or times[0] != 0 or np.any(np.diff(times) <= 0)
            or times[-1] != report['duration'] or report['simulated_seconds'] != times[-1]
            or states.shape != (len(times), 13+2*count)
            or arrays['canonical_states'].shape != (len(times), 2*(6+count))
            or any(array.shape != (len(times),) for name, array in arrays.items()
                   if name not in ('states', 'canonical_states'))
            or not all(np.all(np.isfinite(value)) for value in arrays.values())):
        raise ValueError('complete finite accepted P4 release histories required')
    quaternion_error = np.max(abs(np.linalg.norm(states[:, 3:7], axis=1)-1))
    if quaternion_error > 1e-10:
        raise ValueError('saved physical quaternions violate the unit-quaternion convention')
    if not all(np.isfinite(report[key]) and report[key] > 0 for key in
               ('rtol', 'atol', 'max_step', 'derivative_step')):
        raise ValueError('positive finite integration/shape controls required')
    return report, arrays


def _summary(arrays):
    states = arrays['states']; h = arrays['hamiltonian']; volume = arrays['fluid_volume']
    orientation = Rotation.from_quat(states[:, 3:7], scalar_first=True)
    world_rotation = (orientation*orientation[0].inv()).as_rotvec()
    velocity = states[:, 7:13]
    return dict(accepted_steps=len(states)-1, initial_hamiltonian_joules=float(h[0]),
        final_hamiltonian_change_joules=float(h[-1]-h[0]),
        maximum_accepted_hamiltonian_change_joules=float(np.max(abs(h-h[0]))),
        maximum_accepted_hamiltonian_change_over_initial=(float(np.max(abs(h-h[0]))/abs(h[0])) if h[0] else None),
        final_fluid_volume_change_m3=float(volume[-1]-volume[0]),
        maximum_accepted_fluid_volume_change_m3=float(np.max(abs(volume-volume[0]))),
        maximum_position_change_m=np.max(abs(states[:, :3]-states[0, :3]), axis=0).tolist(),
        maximum_world_relative_rotation_components_rad=np.max(abs(world_rotation), axis=0).tolist(),
        maximum_rotation_geodesic_from_initial_rad=float(np.max(np.linalg.norm(world_rotation, axis=1))),
        initial_world_velocity=velocity[0].tolist(), final_world_velocity=velocity[-1].tolist(),
        minimum_world_velocity=velocity.min(axis=0).tolist(), maximum_world_velocity=velocity.max(axis=0).tolist(),
        maximum_world_velocity_change=np.max(abs(velocity-velocity[0]), axis=0).tolist(),
        world_velocity_units=['m/s']*3+['rad/s']*3,
        maximum_quaternion_norm_error=float(np.max(abs(np.linalg.norm(states[:, 3:7], axis=1)-1))),
        minimum_rim_clearance_lower_bound_m=float(np.min(arrays['rim_clearance_lower_bound'])),
        maximum_weak_normal_relative_residual=float(np.max(arrays['weak_normal_relative_residual'])),
        maximum_fluid_work_identity_absolute_residual_joules=float(np.max(arrays['fluid_work_identity_absolute_residual'])))


def compare(loose_report, loose, tight_report, tight):
    if (loose['times'][-1] != tight['times'][-1]
            or loose['times'][-1] != loose_report['duration']):
        raise ValueError('exact shared final accepted time required')
    identities = ('source_npz_sha256', 'reference_points_sha256', 'elements_sha256',
                  'free_nodes_sha256', 'eta_modes', 'duration', 'derivative_step', 'runner_sha256')
    runtime = ('source_sha256' in loose_report, 'source_sha256' in tight_report)
    if any(runtime) and (not all(runtime) or loose_report['source_sha256'] != tight_report['source_sha256']):
        raise ValueError('runtime numerical source hashes differ or are missing on one side')
    if any(loose_report[key] != tight_report[key] for key in identities):
        raise ValueError('release initial stage, masters, duration, shape step or runner differ')
    if (loose_report['source_stage_metadata'] != tight_report['source_stage_metadata']
            or not np.array_equal(loose['canonical_states'][0], tight['canonical_states'][0])):
        raise ValueError('release stage metadata or exact initial canonical state differs')
    controls = ('rtol', 'atol', 'max_step')
    if (any(tight_report[key] > loose_report[key] for key in controls)
            or not any(tight_report[key] < loose_report[key] for key in controls)):
        raise ValueError('tight run must tighten at least one integration control and loosen none')
    a, b = loose['states'], tight['states']
    count = loose_report['eta_modes']
    attitude = (Rotation.from_quat(a[-1, 3:7], scalar_first=True)
                *Rotation.from_quat(b[-1, 3:7], scalar_first=True).inv()).as_rotvec()
    initial_attitude = (Rotation.from_quat(a[0, 3:7], scalar_first=True)
                       *Rotation.from_quat(b[0, 3:7], scalar_first=True).inv()).magnitude()
    body_difference = a[-1]-b[-1]
    endpoint = dict(time=float(loose['times'][-1]),
        position_difference_m=body_difference[:3].tolist(),
        world_rotation_difference_vector_rad=attitude.tolist(),
        orientation_geodesic_difference_rad=float(np.linalg.norm(attitude)),
        world_velocity_difference=body_difference[7:13].tolist(),
        initial_nonquaternion_component_maximum_absolute_difference=float(np.max(abs(np.r_[
            a[0, :3]-b[0, :3], a[0, 7:]-b[0, 7:]]))),
        initial_orientation_geodesic_difference_rad=float(initial_attitude),
        hamiltonian_difference_joules=float(loose['hamiltonian'][-1]-tight['hamiltonian'][-1]),
        fluid_volume_difference_m3=float(loose['fluid_volume'][-1]-tight['fluid_volume'][-1]))
    for name, difference in (('eta', body_difference[13:13+count]), ('psi', body_difference[13+count:])):
        endpoint[name+'_maximum_absolute_difference'] = float(np.max(abs(difference)))
        endpoint[name+'_nodal_rms_difference'] = float(np.sqrt(np.mean(difference**2)))
    return dict(qualified=False, case='same-case nonlinear P4 release temporal control',
        identity_keys_checked=list(identities)+(['source_sha256'] if all(runtime) else []),
        runtime_numerical_source_hashes_checked=all(runtime),
        controls=dict(loose={key: loose_report[key] for key in controls},
            tight={key: tight_report[key] for key in controls}),
        loose=_summary(loose), tight=_summary(tight), exact_endpoint_comparison=endpoint,
        shape_derivative_step_unchanged=loose_report['derivative_step'],
        scope='Raw accepted-grid drift and exact common-endpoint differences only. All six pose/velocity '
            'components are retained without a fitted nonzero threshold. Quaternion differences are '
            'sign-invariant. Surface endpoint RMS is nodal, not a physical-area norm. No unequal-grid '
            'interpolation, observed order, temporal pass gate, spatial/outgoing or Capytaine qualification.',
        provenance_limit=('Both reports serialize identical runtime numerical source hashes.' if all(runtime) else
            'Legacy release reports pin the initial stage, masters, runner and revision; the '
            'runtime hashes of all loaded numerical modules are not separately serialized. Equal '
            'stage metadata and runner do not prove an absent runtime-source fingerprint.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('loose', type=Path)
    parser.add_argument('tight', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('choose a new comparison report; trajectories are preserved')
    lr, la = read(args.loose); tr, ta = read(args.tight)
    report = compare(lr, la, tr, ta)
    paths = dict(loose_report=args.loose, tight_report=args.tight,
                 loose_npz=args.loose.with_suffix('.npz'), tight_npz=args.tight.with_suffix('.npz'),
                 comparator=Path(__file__))
    report['input_provenance'] = {name: dict(path=str(path.resolve()), sha256=sha256(path.read_bytes()).hexdigest())
                                  for name, path in paths.items()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
