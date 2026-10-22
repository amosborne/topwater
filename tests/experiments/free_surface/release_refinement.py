"""Endpoint refinement diagnostic for completed nonlinear six-DOF releases.

Two trajectories bound a measured difference; they do not establish temporal
order, spatial accuracy, long-time stability or Capytaine qualification.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from .moving_canoe import make_domain
from .surface_damping import matrices


def compare(coarse_path, fine_path):
    paths = [Path(coarse_path), Path(fine_path)]
    reports = [json.loads(path.read_text()) for path in paths]
    coarse, fine = reports
    if any(not np.isfinite(report['rtol']) or not 0 < report['rtol'] < 1
           or not np.isfinite(report['max_step']) or report['max_step'] <= 0
           for report in reports):
        raise ValueError('finite positive timestep and 0 < tolerance < 1 required')
    for name in ('case', 'n', 'wet_subdivisions', 'initial_scale', 'simulated_seconds'):
        if coarse[name] != fine[name]:
            raise ValueError(f'release metadata differs: {name}')
    if (coarse['case'] != 'fully nonlinear variational six-DOF canoe'
            or not fine['max_step'] < coarse['max_step'] or fine['rtol'] > coarse['rtol']):
        raise ValueError('same variational release with smaller step and no looser tolerance required')
    domain = make_domain(coarse['n'], False, coarse['wet_subdivisions'])
    size = len(domain.free_nodes)
    states = []
    for path, report in zip(paths, reports):
        with np.load(path.with_suffix('.npz')) as saved:
            times, values = saved['times'], saved['states']
            if (times.ndim != 1 or len(times) < 2 or not np.all(np.isfinite(times))
                    or values.shape != (len(times), 13+2*size) or not np.all(np.isfinite(values))
                    or times[0] != 0 or times[-1] != report['simulated_seconds']
                    or np.any(np.diff(times) <= 0)):
                raise ValueError('incomplete or incompatible saved release')
            states.append(values[[0, -1]].copy())
    if not np.allclose(states[0][0], states[1][0], rtol=0, atol=1e-14):
        raise ValueError('release initial physical states differ')
    a, b = states[0][-1], states[1][-1]
    if not np.allclose([np.linalg.norm(a[3:7]), np.linalg.norm(b[3:7])], 1., atol=1e-10, rtol=0):
        raise ValueError('saved orientations are not unit quaternions')
    rotation = (Rotation.from_quat(b[3:7], scalar_first=True)
                * Rotation.from_quat(a[3:7], scalar_first=True).inv()).as_rotvec()
    mass, _ = matrices(domain.reference_points, domain.free_faces,
                       domain.free_nodes, np.zeros(size))
    result = dict(qualified=False, sources=[str(path) for path in paths],
        simulated_seconds=fine['simulated_seconds'],
        endpoint_position_difference_m=(b[:3]-a[:3]).tolist(),
        endpoint_world_rotation_difference_radians=rotation.tolist(),
        endpoint_world_velocity_difference=(b[7:13]-a[7:13]).tolist(),
        endpoint_world_velocity_units=['m/s']*3+['rad/s']*3,
        scope='Fine minus coarse at the identical final time; no trajectory interpolation. '
              'Surface L2 uses the flat rest projected consistent mass and shared configuration '
              'node labels. Matching implementation/geometry must be established separately; '
              'these legacy reports have no source fingerprint. Fine is not an exact solution. Two levels do '
              'not establish convergence order or nonlinear/spatial/Capytaine qualification.')
    for name, offset in (('eta', 13), ('psi', 13+size)):
        difference, reference = b[offset:offset+size]-a[offset:offset+size], b[offset:offset+size]
        norm = float(np.sqrt(max(0., difference@(mass@difference))))
        reference_norm = float(np.sqrt(max(0., reference@(mass@reference))))
        result[name+'_endpoint_reference_domain_l2_difference'] = norm
        result[name+'_endpoint_reference_domain_l2_fine'] = reference_norm
        result[name+'_endpoint_relative_l2_difference'] = norm/reference_norm if reference_norm else None
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coarse', type=Path)
    parser.add_argument('fine', type=Path)
    parser.add_argument('--output', type=Path, default=Path('artifacts/nonlinear-release-refinement.json'))
    args = parser.parse_args()
    report = compare(args.coarse, args.fine)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
