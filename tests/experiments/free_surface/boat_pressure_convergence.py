"""Compare spatial pressure evidence only after verifying matched physical state."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from tests.reference_boat import GRAVITY, LENGTH
from .p4_fluid_ordering import fresh


def compare(paths, output):
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    reports, sources, states, provenance = [], [], [], []
    for name in paths:
        path = Path(name)
        report = json.loads(path.read_text())
        if not report.get('complete'):
            raise ValueError('requires terminal complete pressure audits')
        base = path.with_name(path.stem+'-base.json')
        source = json.loads(base.read_text())
        state = dict(np.load(base.with_suffix('.npz'), allow_pickle=False))
        reports.append(report); sources.append(source); states.append(state)
        provenance.append(dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest(),
            state_arrays_sha256=sha256(base.with_suffix('.npz').read_bytes()).hexdigest()))
    for source, state in zip(sources[1:], states[1:]):
        if (source['actual_sea_realization'] != sources[0]['actual_sea_realization']
                or source['body_mass'] != sources[0]['body_mass']
                or source['fixture_sha256'] != sources[0]['fixture_sha256']):
            raise ValueError('different sea or body fixture: not a matched-state comparison')
        for key in ('parent_nodes', 'parents', 'parent_psi', 'pose_position', 'pose_orientation',
                    'source_body_vertices', 'source_body_faces'):
            if state[key].shape != states[0][key].shape or not np.allclose(state[key], states[0][key], rtol=0, atol=1e-13):
                raise ValueError(f'unmatched physical geometry/state: {key}')
    steps = [[d['step'] for d in report['derivatives']] for report in reports]
    if any(step != steps[0] for step in steps[1:]):
        raise ValueError('compare the same derivative steps at each spatial resolution')
    loads = np.array([report['derivatives'][0]['hydrodynamic_load_world'] for report in reports])
    accelerations = np.array([report['derivatives'][0]['acceleration_world'] for report in reports])
    sensitivity = np.array([report['pressure_load_step_difference_world'] for report in reports])
    delta = np.diff(loads, axis=0)
    weight = sources[0]['body_mass']*GRAVITY
    dynamic = loads-np.array([0., 0., weight, 0., 0., 0.])
    axes = []
    for i, name in enumerate(('surge', 'sway', 'heave', 'roll', 'pitch', 'yaw')):
        axes.append(dict(axis=name, load_unit='N' if i < 3 else 'Nm',
            acceleration_unit='m/s^2' if i < 3 else 'rad/s^2',
            loads=loads[:, i].tolist(), dynamic_loads=dynamic[:, i].tolist(),
            accelerations=accelerations[:, i].tolist(),
            spatial_differences=delta[:, i].tolist(),
            absolute_difference_contraction=float(abs(delta[1, i]/delta[0, i])) if delta[0, i] else None,
            last_relative_change=float(abs(delta[-1, i]/loads[-1, i])) if loads[-1, i] else None,
            crosses_zero=bool(np.any(loads[1:, i]*loads[:-1, i] < 0)),
            derivative_step_changes=sensitivity[:, i].tolist(),
            spatial_change_comparable_to_step_sensitivity=(abs(delta[:, i]) <=
                abs(sensitivity[:-1, i])+abs(sensitivity[1:, i])).tolist()))
    result = dict(qualified=False,
        scope='Three matched-state fluid boundary refinements in a finite diagnostic tank; no physical pass inferred from a family floor',
        inputs=provenance, derivative_steps=steps[0], axes=axes,
        force_equivalent_loads=(loads/np.r_[np.ones(3), np.full(3, LENGTH)]).tolist(),
        interpretation='Contraction is an observation, not an error bound; sign crossings and small responses require absolute acceleration/load evidence.',
        levels=[dict(nodes=(r.get('geometry_checkpoints') or [dict(nodes=len(s['points']))])[0]['nodes'],
                     faces=(r.get('geometry_checkpoints') or [dict(faces=len(s['faces']))])[0]['faces'],
                     unresolved_flux=r['surface']['unresolved_kinematic_mass_relative'],
                     audit_seconds=r['seconds'], peak_rss_bytes=r.get('process_peak_rss_bytes'),
                     reserved_gpu_pool_bytes=r.get('gpu_pool_reserved_peak_bytes')) for r, s in zip(reports, states)])
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs=3)
    parser.add_argument('--output', required=True)
    print(json.dumps(compare(**vars(parser.parse_args())), indent=2))
