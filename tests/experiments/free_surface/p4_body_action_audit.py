"""Six finite-pose body directions of the actual nonlinear canoe action.

One fresh nonlinear stage, then 24 geometry-only scalar-action evaluations.
The six directions are world translation/rotation tangents, holding nodal phi
and canonical configuration rates fixed. Their contractions test generalized
action forces; finite-pose rotational values are not raw pressure torques.
"""
import argparse
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from .canonical_canoe import body_map
from .p4_whole_geometry import P4WholeCanoeGeometry
from .p4_canoe_stage import initial_state, fixed_action
from .p4_canoe_field import prepare
from .p4_variational_canoe import encode, evaluate
from .p4_first_variation import _hash, _progress, _revision


def run():
    """Check each force component independently, without perturbed PDE solves."""
    owner = P4WholeCanoeGeometry.make()
    _progress('actual whole-canoe geometry built')
    q, psi, velocity = initial_state(owner)
    size = len(q)
    data = prepare(owner, q)
    if (data['volume_order'], data['boundary_order']) != (8, 9):
        raise ValueError('fixed-action audit requires the verified primal quadrature8/9')
    _progress('one fresh finite-pose current operator prepared')
    canonical = encode(data, psi, velocity)
    stage = evaluate(owner, data, canonical[size:])
    _progress('full nonlinear canonical stage evaluated')
    action = fixed_action(owner, stage['potential'], stage['qdot'])
    directions = np.zeros((size, 6))
    # dr=J_left(r)^-1*e produces a unit WORLD angular tangent at the base pose.
    # All six paths hold canonical qdot fixed, so world twist varies with K(q).
    directions[:6] = np.linalg.solve(body_map(q[3:6]), np.eye(6))
    expected = directions.T@stage['rate'][size:]
    names = ('world x translation', 'world y translation', 'world z translation',
             'world x angular tangent', 'world y angular tangent', 'world z angular tangent')
    arrays = dict(q=q, canonical=canonical, canonical_rate=stage['rate'],
        physical=stage['physical_state'], fixed_potential=stage['potential'],
        fixed_canonical_qdot=stage['qdot'], directions=directions,
        expected_generalized_action_force=expected,
        reference_points=owner.mesh.points, mesh_elements=owner.mesh.elements,
        free_nodes=owner.free_nodes)
    levels = []
    for level, step in enumerate((2e-4, 1e-4)):
        actual, pair_values = [], []
        for axis, name in enumerate(names):
            plus = action(q+step*directions[:, axis])
            minus = action(q-step*directions[:, axis])
            pair_values.append((plus, minus))
            actual.append((plus-minus)/(2*step))
            _progress(f'{name}, step {step:g}, scalar-action pair complete')
        actual = np.asarray(actual)
        error = actual-expected
        arrays[f'level{level}_scalar_action_pairs'] = np.asarray(pair_values)
        arrays[f'level{level}_finite_difference'] = actual
        arrays[f'level{level}_errors'] = error
        levels.append(dict(step=step, components=[dict(direction=name,
            coordinate_step_units='m' if axis < 3 else 'rad',
            generalized_force_units='N' if axis < 3 else 'N m',
            finite_difference=float(actual[axis]), expected=float(expected[axis]),
            signed_error=float(error[axis]), absolute_error=float(abs(error[axis])),
            relative_error=float(abs(error[axis]/expected[axis])) if expected[axis] else None)
            for axis, name in enumerate(names)]))
    ratios = []
    for axis, name in enumerate(names):
        coarse = levels[0]['components'][axis]['absolute_error']
        fine = levels[1]['components'][axis]['absolute_error']
        ratios.append(dict(direction=name, fine_over_coarse=fine/coarse if coarse else None))
    files = ['p4_body_action_audit.py', 'p4_canoe_stage.py', 'p4_canoe_field.py',
        'p4_variational_canoe.py', 'canonical_canoe.py', 'p4_whole_geometry.py',
        'p4_volume_map.py', 'p4_canoe_contact.py', 'p4_geometry.py', 'p4_surface.py',
        'p4_boundary_action.py', 'p4_tetra.py']
    report = dict(qualified=False, case='actual nonlinear P4 canoe six body action directions',
        revision=_revision(), source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                            for name in files},
        package_versions={name: version(name) for name in ('meshpy', 'numpy', 'scipy')},
        master_hashes=dict(reference_points_sha256=_hash(owner.mesh.points, '<f8'),
            elements_sha256=_hash(owner.mesh.elements, '<i8'),
            free_nodes_sha256=_hash(owner.free_nodes, '<i8')),
        body_configuration=q[:6].tolist(), world_velocity=velocity.tolist(), eta_modes=len(psi),
        prepared_current_fluid_operators=1, full_canonical_evaluations=1,
        perturbed_scalar_action_evaluations=24, volume_quadrature_order=8,
        boundary_quadrature_order=9,
        quadrature_scope='Both finite-difference primal action and analytic stage use8/9. '
            'This checks geometry derivatives at that quadrature; it is NOT quadrature refinement '
            'or an independent higher-quadrature physical energy-work audit.',
        stages=dict(weak_normal_relative_residual=stage['full_weak_normal_velocity_relative_residual'],
            fluid_work_identity_absolute_residual_joules=stage['fluid_work_identity_absolute_residual']),
        levels=levels, fine_over_coarse_error_ratios=ratios,
        interpretation='Fixed nodal phi and canonical qdot in L=rigid KE+rho*boundary velocity '
            'action-fluid KE-relative gravity. World angular tangents use inv(body_map(r)); '
            'the reported Lq contractions are generalized action forces, not raw hydrodynamic '
            'pressure torques. Finite-difference rotation paths have the correct tangent at the '
            'base configuration, not an exact finite world-axis rotation at every offset.',
        scope='Six independent body directions at the original finite-pose/nonflat/all-six-velocity '
            'case. No independent fluid solve at perturbed q. Does not establish arbitrary '
            'surface-mode derivatives, longer trajectories, spatial/Capytaine or throughput '
            'qualification. Raw errors/ratios are retained without selecting output-fitted gates.')
    return report, arrays


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-body-action-audit.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        report, arrays = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output.with_suffix('.npz'), **arrays)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
