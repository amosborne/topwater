"""One actual canoe predictor: lagged setup, exact current coupled operator.

This is an operator/solve experiment, not a trajectory. The predictor is the
saved full canonical state plus delta times its actual saved nonlinear rate.
The free-surface potential becomes nonzero. Two old preconditioners are checked
against one fresh current direct elimination, using the same current blocks.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from .p4_canoe_field import assemble
from .p4_canoe_stage import verify_master_identity
from .p4_coupled_solve import solve, body_preconditioner
from .p4_fluid_ordering import load_input, fresh
from .p4_whole_geometry import P4WholeCanoeGeometry


def run(saved_matrix, stage, output, *, delta=1e-4):
    output = fresh(output)
    if not np.isfinite(delta) or delta <= 0:
        raise ValueError('positive finite predictor interval required')
    old_matrix, old, manifest = load_input(saved_matrix)
    stage = Path(stage)
    if sha256(stage.read_bytes()).hexdigest() != manifest['provenance']['stage_sha256']:
        raise ValueError('predictor seed differs from the captured initial matrix state')
    with np.load(stage, allow_pickle=False) as saved:
        canonical, rate = saved['canonical'].copy(), saved['canonical_rate'].copy()
    size = len(canonical)//2
    if (canonical.shape != rate.shape or canonical.shape != (2*size,)
            or not np.all(np.isfinite(canonical)) or not np.all(np.isfinite(rate))
            or not np.array_equal(canonical[:size], old['q'])):
        raise ValueError('finite matching full canonical seed and nonlinear rate required')
    report = dict(qualified=False, complete=False,
        scope='One Euler-predicted operator state, not an accepted integration step',
        delta=delta, input_npz_sha256=manifest['npz_sha256'],
        stage_npz_sha256=sha256(stage.read_bytes()).hexdigest(),
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('p4_operator_reuse.py', 'p4_coupled_solve.py', 'p4_canoe_field.py', 'p4_geometry.py')},
        phase_wall_seconds={}, cases=[])
    fields = {}
    failure = None

    def phase(name, function):
        start = perf_counter()
        try:
            return function()
        finally:
            report['phase_wall_seconds'][name] = perf_counter()-start
            print(json.dumps(dict(phase=name, seconds=report['phase_wall_seconds'][name])), flush=True)

    def lagged_solves(matrix, load, mass, momentum, forcing):
        # The old fluid factor exists only in this frame. Release it before
        # constructing the independent current direct reference below.
        factor = splu(old_matrix, permc_spec='MMD_AT_PLUS_A')
        inverse = body_preconditioner(factor.solve, old['response'], old['rigid_mass']+old['added_mass'])
        results = {}
        for name, preconditioner in (('old_fluid', factor.solve), ('old_coupled', inverse)):
            start = perf_counter()
            current = solve(matrix, load, mass, momentum, forcing, preconditioner,
                            initial=old['reference_potential'], rtol=1e-11)
            report['cases'].append(dict(name=name, wall_seconds=perf_counter()-start,
                iterations=current['iterations'], current_relative_residual=current['coupled_relative_residual']))
            results[name] = current
        return results

    try:
        owner = phase('geometry_make', P4WholeCanoeGeometry.make)
        report['master_identity'] = verify_master_identity(owner, stage)
        predicted = canonical+delta*rate
        data = phase('current_assembly', lambda: assemble(owner, predicted[:size]))
        inside, free = data['inside'], data['free']
        if not np.array_equal(inside, old['inside']) or not np.array_equal(free, old['free']):
            raise ValueError('master IDs changed; old preconditioner reuse is invalid')
        matrix = data['energy'][inside][:, inside].tocsc()
        load, mass = data['load'][inside], data['rigid_mass']
        psi = data['eta_factor'].solve(predicted[size+6:]/1000., trans='T')
        pworld = np.linalg.solve(data['body_map'].T, predicted[size:size+6])
        forcing = -data['energy'][inside][:, free]@psi
        momentum = pworld-1000*(data['load'][free]+data['body_transport']).T@psi
        report.update(surface_potential_l2=float(np.linalg.norm(psi)),
            current_matrix_relative_change=float(np.linalg.norm((matrix-old_matrix).data)/np.linalg.norm(old_matrix.data)))
        results = phase('old_factor_and_current_solves', lambda: lagged_solves(matrix, load, mass, momentum, forcing))
        factor = phase('fresh_current_factor', lambda: splu(matrix, permc_spec='MMD_AT_PLUS_A'))
        response = phase('fresh_six_body_responses', lambda: factor.solve(load))
        added = 1000*load.T@response
        coupling = data['load'][free]-data['energy'][free][:, inside]@response+data['body_transport']
        reference_v = np.linalg.solve(mass+added, pworld-1000*coupling.T@psi)
        reference_phi = response@reference_v+factor.solve(forcing)
        fields.update(predicted_canonical=predicted, psi=psi,
                      reference_velocity=reference_v, reference_potential=reference_phi)
        for case in report['cases']:
            current = results[case['name']]
            case.update(velocity_maximum_absolute_difference=float(np.max(abs(current['velocity']-reference_v))),
                potential_relative_difference=float(np.linalg.norm(current['potential']-reference_phi)/np.linalg.norm(reference_phi)),
                fluid_residual_l2=float(np.linalg.norm(current['fluid_residual'])),
                body_residual_l2=float(np.linalg.norm(current['body_residual'])))
            fields[case['name']+'_potential'] = current['potential']
            fields[case['name']+'_velocity'] = current['velocity']
        report['complete'] = True
    except BaseException as error:
        failure = error
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        try:
            output.parent.mkdir(parents=True, exist_ok=True)
            np.savez(output.with_suffix('.npz'), **fields)
            output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        except BaseException as save_error:
            if failure is not None:
                failure.add_note('Operator-reuse evidence could not be saved: '+repr(save_error))
            else:
                raise
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--stage', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--delta', type=float, default=1e-4)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        print(json.dumps(run(args.input, args.stage, args.output, delta=args.delta), indent=2))
