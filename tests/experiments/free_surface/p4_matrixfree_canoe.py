"""Current-geometry matrix-free field and unchanged nonlinear canoe action.

One-state candidate, not a qualified trajectory. Old setup preconditions the
current equations; the metric cache is rebuilt at the supplied geometry.
No six current auxiliary fields or global volume matrix are constructed.
"""
import argparse
import cProfile
from hashlib import sha256
import json
import platform
from pathlib import Path
from time import perf_counter

import numpy as np
import scipy
from scipy.sparse.linalg import LinearOperator, splu
from threadpoolctl import threadpool_limits, threadpool_info

from .p4_canoe_field import geometry_data
from .p4_canoe_release import _rim_bound
from .p4_canoe_stage import verify_master_identity
from .p4_coupled_solve import solve, body_preconditioner
from .p4_fluid_ordering import load_input, fresh
from .p4_metric_action import MetricAction
from .p4_tetra import face_nodes
from .p4_variational_canoe import action_rate
from .p4_whole_geometry import P4WholeCanoeGeometry


def decode(data, momentum, cache, preconditioner, *, initial=None):
    """Current physical field/twist from canonical momenta, using exact actions."""
    momentum = np.asarray(momentum)
    if momentum.shape != data['q'].shape or not np.all(np.isfinite(momentum)):
        raise ValueError('finite body and surface canonical momenta required')
    inside, free = data['inside'], data['free']
    count = cache.node_count
    action_times = []

    def action(values):
        result, timing = cache.action(values, chunk_size=256)
        action_times.append(timing['total_seconds'])
        return result

    def interior_action(values):
        full = np.zeros(count)
        full[inside] = values
        return action(full)[inside]

    operator = LinearOperator((len(inside),)*2, matvec=interior_action, dtype=float)
    full_operator = LinearOperator((count,)*2, matvec=action, dtype=float)
    psi = data['eta_factor'].solve(momentum[6:]/1000., trans='T')
    trace = np.zeros(count)
    trace[free] = psi
    boundary_action = action(trace) if np.any(psi != 0) else np.zeros(count)
    effective_momentum = (np.linalg.solve(data['body_map'].T, momentum[:6])
        -1000*(data['load'][free]+data['body_transport']).T@psi)
    field = solve(operator, data['load'][inside], data['rigid_mass'],
                  effective_momentum, -boundary_action[inside], preconditioner,
                  initial=initial)
    velocity = field['velocity']
    potential = trace
    potential[inside] = field['potential']
    full_action = action(potential)
    flux = (full_action-data['load']@velocity)[free]
    eta_rate = data['eta_factor'].solve(flux-data['body_transport']@velocity)
    qdot = np.r_[np.linalg.solve(data['body_map'], velocity), eta_rate]
    physical = np.r_[data['q'][:3], data['orientation'], velocity, data['q'][6:], psi]
    hamiltonian = float(500*potential@full_action
                        +.5*velocity@data['rigid_mass']@velocity+data['gravity'])
    data['energy'] = full_operator
    return dict(potential=potential, eta_rate=eta_rate, flux=flux, qdot=qdot,
                physical_state=physical, hamiltonian=hamiltonian,
                operator_action_seconds=action_times,
                solver_iterations=field['iterations'],
                solver_relative_residual=field['coupled_relative_residual'])


def run(stage, old_matrix, output, *, warm_repetitions=0):
    """One actual saved nonlinear state; independently saved full RHS is oracle."""
    output = fresh(output)
    if output.with_suffix('.prof').exists():
        raise FileExistsError('preserve profiling evidence; choose a fresh output stem')
    if not isinstance(warm_repetitions, int) or not 0 <= warm_repetitions <= 3:
        raise ValueError('zero to three warm same-state repetitions required')
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(stage)
    with np.load(stage, allow_pickle=False) as saved:
        canonical = saved['canonical'].copy()
        physical = saved['physical'].copy()
        rate = saved['canonical_rate'].copy()
        reference_phi = saved['fixed_potential'].copy()
    matrix, old, manifest = load_input(old_matrix)
    report = dict(qualified=False, complete=False,
        scope='Full current nonlinear matrix-free RHS at one fixed actual state; no trajectory throughput qualification',
        warm_repetitions=warm_repetitions,
        stage_sha256=sha256(stage.read_bytes()).hexdigest(),
        old_matrix_npz_sha256=manifest['npz_sha256'], phase_wall_seconds={},
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                     native_threadpools=threadpool_info()),
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('p4_matrixfree_canoe.py', 'p4_metric_action.py', 'p4_coupled_solve.py',
                         'p4_canoe_field.py', 'p4_variational_canoe.py')})
    result = None
    failure = None
    phase_times = report['phase_wall_seconds']

    def phase(name, function):
        print(json.dumps(dict(phase=name, event='start')), flush=True)
        start = perf_counter()
        try:
            return function()
        finally:
            phase_times[name] = perf_counter()-start
            print(json.dumps(dict(phase=name, seconds=phase_times[name])), flush=True)

    try:
        owner = phase('geometry_startup', P4WholeCanoeGeometry.make)
        report['master_identity'] = verify_master_identity(owner, stage)
        size = 6+len(owner.free_nodes)
        if (canonical.shape != (2*size,) or rate.shape != canonical.shape
                or physical.shape != (13+2*len(owner.free_nodes),)
                or reference_phi.shape != (len(owner.mesh.points),)
                or not all(np.all(np.isfinite(x)) for x in (canonical, rate, physical, reference_phi))):
            raise ValueError('matching finite actual nonlinear seed/reference arrays required')
        factor = phase('old_fluid_setup', lambda: splu(matrix, permc_spec='MMD_AT_PLUS_A'))
        inverse = body_preconditioner(factor.solve, old['response'], old['rigid_mass']+old['added_mass'])
        report['samples'] = []
        for repetition in range(warm_repetitions+1):
            phase_times = {}
            rhs_started = perf_counter()
            profiler = cProfile.Profile() if repetition == warm_repetitions and repetition > 0 else None
            print(json.dumps(dict(repetition=repetition, profiled=profiler is not None)), flush=True)
            if profiler is not None:
                profiler.enable()
            data = phase('current_boundary_geometry', lambda: geometry_data(owner, canonical[:size]))
            if not np.array_equal(data['inside'], old['inside']) or not np.array_equal(data['free'], old['free']):
                raise ValueError('old preconditioner master ordering differs from the current geometry')
            surface_nodes = face_nodes(owner.mesh, owner.free_faces)
            outer = (owner.mesh.points[:, :2].min(axis=0), owner.mesh.points[:, :2].max(axis=0))
            clearance = _rim_bound(owner, data['state'].points,
                data['q'][:3], data['orientation'], surface_nodes, outer)[0]
            cache = phase('current_metric_cache', lambda: MetricAction(
                data['state'].points[owner.mesh.elements], owner.mesh.elements, len(owner.mesh.points)))
            current = phase('current_coupled_field', lambda: decode(data, canonical[size:], cache,
                inverse, initial=old['reference_potential']))
            kinetic = phase('current_kinetic_energy_and_shape_force', lambda: cache.kinetic(current['potential']))
            result = phase('full_nonlinear_action_rate', lambda: action_rate(owner, data, current,
                kinetic=(kinetic['kinetic_energy'], kinetic['kinetic_gradient'])))
            rhs_wall_seconds = perf_counter()-rhs_started
            sample_report = dict(complete=True, rim_clearance_lower_bound=clearance, solver_iterations=current['solver_iterations'],
                solver_relative_residual=current['solver_relative_residual'],
                metric_persistent_bytes=cache.persistent_array_bytes,
                operator_actions=len(current['operator_action_seconds']),
                operator_action_total_seconds=float(sum(current['operator_action_seconds'])),
                field_relative_difference=float(np.linalg.norm(result['potential']-reference_phi)/np.linalg.norm(reference_phi)),
                physical_state_maximum_absolute_difference=float(np.max(abs(result['physical_state']-physical))),
                canonical_rate_relative_difference=float(np.linalg.norm(result['rate']-rate)/np.linalg.norm(rate)),
                canonical_rate_maximum_absolute_difference=float(np.max(abs(result['rate']-rate))),
                rate_block_maximum_absolute_difference={name: float(np.max(abs(result['rate'][block]-rate[block])))
                    for name, block in [('body_qdot', slice(0, 6)), ('eta_dot', slice(6, size)),
                                        ('body_pdot', slice(size, size+6)), ('surface_pdot', slice(size+6, None))]},
                weak_normal_relative=result['full_weak_normal_velocity_relative_residual'],
                fluid_work_residual=result['fluid_work_identity_absolute_residual'])
            if profiler is not None:
                profiler.disable()
                profiler.dump_stats(output.with_suffix('.prof'))
            sample_report.update(repetition=repetition, profiled=profiler is not None,
                                 phase_wall_seconds=phase_times,
                                 rhs_wall_seconds=rhs_wall_seconds)
            report['samples'].append(sample_report)
        report.update(complete=True)

    except BaseException as error:
        failure = error
        report['incomplete_phase_wall_seconds'] = phase_times
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        try:
            if result is not None:
                np.savez(output.with_suffix('.npz'), potential=result['potential'],
                    physical=result['physical_state'], canonical_rate=result['rate'],
                    qdot=result['qdot'], hamiltonian=result['hamiltonian'])
            output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        except BaseException as save_error:
            if failure is not None:
                failure.add_note('Matrix-free evidence could not be saved: '+repr(save_error))
            else:
                raise
    return report


if __name__ == '__main__':
    # Reuse the existing Windows component-wheel CUDA DLL setup.
    from . import cuda_operator
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', required=True)
    parser.add_argument('--old-matrix', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--warm-repetitions', type=int, default=0)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        print(json.dumps(run(args.stage, args.old_matrix, args.output, warm_repetitions=args.warm_repetitions), indent=2))
