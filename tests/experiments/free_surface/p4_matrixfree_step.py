"""One bounded adaptive nonlinear canoe step and its actual error components.

This uses the complete exact-current matrix-free stage, without a surrogate
operator or changed tolerances. It is not a trajectory or performance gate.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
from time import perf_counter

import numpy as np
import scipy
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits, threadpool_info

from tests.benchmark_rankine import peak_rss_bytes
from .geometry_rk23 import GeometryRK23
from .p4_canoe_stage import verify_master_identity
from .p4_coupled_solve import body_preconditioner
from .p4_fluid_ordering import fresh, load_input
from .p4_matrixfree_canoe import evaluate_current
from .p4_whole_geometry import P4WholeCanoeGeometry


def run(stage, old_matrix, output, *, step=1e-4, rtol=1e-4, atol=1e-10):
    """Advance one accepted step, at most 13 RHS calls and a soft 120 s budget."""
    output = fresh(output)
    if not all(np.isfinite(v) and v > 0 for v in (step, rtol, atol)) or rtol >= 1:
        raise ValueError('positive finite controls and relative tolerance below one required')
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(stage)
    with np.load(stage, allow_pickle=False) as saved:
        initial = saved['canonical'].copy()
        saved_rate = saved['canonical_rate'].copy()
        saved_physical = saved['physical'].copy()
    matrix, old, manifest = load_input(old_matrix)
    report = dict(qualified=False, complete=False,
        scope='One accepted full nonlinear six-DOF RK23 step; not trajectory throughput or temporal validation',
        stage_sha256=sha256(stage.read_bytes()).hexdigest(),
        old_matrix_npz_sha256=manifest['npz_sha256'],
        controls=dict(first_step=step, rtol=rtol, scalar_atol=atol, max_rhs=13, rhs_budget_seconds=120.),
        budget_scope='Elapsed budget checked before calls; an in-flight RHS is not preempted.',
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                     native_threadpools=threadpool_info()),
        rhs=[], embedded_errors=[], startup_seconds={},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('p4_matrixfree_step.py', 'p4_matrixfree_canoe.py', 'p4_variational_canoe.py',
             'p4_whole_geometry.py', 'p4_metric_action.py', 'p4_canoe_field.py',
             'p4_coupled_solve.py', 'geometry_rk23.py')})
    arrays = dict(initial_canonical=initial)
    solver = None
    try:
        before = perf_counter()
        owner = P4WholeCanoeGeometry.make()
        report['startup_seconds']['geometry'] = perf_counter()-before
        report['master_identity'] = verify_master_identity(owner, stage)
        size = 6+len(owner.free_nodes)
        if (initial.shape != (2*size,) or saved_rate.shape != initial.shape
                or saved_physical.shape != (13+2*len(owner.free_nodes),)
                or not all(np.all(np.isfinite(v)) for v in (initial, saved_rate, saved_physical))):
            raise ValueError('matching finite actual canoe seed arrays required')
        inside = np.setdiff1d(np.unique(owner.mesh.elements), owner.free_nodes)
        if not np.array_equal(inside, old['inside']) or not np.array_equal(owner.free_nodes, old['free']):
            raise ValueError('cached preconditioner master IDs differ')
        before = perf_counter()
        factor = splu(matrix, permc_spec='MMD_AT_PLUS_A')
        inverse = body_preconditioner(factor.solve, old['response'], old['rigid_mass']+old['added_mass'])
        report['startup_seconds']['old_fluid_factor'] = perf_counter()-before
        latest = {}
        previous_potential = old['reference_potential']
        rhs_started = perf_counter()

        def rhs(time, canonical):
            nonlocal previous_potential
            if len(report['rhs']) >= 13 or perf_counter()-rhs_started > 120:
                raise RuntimeError('bounded adaptive-step diagnostic exhausted its RHS budget')
            row = dict(time=float(time), complete=False)
            report['rhs'].append(row)
            print(json.dumps(dict(rhs=len(report['rhs']), time=time, event='start')), flush=True)
            current = evaluate_current(owner, canonical, inverse, initial=previous_potential)
            previous_potential = current['potential'][inside].copy()
            if len(report['rhs']) == 1:
                report['initial_physical_maximum_absolute_difference'] = float(np.max(abs(current['physical_state']-saved_physical)))
                report['initial_rate_maximum_absolute_difference'] = float(np.max(abs(current['rate']-saved_rate)))
                if report['initial_physical_maximum_absolute_difference'] > 1e-9:
                    raise ValueError('initial decoded physical state differs from saved actual seed')
                arrays['initial_physical'] = current['physical_state'].copy()
                arrays['initial_hamiltonian'] = np.asarray(current['hamiltonian'])
                arrays['initial_fluid_volume'] = np.asarray(current['fluid_volume'])
            row.update(complete=True, wall_seconds=current['rhs_wall_seconds'],
                metric_persistent_bytes=current['metric_persistent_bytes'],
                solver_iterations=current['solver_iterations'], phase_wall_seconds=current['phase_wall_seconds'],
                weak_normal_relative=current['full_weak_normal_velocity_relative_residual'],
                work_residual=current['fluid_work_identity_absolute_residual'])
            latest.update(time=float(time), canonical=canonical.copy(), result=current)
            print(json.dumps(row), flush=True)
            return current['rate']

        blocks = [('position_m', slice(0, 3)), ('rotation_rad', slice(3, 6)),
            ('height_m', slice(6, size)), ('linear_momentum_kg_m_per_s', slice(size, size+3)),
            ('angular_momentum_kg_m2_per_s', slice(size+3, size+6)),
            ('surface_momentum_kg_m_per_s', slice(size+6, None))]
        contact = np.isin(owner.free_nodes, owner.chart_to_global[owner.chart.contact_nodes])

        class AuditedRK23(GeometryRK23):
            def _estimate_error_norm(self, k, h, scale):
                error = self._estimate_error(k, h)
                squared = (error/scale)**2
                total = float(squared.sum())
                value = super()._estimate_error_norm(k, h, scale)
                report['embedded_errors'].append(dict(time=float(self.t), attempted_step=float(h), norm=float(value),
                    contact_height_scaled_squared_norm_fraction=float(squared[6:size][contact].sum()/total) if total else 0.,
                    contact_surface_momentum_scaled_squared_norm_fraction=float(squared[size+6:][contact].sum()/total) if total else 0.,
                    blocks={name: dict(scaled_squared_norm_fraction=float(squared[block].sum()/total) if total else 0.,
                        maximum_absolute_error=float(np.max(abs(error[block]))),
                        minimum_scale=float(np.min(scale[block])), maximum_scale=float(np.max(scale[block])))
                        for name, block in blocks}))
                arrays['last_embedded_error'] = error.copy()
                arrays['last_error_scale'] = scale.copy()
                return value

        solver = AuditedRK23(rhs, 0., initial, step, first_step=step, max_step=step, rtol=rtol, atol=atol)
        message = solver.step()
        if solver.status == 'failed':
            raise RuntimeError(message)
        if latest['time'] != solver.t or not np.array_equal(latest['canonical'], solver.y):
            raise RuntimeError('accepted endpoint has no matching exact-current diagnostics')
        current = latest['result']
        arrays.update(accepted_time=np.asarray(solver.t), accepted_canonical=solver.y.copy(),
            accepted_physical=current['physical_state'], accepted_hamiltonian=np.asarray(current['hamiltonian']),
            accepted_fluid_volume=np.asarray(current['fluid_volume']))
        report.update(complete=True, accepted_steps=1, accepted_time=float(solver.t),
            rhs_wall_seconds=perf_counter()-rhs_started,
            invalid_geometry_trial_rejections=solver.invalid_geometry_rejections,
            hamiltonian_change=float(current['hamiltonian']-arrays['initial_hamiltonian']),
            fluid_volume_change=float(current['fluid_volume']-arrays['initial_fluid_volume']))
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['peak_rss_bytes'] = peak_rss_bytes()
        if solver is not None:
            report['invalid_geometry_trial_rejections'] = solver.invalid_geometry_rejections
        np.savez(output.with_suffix('.npz'), **arrays)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    from . import cuda_operator
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', required=True)
    parser.add_argument('--old-matrix', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--step', type=float, default=1e-4)
    args = parser.parse_args()
    with threadpool_limits(1):
        print(json.dumps(run(args.stage, args.old_matrix, args.output, step=args.step), indent=2))
