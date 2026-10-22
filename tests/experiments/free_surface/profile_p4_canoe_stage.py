"""ONE actual saved-state nonlinear RHS cost diagnostic; no integration.

Cold geometry construction, current prepare and full evaluate are timed
separately. cProfile retains subroutine costs including assembly/factor/shape
work. This measures the present reference, not a benchmark or accuracy waiver.
"""
import argparse
import cProfile
from hashlib import sha256
import json
from pathlib import Path
import pstats
import platform
from time import perf_counter

import numpy as np
import scipy
from threadpoolctl import threadpool_limits, threadpool_info

from .p4_canoe_field import prepare
from .p4_canoe_stage import verify_master_identity
from .p4_variational_canoe import evaluate
from .p4_whole_geometry import P4WholeCanoeGeometry
from .p4_fluid_ordering import fresh, digest, save_input


def run(stage, output, *, derivative_step=1e-4, prepare_only=False, save_fluid_matrix=None):
    stage = Path(stage);output = Path(output).with_suffix('.json')
    if any(output.with_suffix(suffix).exists() for suffix in ('.json', '.npz', '.prof')):
        raise FileExistsError('preserve existing profiling output; choose a fresh stem')
    if save_fluid_matrix is not None:
        save_fluid_matrix = fresh(save_fluid_matrix)
        if save_fluid_matrix.resolve() == output.resolve():
            raise ValueError('fluid-matrix and profile output stems must differ')
    output.parent.mkdir(parents=True, exist_ok=True)
    if not np.isfinite(derivative_step) or derivative_step <= 0:
        raise ValueError('positive finite derivative step required')
    with np.load(stage, allow_pickle=False) as saved:
        canonical = saved['canonical'].copy();physical = saved['physical'].copy()
        saved_rate = saved['canonical_rate'].copy()
    if canonical.ndim != 1 or len(canonical) % 2 or not np.all(np.isfinite(canonical)):
        raise ValueError('finite saved canonical state required')
    report = dict(scope='One fresh current prepare'+(' only' if prepare_only else ' and full nonlinear RHS')+
                  '; no trajectory/timing qualification', prepare_only=bool(prepare_only),
                  complete=False, stage_path=str(stage), stage_sha256=sha256(stage.read_bytes()).hexdigest(),
                  derivative_step=derivative_step, phase_wall_seconds={},
                  profiling_overhead='Cold cProfile-instrumented diagnostic; Python call costs include profiler overhead',
                  runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                               native_threadpools=threadpool_info()))
    paths = ['profile_p4_canoe_stage.py', 'p4_canoe_field.py', 'p4_variational_canoe.py',
             'p4_geometry.py', 'p4_boundary_action.py', 'p4_surface.py',
             'p4_whole_geometry.py', 'p4_volume_map.py', 'p4_canoe_contact.py', 'p4_fluid_ordering.py']
    report['source_sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in paths}
    profiler = cProfile.Profile();result = None;failure = None;started = perf_counter()

    def phase(name, function):
        print(json.dumps(dict(progress=True, phase=name, event='start')), flush=True)
        start = perf_counter()
        try:
            return function()
        finally:
            report['phase_wall_seconds'][name] = perf_counter()-start
            print(json.dumps(dict(progress=True, phase=name, event='finished_or_failed',
                                  seconds=report['phase_wall_seconds'][name])), flush=True)

    profiler.enable()
    try:
        owner = phase('geometry_make', P4WholeCanoeGeometry.make)
        report['master_identity'] = verify_master_identity(owner, stage)
        size = 6+len(owner.free_nodes)
        if canonical.shape != (2*size,):
            raise ValueError('saved state differs from current geometry master space')
        report.update(nodes=len(owner.mesh.points), cells=len(owner.mesh.elements), free_surface_nodes=len(owner.free_nodes))
        data = phase('prepare', lambda: prepare(owner, canonical[:size]))
        if save_fluid_matrix is not None:
            provenance = dict(stage_path=str(stage), stage_sha256=report['stage_sha256'],
                source_sha256=report['source_sha256'], master_identity=report['master_identity'],
                reference_points_sha256=digest(owner.mesh.points), elements_sha256=digest(owner.mesh.elements),
                free_nodes_sha256=digest(owner.free_nodes), configuration_sha256=digest(data['q']),
                volume_order=data['volume_order'], boundary_order=data['boundary_order'])
            def capture():
                # All six body responses already exist. Only a nonzero surface
                # trace needs one additional RHS through that existing factor.
                pbody = canonical[size:size+6]
                psi = data['eta_factor'].solve(canonical[size+6:]/1000., trans='T')
                surface_rhs = -data['energy'][data['inside']][:, data['free']]@psi
                effective = (np.linalg.solve(data['body_map'].T, pbody)
                             -1000*(data['load'][data['free']]+data['body_transport']).T@psi)
                velocity = np.linalg.solve(data['effective_mass'],
                    np.linalg.solve(data['body_map'].T, pbody)-1000*data['coupling'].T@psi)
                extra_solve = bool(np.any(surface_rhs != 0))
                surface_response = data['factor'].solve(surface_rhs) if extra_solve else np.zeros(len(data['inside']))
                potential = data['response']@velocity+surface_response
                provenance['nonzero_surface_extra_fluid_rhs_solve'] = extra_solve
                coupled = dict(rigid_mass=data['rigid_mass'], body_map=data['body_map'], body_momentum=pbody,
                    psi=psi, fluid_surface_rhs=surface_rhs, p_effective=effective,
                    reference_velocity=velocity, reference_potential=potential)
                return save_input(save_fluid_matrix,
                    data['energy'][data['inside']][:, data['inside']].tocsc(),
                    data['load'][data['inside']], data['response'], data['added_mass'], data['q'],
                    data['inside'], data['free'], provenance, coupled=coupled)
            report['fluid_matrix'] = phase('persist_existing_fluid_matrix', capture)
        if prepare_only:
            report.update(complete=True, evaluate_executed=False,
                          saved_physical_and_rate_checks='Not evaluated in prepare-only mode.')
        else:
            result = phase('evaluate', lambda: evaluate(owner, data, canonical[size:], derivative_step=derivative_step))
            error = float(np.max(abs(result['physical_state']-physical)))
            report['saved_physical_state_max_error'] = error
            if saved_rate.shape != result['rate'].shape or not np.all(np.isfinite(saved_rate)):
                raise ValueError('finite saved canonical_rate matching current master space required')
            report['saved_canonical_rate_max_error'] = float(np.max(abs(result['rate']-saved_rate)))
            report['rate_comparison'] = 'Raw difference reported, not an accuracy gate; source/derivative controls may differ from the saved stage'
            if error > 1e-9:
                raise ValueError('fresh state differs from saved physical state')
            report.update(complete=True, evaluate_executed=True,
                          weak_normal_relative=result['full_weak_normal_velocity_relative_residual'],
                          fluid_work_residual=result['fluid_work_identity_absolute_residual'])
    except BaseException as error:
        failure = error
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        profiler.disable();report['total_wall_seconds'] = perf_counter()-started
        try:
            profiler.dump_stats(str(output.with_suffix('.prof')))
            stats = pstats.Stats(profiler)
            ranked = sorted(stats.stats.items(), key=lambda entry: entry[1][3], reverse=True)[:50]
            report['top_cumulative_calls'] = [dict(file=key[0], line=key[1], function=key[2],
                primitive_calls=value[0], total_calls=value[1], self_seconds=value[2], cumulative_seconds=value[3])
                for key, value in ranked]
            if result is not None:
                np.savez(output.with_suffix('.npz'), canonical=canonical, canonical_rate=result['rate'],
                         physical=result['physical_state'], qdot=result['qdot'], saved_canonical_rate=saved_rate)
            output.write_text(json.dumps(report, indent=2)+'\n')
        except BaseException as save_error:
            # Report I/O failure without masking an already-active numeric error.
            if failure is not None:
                failure.add_note('Profiling evidence could not be saved: '+repr(save_error))
            else:
                raise
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--derivative-step', type=float, default=1e-4)
    parser.add_argument('--prepare-only', action='store_true', help='omit evaluate; no full RHS timing/check')
    parser.add_argument('--save-fluid-matrix', type=Path, help='fresh separate CSC/baseline output stem')
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        print(json.dumps(run(args.stage, args.output, derivative_step=args.derivative_step,
                             prepare_only=args.prepare_only, save_fluid_matrix=args.save_fluid_matrix), indent=2))
