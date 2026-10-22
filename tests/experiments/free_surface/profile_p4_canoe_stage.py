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


def run(stage, output, *, derivative_step=1e-4):
    stage = Path(stage);output = Path(output).with_suffix('.json')
    if any(output.with_suffix(suffix).exists() for suffix in ('.json', '.npz', '.prof')):
        raise FileExistsError('preserve existing profiling output; choose a fresh stem')
    output.parent.mkdir(parents=True, exist_ok=True)
    if not np.isfinite(derivative_step) or derivative_step <= 0:
        raise ValueError('positive finite derivative step required')
    with np.load(stage, allow_pickle=False) as saved:
        canonical = saved['canonical'].copy();physical = saved['physical'].copy()
        saved_rate = saved['canonical_rate'].copy()
    if canonical.ndim != 1 or len(canonical) % 2 or not np.all(np.isfinite(canonical)):
        raise ValueError('finite saved canonical state required')
    report = dict(scope='One fresh full nonlinear saved-state RHS; no trajectory/timing qualification',
                  complete=False, stage_path=str(stage), stage_sha256=sha256(stage.read_bytes()).hexdigest(),
                  derivative_step=derivative_step, phase_wall_seconds={},
                  profiling_overhead='Cold cProfile-instrumented diagnostic; Python call costs include profiler overhead',
                  runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                               native_threadpools=threadpool_info()))
    paths = ['profile_p4_canoe_stage.py', 'p4_canoe_field.py', 'p4_variational_canoe.py',
             'p4_geometry.py', 'p4_boundary_action.py', 'p4_surface.py',
             'p4_whole_geometry.py', 'p4_volume_map.py', 'p4_canoe_contact.py']
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
        result = phase('evaluate', lambda: evaluate(owner, data, canonical[size:], derivative_step=derivative_step))
        error = float(np.max(abs(result['physical_state']-physical)))
        report['saved_physical_state_max_error'] = error
        if saved_rate.shape != result['rate'].shape or not np.all(np.isfinite(saved_rate)):
            raise ValueError('finite saved canonical_rate matching current master space required')
        report['saved_canonical_rate_max_error'] = float(np.max(abs(result['rate']-saved_rate)))
        report['rate_comparison'] = 'Raw difference reported, not an accuracy gate; source/derivative controls may differ from the saved stage'
        if error > 1e-9:
            raise ValueError('fresh state differs from saved physical state')
        report.update(complete=True, weak_normal_relative=result['full_weak_normal_velocity_relative_residual'],
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
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        print(json.dumps(run(args.stage, args.output, derivative_step=args.derivative_step), indent=2))
