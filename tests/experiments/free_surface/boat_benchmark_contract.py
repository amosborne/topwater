"""Strict performance-comparison gate; implementation hashes are provenance.

Accuracy experiments may intentionally differ. They must not produce a speedup
ratio merely by supplying an explanation for changed numerical controls.
"""
import argparse
from hashlib import sha256
import json
import math
from pathlib import Path


def array_identity(arrays):
    """Hash named numerical data, independent of NPZ timestamps/compression."""
    import numpy as np
    digest = sha256()
    for name, value in sorted(arrays.items()):
        original = np.asarray(value)
        array = np.ascontiguousarray(original)
        if array.dtype.hasobject:
            raise ValueError('benchmark data must not contain Python objects')
        digest.update(name.encode()+b'\0'+array.dtype.str.encode()+b'\0')
        digest.update(np.int64(original.ndim).tobytes())
        digest.update(np.asarray(original.shape, dtype=np.int64).tobytes())
        digest.update(array.tobytes())
    return digest.hexdigest()


def performance_comparison(reference, candidate):
    reasons = []
    for label, report in [('reference', reference), ('candidate', candidate)]:
        if report.get('complete') is not True:
            reasons.append(label+': incomplete run')
        if not report.get('benchmark_contract'):
            reasons.append(label+': missing comparison contract; no historical reconstruction')
        contract = report.get('benchmark_contract', {})
        required = {'physics', 'input_sha256', 'state_layout', 'state_count', 'error_norm', 'duration_seconds',
                    'integrator', 'tolerances', 'spatial', 'max_step_seconds', 'first_step_seconds', 'solve_controls'}
        if required-set(contract):
            reasons.append(label+': incomplete comparison contract: '+','.join(sorted(required-set(contract))))
        tolerances = contract.get('tolerances', {})
        count = contract.get('state_count')
        atol = tolerances.get('absolute_by_component', [])
        if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
            reasons.append(label+': invalid state count')
        if not isinstance(atol, list) or len(atol) != count:
            reasons.append(label+': absolute tolerances do not own every state component')
        if (not positive(tolerances.get('rtol')) or not isinstance(atol, list)
                or not all(positive(value) for value in atol) or not positive(contract.get('duration_seconds'))):
            reasons.append(label+': invalid effective tolerances or duration')
        environment = report.get('measurement_environment', {})
        if environment.get('resource_status') != 'uncontended':
            reasons.append(label+': resource load is not declared uncontended')
        if not environment.get('hardware'):
            reasons.append(label+': missing hardware identity')
        if environment.get('cache_policy') not in ['isolated cold caches', 'prewarmed current producers']:
            reasons.append(label+': startup cache state is not controlled')
        if report.get('achieved_duration_seconds') != report.get('benchmark_contract', {}).get('duration_seconds'):
            reasons.append(label+': requested/achieved duration differs')
        seconds = report.get('benchmark_seconds')
        if isinstance(seconds, bool) or not isinstance(seconds, (float, int)) or not math.isfinite(seconds) or seconds <= 0:
            reasons.append(label+': missing or invalid benchmark time')
    a, b = reference.get('benchmark_contract', {}), candidate.get('benchmark_contract', {})
    for key in sorted(set(a)|set(b)):
        if a.get(key) != b.get(key):
            reasons.append('contract differs: '+key)
    ea, eb = reference.get('measurement_environment', {}), candidate.get('measurement_environment', {})
    for key in ['hardware', 'software', 'thread_pools', 'thread_environment',
                'timing_scope', 'output_policy', 'startup_treatment', 'cache_policy']:
        if ea.get(key) != eb.get(key) or not ea.get(key):
            reasons.append('measurement environment differs or is missing: '+key)
    if reasons:
        return dict(comparable=False, reasons=reasons)
    return dict(comparable=True, speedup=reference['benchmark_seconds']/candidate['benchmark_seconds'],
        scope='Matched reported problem/numerics/environment; physical accuracy qualification remains separate')


def positive(value):
    return not isinstance(value, bool) and isinstance(value, (float, int)) and math.isfinite(value) and value > 0


def measurement_environment(cp, resource_status, description, output_policy):
    import platform
    import os
    import numba
    import scipy
    import numpy as np
    from threadpoolctl import threadpool_info
    properties = cp.cuda.runtime.getDeviceProperties(0)
    name = properties['name']
    return dict(resource_status=resource_status, resource_context=description,
        hardware=dict(cpu=platform.processor(), gpu=name.decode() if isinstance(name, bytes) else name,
                      gpu_memory_bytes=int(properties['totalGlobalMem'])),
        software=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                      numba=numba.__version__, cupy=cp.__version__, os=platform.platform(),
                      cuda_driver=cp.cuda.runtime.driverGetVersion(), cuda_runtime=cp.cuda.runtime.runtimeGetVersion()),
        thread_pools=sorted([{k:p.get(k) for k in ['user_api','internal_api','prefix','version','num_threads',
                                                   'threading_layer','architecture']} for p in threadpool_info()],
                            key=lambda p: str(sorted(p.items()))),
        thread_environment={k:os.environ.get(k, 'unset') for k in
            ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']},
        timing_scope='Configured-driver setup plus integration, including in-loop artifacts; final output/audits/archive excluded',
        output_policy=output_policy, cache_policy='uncontrolled persistent caches',
        startup_treatment='Fresh process; first-use numerical compilation/cache/context work included; imports/environment launcher excluded')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    args = parser.parse_args()
    result = performance_comparison(json.loads(args.reference.read_text()), json.loads(args.candidate.read_text()))
    print(json.dumps(result, indent=2))
    if not result['comparable']:
        raise SystemExit(2)
