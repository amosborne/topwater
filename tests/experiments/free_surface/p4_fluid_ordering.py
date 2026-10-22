"""Exact saved real-fluid CSC ordering comparison; no geometry or PDE assembly.

The baseline is the existing prepare's COLAMD factor/response. This reader
factors the identical matrix once with MMD_AT_PLUS_A, using the identical six
body-load columns. No L/U arrays are materialized for measurement. Lifetime
peak RSS is not incremental factor memory or a bound on factor fill.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
from time import perf_counter

import numpy as np
import scipy
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_info, threadpool_limits

from tests.benchmark_rankine import peak_rss_bytes


def digest(values):
    values = np.ascontiguousarray(values)
    header = f'{values.dtype.str}:{values.shape}:'.encode()
    return sha256(header+values.tobytes()).hexdigest()


def fresh(path):
    path = Path(path).with_suffix('.json')
    if path.exists() or path.with_suffix('.npz').exists():
        raise FileExistsError('preserve saved matrix/order evidence; choose a fresh stem')
    return path


def _validate(matrix, load, response, added_mass):
    if (not isinstance(matrix, csc_matrix) or matrix.shape[0] != matrix.shape[1]
            or not 0 < matrix.shape[0] <= 100000 or matrix.nnz > 12000000
            or sum(a.nbytes for a in (matrix.data, matrix.indices, matrix.indptr)) > 300000000
            or matrix.dtype != np.float64 or not matrix.has_canonical_format
            or not np.all(np.isfinite(matrix.data))):
        raise ValueError('finite canonical real CSC within 100k rows/12m nnz/300MB required')
    for name, values, shape in (('load', load, (matrix.shape[0], 6)),
                               ('response', response, (matrix.shape[0], 6)),
                               ('added_mass', added_mass, (6, 6))):
        if values.shape != shape or not np.all(np.isfinite(values)) or np.iscomplexobj(values):
            raise ValueError('invalid complete real '+name)


def save_input(path, matrix, load, response, added_mass, q, inside, free, provenance, *, coupled):
    """Persist existing prepare results, without solving or modifying them."""
    path = fresh(path)
    _validate(matrix, load, response, added_mass)
    if (inside.shape != (matrix.shape[0],) or inside.dtype.kind not in 'iu'
            or free.ndim != 1 or free.dtype.kind not in 'iu' or q.ndim != 1
            or not np.all(np.isfinite(q)) or np.any(np.diff(inside) <= 0)
            or np.any(np.diff(free) <= 0) or np.intersect1d(inside, free).size):
        raise ValueError('ordered disjoint master IDs and finite configuration required')
    arrays = dict(data=matrix.data, indices=matrix.indices, indptr=matrix.indptr,
                  shape=np.asarray(matrix.shape, dtype=np.int64), body_load=load,
                  response=response, added_mass=added_mass, q=q, inside=inside, free=free)
    shapes = dict(rigid_mass=(6, 6), body_map=(6, 6), body_momentum=(6,), psi=(len(free),),
                  fluid_surface_rhs=(len(inside),), p_effective=(6,),
                  reference_velocity=(6,), reference_potential=(len(inside),))
    if coupled.keys() != shapes.keys() or any(
            np.shape(coupled[key]) != shape or not np.all(np.isfinite(coupled[key]))
            or np.iscomplexobj(coupled[key]) for key, shape in shapes.items()):
        raise ValueError('finite current coupled-body capture with the declared shapes required')
    arrays.update(coupled)
    report = dict(scope='Exact saved current real P4 fluid operator; no new assembly/solve',
        complete=True, baseline_ordering='COLAMD (unchanged scipy splu default)', density=1000.,
        mode_order=['surge', 'sway', 'heave', 'roll', 'pitch', 'yaw'],
        rows=matrix.shape[0], nonzeros=matrix.nnz,
        array_sha256={key: digest(value) for key, value in arrays.items()}, provenance=provenance,
        baseline_relative_residual=float(np.linalg.norm(matrix@response-load)/max(np.linalg.norm(load), np.finfo(float).tiny)),
        writer_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path.with_suffix('.npz'), **arrays)
    report['npz_sha256'] = sha256(path.with_suffix('.npz').read_bytes()).hexdigest()
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return report


def load_input(path):
    """Read and verify one exact saved operator and its reference fields."""
    path = Path(path).with_suffix('.json')
    manifest = json.loads(path.read_text())
    archive = path.with_suffix('.npz')
    if archive.stat().st_size > 400000000 or sha256(archive.read_bytes()).hexdigest() != manifest['npz_sha256']:
        raise ValueError('saved CSC archive exceeds cap or differs from its manifest')
    with np.load(archive, allow_pickle=False) as saved:
        arrays = {key: saved[key] for key in manifest['array_sha256']}
    if any(digest(arrays[key]) != value for key, value in manifest['array_sha256'].items()):
        raise ValueError('saved array hash differs from its manifest')
    shape = arrays['shape']
    if (shape.shape != (2,) or not np.array_equal(shape, [manifest['rows']]*2)
            or not 0 < shape[0] <= 100000 or arrays['data'].shape != (manifest['nonzeros'],)
            or len(arrays['data']) > 12000000 or arrays['indices'].shape != arrays['data'].shape
            or arrays['indptr'].shape != (shape[0]+1,) or arrays['indptr'][0] != 0
            or arrays['indptr'][-1] != len(arrays['data']) or np.any(np.diff(arrays['indptr']) < 0)
            or np.any(arrays['indices'] < 0) or np.any(arrays['indices'] >= shape[0])):
        raise ValueError('saved CSC shape/index sanity failed before factorization')
    matrix = csc_matrix((arrays['data'], arrays['indices'], arrays['indptr']), shape=tuple(shape), copy=False)
    load, baseline, baseline_added = arrays['body_load'], arrays['response'], arrays['added_mass']
    _validate(matrix, load, baseline, baseline_added)
    if manifest['density'] != 1000. or manifest['baseline_ordering'] != 'COLAMD (unchanged scipy splu default)':
        raise ValueError('baseline density/ordering differs from the prepared operator contract')
    return matrix, arrays, manifest


def compare(path, output):
    """One MMD factor on saved input; input caps do not bound factor fill."""
    path, output = Path(path).with_suffix('.json'), fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    matrix, arrays, manifest = load_input(path)
    load, baseline, baseline_added = arrays['body_load'], arrays['response'], arrays['added_mass']
    report = dict(qualified=False, complete=False, scope='Exact real-fluid ordering diagnostic only; no throughput qualification',
        input_json_sha256=sha256(path.read_bytes()).hexdigest(), input_npz_sha256=manifest['npz_sha256'],
        input_manifest=manifest, ordering='MMD_AT_PLUS_A',
        reader_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                     threadpools=threadpool_info()),
        memory_scope='Lifetime process peak RSS, not incremental factor memory; input caps do not bound LU fill.',
        peak_rss_bytes_before=peak_rss_bytes())
    factor = None
    failure = None
    try:
        start = perf_counter()
        try:
            factor = splu(matrix, permc_spec='MMD_AT_PLUS_A')
        finally:
            report['factor_wall_seconds'] = perf_counter()-start
        start = perf_counter()
        try:
            response = factor.solve(load)
        finally:
            report['solve_wall_seconds'] = perf_counter()-start
        del factor
        factor = None
        if not np.all(np.isfinite(response)):
            raise FloatingPointError('nonfinite saved-matrix solve')
        added = 1000*load.T@response
        residual = matrix@response-load
        report.update(complete=True,
            relative_residual=float(np.linalg.norm(residual)/max(np.linalg.norm(load), np.finfo(float).tiny)),
            response_relative_difference=float(np.linalg.norm(response-baseline)/max(np.linalg.norm(baseline), np.finfo(float).tiny)),
            response_maximum_absolute_difference=float(np.max(abs(response-baseline))),
            added_mass=added.tolist(), added_mass_difference=(added-baseline_added).tolist(),
            added_mass_maximum_absolute_difference=float(np.max(abs(added-baseline_added))))
        np.savez(output.with_suffix('.npz'), response=response, added_mass=added,
                 response_difference=response-baseline, added_mass_difference=added-baseline_added)
    except BaseException as error:
        failure = error
        report['complete'] = False
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        del factor
        try:
            report['peak_rss_bytes_after'] = peak_rss_bytes()
            output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        except BaseException as save_error:
            if failure is not None:
                failure.add_note('Ordering evidence could not be saved: '+repr(save_error))
            else:
                raise
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        print(json.dumps(compare(args.input, args.output), indent=2))
