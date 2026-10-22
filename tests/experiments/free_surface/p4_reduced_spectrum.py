"""Offline fixed-reference cosine restriction of the saved actual canoe DNO.

One saved interior LU; at most 326 restricted Schur columns, in blocks of 16.
No geometry owner, volume assembly, full surface Schur matrix or trajectory.
This changes the resolved surface space; the incident cutoff does not justify it.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
from time import perf_counter

import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import norm, splu
from threadpoolctl import threadpool_info, threadpool_limits

from tests.benchmark_rankine import peak_rss_bytes
from .p4_fluid_ordering import fresh, load_input
from .p4_surface_basis import rectangle_cosines


class ReducedBudget(RuntimeError):
    """The single diagnostic ended; partial projected columns are retained."""


def projected_columns(kif, kff, basis, solve):
    """Actual Schur action and columnwise interior residual; no dense full Schur."""
    right = kif@basis
    extension = solve(right)
    return kff@basis-kif.T@extension, right, extension


def _sparse(saved, name, shape):
    if not np.array_equal(saved[name+'_shape'], shape):
        raise ValueError('saved surface block shape differs: '+name)
    block = csr_matrix((saved[name+'_data'], saved[name+'_indices'], saved[name+'_indptr']), shape=shape)
    block.check_format(full_check=True)
    if not block.has_canonical_format or not np.all(np.isfinite(block.data)):
        raise ValueError('finite canonical surface CSR required: '+name)
    return block


def run(stage_path, matrix_path, surface_path, output, *, seconds=120., block_size=16):
    """Bounded single-factor frozen-body diagnostic, preserving partial results."""
    if not np.isfinite(seconds) or seconds <= 0 or not 1 <= block_size <= 16:
        raise ValueError('positive budget and 1..16 RHS columns per block required')
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    report = dict(qualified=False, complete=False, status='initializing',
        scope='Frozen-body calm actual-canoe DNO restricted to fixed-reference P4 nodal cosine traces. '
              'Not a nonlinear-six-DOF trajectory, force, wave-bandwidth or throughput qualification.',
        missing_validation='Nonzero surface potential, generated harmonics, finite-pose physical bandwidth, '
                           'and force/motion convergence under surface-space enrichment.',
        budget_seconds=seconds, block_size=block_size, completed_columns=0,
        budget_scope='Checked between blocks and major operations; an individual LU/solve is not preempted.',
        body_correction='Body held fixed; no six-DOF low-rank correction claimed.',
        reference_basis='Fixed on saved reference master XY; never reevaluated at current XY.',
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), phase_seconds={},
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                     native_threadpools=threadpool_info()),
        memory_scope='Lifetime process peak RSS, not incremental memory or a factor-fill bound.')
    arrays = {}

    def checkpoint():
        report.update(elapsed_seconds=perf_counter()-start, peak_rss_bytes=peak_rss_bytes())
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        if arrays:
            np.savez(output.with_suffix('.npz'), **arrays)

    def check():
        if perf_counter()-start > seconds:
            raise ReducedBudget('wall budget exhausted between operations; no restart')

    def timed(name, function):
        print(json.dumps(dict(phase=name, event='start')), flush=True)
        before = perf_counter()
        try:
            return function()
        finally:
            report['phase_seconds'][name] = perf_counter()-before
            print(json.dumps(dict(phase=name, seconds=report['phase_seconds'][name])), flush=True)

    try:
        stage_path, matrix_path, surface_path = map(Path, (stage_path, matrix_path, surface_path))
        block_path = surface_path.with_suffix('.blocks.npz')
        source = json.loads(surface_path.read_text())
        stage_hash, block_hash = [sha256(path.read_bytes()).hexdigest() for path in (stage_path, block_path)]
        if (not source['complete'] or (source['volume_order'], source['boundary_order']) != (8, 9)
                or source['stage_sha256'] != stage_hash
                or source['surface_blocks_sha256'] != block_hash
                or source['capture_checks']['maximum_free_surface_height'] != 0.
                or source['capture_checks']['eta_map_mass_relative_difference'] > 1e-11):
            raise ValueError('complete hash-matched flat Ceta=M surface capture required')
        kii, captured, manifest = timed('verified_interior_capture', lambda: load_input(matrix_path))
        with np.load(stage_path, allow_pickle=False) as saved:
            q, reference, free = [saved[name] for name in ('q', 'reference_points', 'free_nodes')]
            elements = saved['mesh_elements']
        with np.load(block_path, allow_pickle=False) as saved:
            if (str(saved['stage_sha256']) != stage_hash
                    or str(saved['interior_archive_sha256']) != manifest['npz_sha256']
                    or source['matrix_manifest_npz_sha256'] != manifest['npz_sha256']
                    or manifest['provenance']['stage_sha256'] != stage_hash
                    or not np.array_equal(saved['q'], q) or not np.array_equal(captured['q'], q)
                    or not np.array_equal(saved['free'], free) or not np.array_equal(captured['free'], free)
                    or not np.array_equal(saved['inside'], captured['inside'])
                    or not np.array_equal(np.setdiff1d(np.unique(elements), free), captured['inside'])
                    or not np.all(q[6:] == 0.) or np.max(abs(saved['current_points'][free, 2])) != 0.):
                raise ValueError('saved stage, blocks and interior master identity differs')
            kif = _sparse(saved, 'kif', (len(captured['inside']), len(free)))
            kff = _sparse(saved, 'kff', (len(free),)*2)
            mass = _sparse(saved, 'mass', (len(free),)*2)
        if max(norm(kii-kii.T)/norm(kii), norm(kff-kff.T)/norm(kff), norm(mass-mass.T)/norm(mass)) > 1e-11:
            raise ValueError('saved energy or consistent mass symmetry failed')
        bounds = np.array([reference[:, :2].min(axis=0), reference[:, :2].max(axis=0)])
        expected = np.array([[-2*np.pi, -np.pi], [2*np.pi, np.pi]])
        if not np.allclose(bounds, expected, rtol=0, atol=2e-14):
            raise ValueError('this bounded actual-canoe diagnostic requires the original tank rectangle')
        basis, pairs = rectangle_cosines(reference[free, :2], .9, enrichment=2.)
        wavevectors = pairs*np.pi/(bounds[1]-bounds[0])
        radii = np.linalg.norm(wavevectors, axis=1)
        counts = [int(np.count_nonzero(radii <= factor*np.pi/.9*(1+32*np.finfo(float).eps)))
                  for factor in (1., 1.5, 2.)]
        if counts != [85, 189, 326] or basis.shape != (3168, 326):
            raise ValueError('actual reference-canoe cosine master dimensions differ')
        gram = basis.T@(mass@basis)
        stiffness = np.full((len(pairs),)*2, np.nan)
        arrays.update(projected_stiffness=stiffness, projected_mass=gram, mode_pairs=pairs,
                      wavevectors=wavevectors, counts=np.asarray(counts), free=free,
                      reference_free_xy=reference[free, :2], bounds=bounds)
        report.update(beam=.9, rectangle_lengths=(bounds[1]-bounds[0]).tolist(), counts=counts,
            flat_eta_map_mass_relative_difference=source['capture_checks']['eta_map_mass_relative_difference'],
            basis_source_sha256=sha256(Path(__file__).with_name('p4_surface_basis.py').read_bytes()).hexdigest(),
            input_hashes=dict(stage=stage_hash, surface_json=sha256(surface_path.read_bytes()).hexdigest(),
                              surface_blocks=block_hash, interior_json=sha256(matrix_path.read_bytes()).hexdigest(),
                              interior_archive=manifest['npz_sha256']),
            input_paths=dict(stage=str(stage_path), surface=str(surface_path), interior=str(matrix_path)),
            basis_bytes=basis.nbytes, interior_rows=kii.shape[0], surface_rows=len(free),
            column_residual_maximum=0., surface_results=[])
        checkpoint()
        check()
        factor = timed('single_interior_mmd_factor', lambda: splu(kii, permc_spec='MMD_AT_PLUS_A'))
        check()
        phase = perf_counter()
        for first in range(0, len(pairs), block_size):
            check()
            stop = min(first+block_size, len(pairs))
            action, right, extension = projected_columns(kif, kff, basis[:, first:stop], factor.solve)
            residuals = np.linalg.norm(kii@extension-right, axis=0)/np.maximum(
                np.linalg.norm(right, axis=0), np.finfo(float).tiny)
            report['column_residual_maximum'] = max(report['column_residual_maximum'], float(residuals.max()))
            if report['column_residual_maximum'] > 1e-9:
                raise ValueError('actual interior solve residual failed')
            stiffness[:, first:stop] = basis.T@action
            if first == 0:
                report['constant_gauge_relative_residual'] = float(np.linalg.norm(action[:, 0])/
                    max(np.linalg.norm(kff@basis[:, 0]), np.finfo(float).tiny))
                if report['constant_gauge_relative_residual'] > 1e-9:
                    raise ValueError('actual constant gauge failed')
            report['completed_columns'] = stop
            checkpoint()
            print(json.dumps(dict(completed_columns=stop, elapsed_seconds=perf_counter()-start)), flush=True)
        report['phase_seconds']['restricted_schur_columns'] = perf_counter()-phase
        for enrichment, count in zip((1., 1.5, 2.), counts):
            check()
            small_s, small_m = stiffness[:count, :count], gram[:count, :count]
            asymmetry = float(np.linalg.norm(small_s-small_s.T)/np.linalg.norm(small_s))
            mass_asymmetry = float(np.linalg.norm(small_m-small_m.T)/np.linalg.norm(small_m))
            mass_values = eigh(small_m, eigvals_only=True)
            threshold = 1e-12*mass_values[-1]
            rank = int(np.count_nonzero(mass_values > threshold))
            if max(asymmetry, mass_asymmetry) > 1e-9 or rank != count:
                raise ValueError('restricted stiffness symmetry or mass rank failed')
            values, vectors = eigh((small_s+small_s.T)/2, (small_m+small_m.T)/2)
            if values[0] < -1e-8*values[-1] or values[-4] <= 0:
                raise ValueError('positive top modes and nonnegative frozen energy required')
            rows = []
            for index in range(count-1, count-5, -1):
                check()
                coefficient, eigenvalue = vectors[:, index], float(values[index])
                lifted = basis[:, :count]@coefficient
                action, right, extension = projected_columns(kif, kff, lifted[:, None], factor.solve)
                action = action[:, 0]
                residual = basis[:, :count].T@(action-eigenvalue*(mass@lifted))
                scale = max(np.linalg.norm(basis[:, :count].T@action),
                            abs(eigenvalue)*np.linalg.norm(small_m@coefficient), np.finfo(float).tiny)
                relative = float(np.linalg.norm(residual)/scale)
                interior_residual = float(np.linalg.norm(kii@extension-right)/
                    max(np.linalg.norm(right), np.finfo(float).tiny))
                if relative > 1e-8 or interior_residual > 1e-9:
                    raise ValueError('independent actual projected eigen residual failed')
                omega = float(np.sqrt(9.81*eigenvalue))
                rows.append(dict(lambda_value=eigenvalue, omega=omega, period=2*np.pi/omega,
                    rk23_linear_imaginary_axis_step_ceiling=np.sqrt(3)/omega,
                    actual_projected_generalized_relative_residual=relative,
                    actual_interior_relative_residual=interior_residual,
                    full_surface_weak_residual_euclidean_norm=float(np.linalg.norm(action-eigenvalue*(mass@lifted)))))
            arrays[f'eigenvalues_{count}'] = values
            arrays[f'coefficients_{count}'] = vectors[:, -4:][:, ::-1]
            report['surface_results'].append(dict(enrichment=enrichment, count=count,
                shortest_reference_wavelength=1.8/enrichment, mass_rank=rank,
                mass_relative_rank_threshold=1e-12, mass_condition=float(mass_values[-1]/mass_values[0]),
                stiffness_relative_asymmetry=asymmetry, mass_relative_asymmetry=mass_asymmetry,
                smallest_eigenvalue=float(values[0]), modes=rows))
            checkpoint()
        report.update(complete=True, status='three nested frozen-body restricted spectra verified',
            residual_scope='Only R.T residual is constrained; full surface residual is intentionally unrestricted.')
    except ReducedBudget as error:
        report.update(status='budget exhausted; incomplete', failure=str(error))
    except BaseException as error:
        report.update(status='failed', failure=dict(type=type(error).__name__, message=str(error)))
        raise
    finally:
        checkpoint()
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('stage', 'matrix', 'surface', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=120.)
    parser.add_argument('--block-size', type=int, default=16)
    args = parser.parse_args()
    with threadpool_limits(1):
        report = run(args.stage, args.matrix, args.surface, args.output,
                     seconds=args.seconds, block_size=args.block_size)
    print(json.dumps(report, indent=2), flush=True)
    raise SystemExit(0 if report['complete'] else 2)
