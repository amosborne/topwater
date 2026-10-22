"""Bounded frozen-body wave spectrum on the actual saved calm canoe.

Only surface-touching cell blocks are assembled; the exact saved interior
matrix is reused. One fluid MMD factor, no dense surface Schur matrix.
This is not the linearization of the moving six-DOF nonlinear trajectory.
No numerical work runs on import.
"""
import argparse
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import platform
from time import perf_counter

import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigsh, norm, splu
from threadpoolctl import threadpool_info, threadpool_limits

from tests.benchmark_rankine import peak_rss_bytes
from .p4_canoe_field import geometry_data
from .p4_fluid_ordering import load_input
from .p4_geometry import metrics
from .p4_surface import quadrature
from .p4_tetra import face_nodes
from .p4_whole_geometry import P4WholeCanoeGeometry


class SpectrumBudget(RuntimeError):
    """The declared diagnostic budget ended; no convergence is inferred."""


def surface_blocks(owner, points, inside, free, *, order=8):
    """Assemble K_IF/K_FF only; cells without a free node cannot contribute."""
    interior = np.full(len(points), -1, dtype=int)
    exterior = interior.copy()
    interior[inside] = np.arange(len(inside))
    exterior[free] = np.arange(len(free))
    elements = owner.mesh.elements
    selected = np.flatnonzero(np.any(exterior[elements] >= 0, axis=1))
    entries = [([], [], []), ([], [], [])]
    for element in selected:
        nodes = elements[element]
        on_free = exterior[nodes] >= 0
        _, gradient, measure, _ = metrics(points[nodes], order)
        gfree = gradient[:, on_free]
        for (rows, columns, values), mask, lookup in zip(entries,
                (~on_free, on_free), (interior, exterior)):
            local = np.einsum('qai,qbi,q->ab', gradient[:, mask], gfree, measure, optimize=True)
            rows.append(np.repeat(lookup[nodes[mask]], on_free.sum()))
            columns.append(np.tile(exterior[nodes[on_free]], mask.sum()))
            values.append(local.ravel())
    blocks = []
    for (rows, columns, values), count in zip(entries, (len(inside), len(free))):
        blocks.append(coo_matrix((np.concatenate(values),
            (np.concatenate(rows), np.concatenate(columns))), shape=(count, len(free))).tocsr())
    return *blocks, len(selected)


def localization(owner, points, mode, initial_eta_rate, mass):
    """Positive physical face integrals, rather than sliced mass products."""
    nodes = face_nodes(owner.mesh, owner.free_faces)
    lookup = np.full(len(points), -1, dtype=int)
    lookup[owner.free_nodes] = np.arange(len(owner.free_nodes))
    values, derivative, weights = quadrature(9)
    coordinates = points[nodes]
    tangent = np.einsum('qaj,fai->fqji', derivative, coordinates)
    projected = np.cross(tangent[:, :, 0], tangent[:, :, 1])[:, :, 2]
    if np.min(projected) <= 0 or np.min(weights) <= 0:
        raise ValueError('positive physical projected quadrature required')
    samples = np.einsum('qa,fa->fq', values, mode[lookup[nodes]])
    face_mass = np.einsum('q,fq,fq->f', weights, projected, samples**2)
    total = float(face_mass.sum())
    expected = float(mode@mass@mode)
    if not np.isclose(total, expected, rtol=1e-10, atol=1e-12):
        raise ValueError('positive face integrals differ from consistent mass')
    contact_nodes = owner.chart_to_global[owner.chart.contact_nodes]
    contact = np.any(np.isin(nodes, contact_nodes), axis=1)
    contact_vertices = np.unique(owner.free_faces[contact])
    adjacent = np.any(np.isin(owner.free_faces, contact_vertices), axis=1) & ~contact
    order = np.argsort(face_mass)[::-1]
    cumulative = np.cumsum(face_mass[order])/total
    selected = order[:np.searchsorted(cumulative, .9)+1]
    xy = coordinates[selected, :, :2]
    distance = np.linalg.norm(xy[:, :, None]-xy[:, None, :], axis=-1)
    distance[:, np.arange(15), np.arange(15)] = np.inf
    edge = points[owner.free_faces[selected]]
    edge_lengths = np.linalg.norm(edge-np.roll(edge, 1, axis=1), axis=-1)
    horizontal = tangent[selected, :, :, :2]
    singular = np.linalg.svd(horizontal, compute_uv=False)
    eta_norm = float(initial_eta_rate@mass@initial_eta_rate)
    projection = float(mode@mass@initial_eta_rate)
    return dict(contact_face_mass_fraction=float(face_mass[contact].sum()/total),
        adjacent_face_mass_fraction=float(face_mass[adjacent].sum()/total),
        faces_for_half_mass=int(np.searchsorted(cumulative, .5)+1),
        faces_for_ninety_percent_mass=len(selected),
        ninety_percent_faces_minimum_corner_edge=float(edge_lengths.min()),
        ninety_percent_faces_minimum_p4_node_spacing=float(distance.min()),
        ninety_percent_faces_minimum_horizontal_tangent_singular_value=float(singular.min()),
        ninety_percent_faces_maximum_tangent_condition=float((singular[..., 0]/singular[..., 1]).max()),
        tangent_convention='Rows are derivatives in the unit reference triangle coordinates.',
        initial_eta_rate_mass_norm_squared=eta_norm,
        initial_eta_rate_modal_projection=projection,
        initial_eta_rate_modal_fraction=projection**2/(expected*eta_norm) if eta_norm else None), face_mass


def run(stage_path, matrix_path, output, *, max_actions=80, setup_seconds=120., eigen_seconds=60., modes=4):
    """Capture one actual configuration and preserve partial/failure evidence."""
    if (not 20 <= max_actions <= 100 or modes not in (1, 2, 4)
            or setup_seconds <= 0 or eigen_seconds <= 0):
        raise ValueError('20..100 actions, 1/2/4 modes and positive wall budget required')
    output = Path(output).with_suffix('.json')
    block_path = output.with_suffix('.blocks.npz')
    if output.exists() or output.with_suffix('.npz').exists() or block_path.exists():
        raise FileExistsError('preserve diagnostic evidence; choose a fresh output stem')
    output.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    report = dict(qualified=False, complete=False, status='initializing',
        scope='Frozen-body current-geometry flat-surface gravity-wave diagnostic; '
            'omits actual body twist, evolving contact, full free-body restoring coupling '
            'and nonlinear terms. Not a coupled stability or throughput qualification.',
        max_schur_actions=max_actions, setup_budget_seconds=setup_seconds,
        eigen_budget_seconds=eigen_seconds,
        wall_budget_scope='Checked between operations; an individual assembly/factor/action is not preempted.',
        requested_modes=modes, source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        stage_path=str(stage_path), matrix_path=str(matrix_path), phase_seconds={},
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                     meshpy=version('meshpy'),
                     native_threadpools=threadpool_info()))
    arrays = {}
    calls = 0
    best = {}
    phase = 'capture'
    eigen_start = None

    def checkpoint():
        report.update(elapsed_seconds=perf_counter()-start, schur_actions=calls,
                      peak_rss_bytes=peak_rss_bytes())
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')

    def clock_check():
        elapsed = perf_counter()-(start if eigen_start is None else eigen_start)
        if elapsed > (setup_seconds if eigen_start is None else eigen_seconds):
            raise SpectrumBudget(('setup' if eigen_start is None else 'eigen')+' budget exhausted between operations')

    def timed(name, function):
        print(json.dumps(dict(phase=name, event='start')), flush=True)
        before = perf_counter()
        try:
            return function()
        finally:
            report['phase_seconds'][name] = perf_counter()-before
            print(json.dumps(dict(phase=name, seconds=report['phase_seconds'][name])), flush=True)

    try:
        before = perf_counter()
        stage_path = Path(stage_path)
        report['stage_sha256'] = sha256(stage_path.read_bytes()).hexdigest()
        with np.load(stage_path, allow_pickle=False) as saved:
            q, rate = saved['q'], saved['canonical_rate']
            reference = {name: saved[name] for name in ('reference_points', 'mesh_elements', 'free_nodes')}
        report['phase_seconds']['read_saved_stage'] = perf_counter()-before
        kii, captured, manifest = timed('load_verified_interior_capture', lambda: load_input(matrix_path))
        report['matrix_manifest_npz_sha256'] = manifest['npz_sha256']
        if not np.array_equal(q, captured['q']) or not np.all(q[6:] == 0):
            raise ValueError('exact saved calm configuration and cached-matrix q must match')
        provenance = manifest['provenance']
        if (provenance['stage_sha256'] != report['stage_sha256']
                or (provenance['volume_order'], provenance['boundary_order']) != (8, 9)):
            raise ValueError('cached interior must identify the same stage and primal quadrature8/9')
        owner = timed('geometry_make', P4WholeCanoeGeometry.make)
        for name, actual in (('reference_points', owner.mesh.points),
                             ('mesh_elements', owner.mesh.elements), ('free_nodes', owner.free_nodes)):
            if not np.array_equal(reference[name], actual):
                raise ValueError('saved/regenerated master identity differs: '+name)
        clock_check()
        data = timed('current_boundary_geometry', lambda: geometry_data(owner, q))
        inside, free, mass = [data[name] for name in ('inside', 'free', 'mass')]
        if not np.array_equal(inside, captured['inside']) or not np.array_equal(free, captured['free']):
            raise ValueError('current/cached matrix master ordering differs')
        kif, kff, touched = timed('surface_touching_blocks',
            lambda: surface_blocks(owner, data['state'].points, inside, free, order=data['volume_order']))
        report.update(surface_touching_cells=touched, total_cells=len(owner.mesh.elements),
            interior_matrix_scope='Exact saved calm K_II reused; no current interior or whole-volume assembly.',
            volume_order=data['volume_order'], boundary_order=data['boundary_order'])
        kfi = kif.T.tocsr()
        stage = perf_counter()
        flat = float(np.max(abs(data['state'].points[free, 2])))
        mass_difference = float(norm(data['eta_map']-mass)/norm(mass))
        symmetry = float(norm(kii-kii.T)/norm(kii))
        mass_symmetry = float(norm(mass-mass.T)/norm(mass))
        report['capture_checks'] = dict(
            maximum_free_surface_height=flat, eta_map_mass_relative_difference=mass_difference,
            fluid_matrix_relative_asymmetry=symmetry, mass_relative_asymmetry=mass_symmetry)
        if max(flat, mass_difference, symmetry, mass_symmetry) > 1e-11:
            raise ValueError('flat, consistent-mass, symmetry or cached-matrix identity check failed')
        report['phase_seconds']['capture_checks'] = perf_counter()-stage
        # Preserve the expensive missing blocks BEFORE any factor/eigen budget.
        surface_nodes = face_nodes(owner.mesh, owner.free_faces)
        contact_nodes = owner.chart_to_global[owner.chart.contact_nodes]
        block_arrays = dict(q=q, inside=inside, free=free, current_points=data['state'].points,
            free_faces=owner.free_faces, surface_nodes=surface_nodes, contact_nodes=contact_nodes,
            stage_sha256=np.asarray(report['stage_sha256']),
            interior_archive_sha256=np.asarray(manifest['npz_sha256']))
        for name, block in (('kif', kif), ('kff', kff), ('mass', mass)):
            block_arrays.update({name+'_data': block.data, name+'_indices': block.indices,
                name+'_indptr': block.indptr, name+'_shape': np.asarray(block.shape)})
        timed('save_surface_blocks', lambda: np.savez(block_path, **block_arrays))
        report.update(surface_blocks_path=str(block_path),
                      surface_blocks_sha256=sha256(block_path.read_bytes()).hexdigest())
        checkpoint()
        clock_check()
        fluid_factor = timed('interior_mmd_factor', lambda: splu(kii, permc_spec='MMD_AT_PLUS_A'))
        mass_factor = timed('consistent_mass_factor', lambda: splu(mass.tocsc(), permc_spec='MMD_AT_PLUS_A'))
        stage = perf_counter()
        initial_eta_rate = rate[6:len(q)]
        if initial_eta_rate.shape != (len(free),):
            raise ValueError('saved canonical rate has no matching physical eta-rate block')
        saved_phi, saved_psi, saved_velocity = [captured[name] for name in
            ('reference_potential', 'psi', 'reference_velocity')]
        fluid_right = data['load'][inside]@saved_velocity-kif@saved_psi
        saved_field_residual = float(np.linalg.norm(kii@saved_phi-fluid_right)/max(
            np.linalg.norm(fluid_right), np.finfo(float).tiny))
        eta_current = data['eta_factor'].solve(kfi@saved_phi+kff@saved_psi
            -(data['load'][free]+data['body_transport'])@saved_velocity)
        eta_difference = float(np.linalg.norm(eta_current-initial_eta_rate)/max(
            np.linalg.norm(initial_eta_rate), np.finfo(float).tiny))
        report['initial_rate_checks'] = dict(saved_actual_field_relative_residual=saved_field_residual,
            saved_actual_eta_rate_relative_difference=eta_difference,
            scope='Saved actual moving-body seed, checked with current blocks; no new auxiliary body solves.')
        if max(saved_field_residual, eta_difference) > 1e-9:
            raise ValueError('saved actual seed field or elevation rate differs from current blocks')
        arrays.update(free_points=data['state'].points[free], free_nodes=free,
                      initial_eta_rate=initial_eta_rate, q=q)
        report['phase_seconds']['initial_rate_checks'] = perf_counter()-stage
        clock_check()
        report['setup_seconds'] = perf_counter()-start
        eigen_start = perf_counter()

        def schur(vector):
            nonlocal calls
            clock_check()
            limit = max_actions-modes if phase == 'eigen' else max_actions
            if calls >= limit:
                raise SpectrumBudget('Schur action budget exhausted; final residual actions reserved')
            vector = np.asarray(vector).reshape(len(free))
            right = kif@vector
            extension = fluid_factor.solve(right)
            result = kff@vector-kfi@extension
            if not np.all(np.isfinite(result)):
                raise FloatingPointError('nonfinite frozen Schur action')
            calls += 1
            residual = float(np.linalg.norm(kii@extension-right)/max(np.linalg.norm(right), np.finfo(float).tiny))
            report['maximum_interior_solve_relative_residual'] = max(
                residual, report.get('maximum_interior_solve_relative_residual', 0.))
            denominator = float(vector@mass@vector)
            if denominator > 0:
                quotient = float(vector@result/denominator)
                if quotient > best.get('lambda', -np.inf):
                    best.update({'lambda': quotient, 'vector': vector.copy(), 'action': result.copy()})
            if calls % 10 == 0:
                print(json.dumps(dict(schur_actions=calls, elapsed_seconds=perf_counter()-start,
                                      maximum_rayleigh_quotient=best.get('lambda'))), flush=True)
            return result

        phase = 'checks'
        ones = np.ones(len(free))
        gauge = float(np.linalg.norm(schur(ones))/np.linalg.norm(kff@ones))
        random = np.random.default_rng(20261022)
        x, y = random.normal(size=(2, len(free)))
        sx, sy = schur(x), schur(y)
        scale = np.sqrt(abs((x@sx)*(y@sy)))
        bilinear = float(abs(x@sy-y@sx)/max(scale, np.finfo(float).tiny))
        report['schur_checks'] = dict(gauge_relative_residual=gauge,
            bilinear_relative_asymmetry=bilinear, sampled_energy_positive=bool(x@sx > 0 and y@sy > 0),
            consistent_mass_sample_positive=bool(x@mass@x > 0 and y@mass@y > 0), seed=20261022)
        if (max(gauge, bilinear, report['maximum_interior_solve_relative_residual']) > 1e-9
                or min(x@sx, y@sy, x@mass@x, y@mass@y) <= 0):
            raise ValueError('frozen Schur gauge, symmetry, solve or positivity check failed')
        checkpoint()
        phase = 'eigen'
        operator = LinearOperator((len(free),)*2, matvec=schur, dtype=np.float64)
        inverse = LinearOperator(mass.shape, matvec=mass_factor.solve, dtype=np.float64)
        initial = x.copy()
        contact = owner.chart_to_global[owner.chart.contact_nodes]
        initial[np.isin(free, contact)] += np.where(np.arange(len(free))[np.isin(free, contact)] % 2, 1., -1.)
        converged = True
        try:
            eigenvalues, vectors = eigsh(operator, k=modes, which='LA', M=mass,
                Minv=inverse, v0=initial, ncv=20, tol=1e-6, maxiter=80)
        except ArpackNoConvergence as error:
            eigenvalues, vectors = error.eigenvalues, error.eigenvectors
            converged = False
        phase = 'residuals'
        rows, face_masses = [], []
        order = np.argsort(eigenvalues)[::-1]
        vectors = vectors[:, order]
        for index, eigenvalue in enumerate(eigenvalues[order]):
            mode = vectors[:, index]
            mode /= np.sqrt(mode@mass@mode)
            action = schur(mode)
            residual = float(np.linalg.norm(action-eigenvalue*(mass@mode))/max(
                np.linalg.norm(action), abs(eigenvalue)*np.linalg.norm(mass@mode), np.finfo(float).tiny))
            if eigenvalue <= 0 or residual > 1e-5:
                raise ValueError('positive eigenvalue with verified generalized residual required')
            location, face_mass = localization(owner, data['state'].points, mode, initial_eta_rate, mass)
            omega = float(np.sqrt(9.81*eigenvalue))
            rows.append(dict(lambda_value=float(eigenvalue), omega=omega,
                period=2*np.pi/omega, rk23_linear_imaginary_axis_step_ceiling=np.sqrt(3)/omega,
                generalized_relative_residual=residual, **location))
            face_masses.append(face_mass)
        arrays.update(modes=vectors, eigenvalues=eigenvalues[order], mode_face_masses=np.asarray(face_masses))
        report.update(complete=bool(converged and len(rows) == modes),
            status='converged requested largest computed modes' if converged else 'partial ARPACK convergence',
            modes=rows, eigenvalue_order='descending; no full nonlinear-six-DOF spectrum claimed')
    except SpectrumBudget as error:
        report.update(status='budget exhausted; spectrum incomplete', failure=str(error))
    except BaseException as error:
        report.update(status='failed', failure=dict(type=type(error).__name__, message=str(error)))
        raise
    finally:
        if eigen_start is not None:
            report['eigen_seconds'] = perf_counter()-eigen_start
        if best:
            report['maximum_observed_rayleigh_quotient'] = best['lambda']
            report['rayleigh_scope'] = 'Lower bound on lambda_max of the frozen discrete map; not a converged eigenpair.'
            arrays.update(rayleigh_vector=best['vector'], rayleigh_action=best['action'])
        if arrays:
            np.savez(output.with_suffix('.npz'), **arrays)
        checkpoint()
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--matrix', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-actions', type=int, default=80)
    parser.add_argument('--setup-seconds', type=float, default=120.)
    parser.add_argument('--eigen-seconds', type=float, default=60.)
    parser.add_argument('--modes', type=int, choices=(1, 2, 4), default=4)
    args = parser.parse_args()
    with threadpool_limits(1):
        result = run(args.stage, args.matrix, args.output, max_actions=args.max_actions,
                     setup_seconds=args.setup_seconds, eigen_seconds=args.eigen_seconds, modes=args.modes)
    print(json.dumps(result, indent=2), flush=True)
    raise SystemExit(0 if result['complete'] else 2)
