"""One owned full second-order spectrum engine for two verification callers."""
from hashlib import sha256
import json
from time import perf_counter
import numpy as np
from scipy.linalg import eig, lu_factor, lu_solve
from topwater.constants import GRAVITY
from tests.benchmark_rankine import peak_rss_bytes


def checkpoint(args, report, arrays):
    np.savez(args.output.with_suffix('.npz'), **arrays)
    report['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))


def relative_residual(numerator, denominator):
    if denominator == 0:
        return 0. if numerator == 0 else None
    return float(numerator/denominator)


def run_spectrum(snapshot, nmap, jd, F4, F8, psi, args, report, *, reference_rate, arrays=None):
    mass = snapshot.patch.mass
    count = len(psi)
    arrays = {} if arrays is None else arrays
    # Eight real N² arrays cover input, LAPACK copies/work, real+complex right
    # vectors and conversion overlap. Batch and snapshot/input allowances are
    # separate. This is an array/workspace bound, not a total-process RSS cap.
    fft_cells = int(np.prod(snapshot.flux.fft_shape))
    retained_inputs = (F4.nbytes+F8.nbytes+jd.nbytes+nmap.nbytes
        +snapshot.coarse['cross'].nbytes+snapshot.flux.contact.nbytes
        +snapshot.flux.spectrum.nbytes)
    bound = 64*count*count+8*args.batch_columns*(10*fft_cells+12*count)+retained_inputs+32*1024**2
    if bound > args.maximum_bytes:
        raise MemoryError('full verification spectrum exceeds conservative array/workspace preflight')
    factor = lu_factor(jd)
    scaling = np.sqrt(mass.diagonal())
    report.update(dimension=count, verification_array_workspace_bound_bytes=bound,
        maximum_bytes=args.maximum_bytes, process_peak_rss_before=peak_rss_bytes(),
        memory_scope='Conservative NumPy/LAPACK/FFT array allowance; excludes Python/JIT/library baseline and is not a process-RSS guarantee',
        batch_columns=args.batch_columns, stage='assembling diagonal-similar D', qualified=False)

    def action(values, F=F8):
        sigma = -lu_solve(factor, nmap@values)
        return snapshot.mass_solve(snapshot.flux(values)+F@sigma)

    np.testing.assert_allclose(action(psi), reference_rate, rtol=1e-11, atol=1e-12)
    dense = np.empty((count, count))
    started = perf_counter()
    for first in range(0, count, args.batch_columns):
        last = min(first+args.batch_columns, count)
        columns = np.zeros((count, last-first))
        columns[np.arange(first, last), np.arange(last-first)] = 1/scaling[first:last]
        dense[:, first:last] = scaling[:, None]*action(columns)
        if first == 0 or last == count or first//512 != (last-1)//512:
            print(json.dumps(dict(assembled_columns=last, total=count)), flush=True)
    if not np.isfinite(dense).all():
        raise FloatingPointError('nonfinite assembled spectrum matrix')
    report['assembly_seconds'] = perf_counter()-started
    rng = np.random.default_rng(20261023)
    audit = np.column_stack((rng.standard_normal((count, 2)), psi))
    original = action(audit)
    assembled = (dense@(scaling[:, None]*audit))/scaling[:, None]
    error = np.linalg.norm(assembled-original)/np.linalg.norm(original)
    np.testing.assert_allclose(assembled, original, rtol=2e-10, atol=2e-10)
    report['matrixfree_action_relative_error'] = float(error)
    report['matrixfree_action_max_error'] = float(abs(assembled-original).max())
    report['stage'] = 'assembled; starting full nonsymmetric eigensolve'
    arrays.update(diagonal_similar_D=dense, mass_diagonal_scaling=scaling)
    checkpoint(args, report, arrays)
    print(json.dumps(dict(full_eigensolve=True, dimension=count, array_bound_bytes=bound)), flush=True)
    started = perf_counter()
    mu, vectors = eig(dense, left=False, right=True, overwrite_a=False, check_finite=False)
    report['eigensolve_seconds'] = perf_counter()-started
    if not np.isfinite(mu).all() or not np.isfinite(vectors).all():
        raise FloatingPointError('nonfinite full eigensystem')
    lam = np.sqrt(-GRAVITY*mu.astype(complex))
    arrays.update(mu=mu, lambda_plus=lam)
    report['stage'] = 'raw full eigenspectrum saved; selected diagnostics pending'
    checkpoint(args, report, arrays)
    chosen = np.argsort(lam.real)[-min(8, count):][::-1]
    infinity_norm = np.linalg.norm(dense, ord=np.inf)
    records, selected = [], []
    area = float(mass.sum())
    for index in chosen:
        vector = vectors[:, index]
        residual = dense@vector-mu[index]*vector
        physical = vector/scaling
        physical /= np.sqrt(np.vdot(physical, mass@physical).real/area)
        eta = -lam[index]*physical/GRAVITY
        sigma = -lu_solve(factor, nmap@physical)
        weak8 = snapshot.flux(physical)+F8@sigma
        weak4 = snapshot.flux(physical)+F4@sigma
        qchange = snapshot.mass_solve(weak4-weak8)
        rate8 = snapshot.mass_solve(weak8)
        norm_qchange = np.sqrt(np.vdot(qchange, mass@qchange).real/area)
        norm_rate = np.sqrt(np.vdot(rate8, mass@rate8).real/area)
        kinematic = weak8-lam[index]*(mass@eta)
        dynamic = -GRAVITY*(mass@eta)-lam[index]*(mass@physical)
        denominator = np.vdot(physical, mass@physical)
        ray4 = np.vdot(physical, weak4)/denominator
        ray8 = np.vdot(physical, weak8)/denominator
        dense_den = (infinity_norm+abs(mu[index]))*np.linalg.norm(vector)
        kin_den = np.linalg.norm(weak8)+abs(lam[index])*np.linalg.norm(mass@eta)
        dyn_den = GRAVITY*np.linalg.norm(mass@eta)+abs(lam[index])*np.linalg.norm(mass@physical)
        records.append(dict(index=int(index), mu_real=float(mu[index].real), mu_imag=float(mu[index].imag),
            lambda_plus_real=float(lam[index].real), lambda_plus_imag=float(lam[index].imag),
            dense_right_relative_residual=relative_residual(np.linalg.norm(residual), dense_den),
            original_weak_kinematic_relative_residual=relative_residual(np.linalg.norm(kinematic), kin_den),
            original_weak_dynamic_relative_residual=relative_residual(np.linalg.norm(dynamic), dyn_den),
            body_collocation_residual_max=float(abs(nmap@physical+jd@sigma).max()),
            q4_q8_action_difference_mass_norm=float(norm_qchange),
            q4_q8_action_difference_relative=relative_residual(norm_qchange, norm_rate),
            zero_denominators=[name for name, value in [('dense_residual', dense_den),
                ('weak_kinematic', kin_den), ('weak_dynamic', dyn_den), ('q_action', norm_rate)] if value == 0],
            q4_rayleigh_real=float(ray4.real), q4_rayleigh_imag=float(ray4.imag),
            q8_rayleigh_real=float(ray8.real), q8_rayleigh_imag=float(ray8.imag)))
        selected.append(physical)
    report.update(complete=True, stage='full finite-dimensional spectrum complete',
        raw_largest_lambda_plus_real=float(lam.real.max()), raw_smallest_mu_real=float(mu.real.min()),
        raw_largest_absolute_mu_imag=float(abs(mu.imag).max()),
        raw_positive_lambda_real_count=int(np.count_nonzero(lam.real > 0)),
        selected_modes=records, process_peak_rss_after=peak_rss_bytes(),
        limits=['All raw eigenvalues retained; no threshold/filter/symmetry repair',
            'Small complex growth must be judged against original residuals and q4/q8 action sensitivity',
            'No second full q4 eigensolve; modal q4 Rayleigh values are not q4 eigenvalues',
            'Full finite-dimensional spectrum is not continuum passivity, open-water or nonlinear6DOF qualification',
            'Runtime physical flux remains matrix-free; dense matrix is a verification artifact only'])
    arrays.update(mu=mu, lambda_plus=lam, selected_indices=chosen,
                  selected_physical_psi=np.column_stack(selected))
    checkpoint(args, report, arrays)

