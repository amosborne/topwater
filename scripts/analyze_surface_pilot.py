"""Diagnose retained surface transfers; optional bounded direct refinement.

Run only under coordinator numerical allocation. The default mode performs
stored-data analysis, no BEM. --refine performs about3000 tiny retained-source
solves at the existing three exterior points, bounded by300s between solves.
No body trajectory, hydrodynamic force or startup definition is changed.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.interpolate import interp1d

from topwater._capytaine_surface import (
    build_hydrodynamics, cosine_kernel, radiation_velocity_transfer, recover_transfer,
)


def endpoint_fit(omega, transfer, indices):
    """Diagnostic real-coefficient D+i*A/w+B/w^2 fit; not a tail prescription."""
    w = omega[indices]
    basis = np.stack([np.ones_like(w), 1j/w, 1/w**2], axis=1)
    design = np.r_[basis.real, basis.imag]
    target = transfer[indices].reshape(len(indices), -1)
    coefficients = np.linalg.lstsq(design, np.r_[target.real, target.imag], rcond=None)[0]
    prediction = basis @ coefficients
    scale = max(1e-30, float(np.max(abs(target))))
    return dict(real_coefficients_D_A_B=coefficients.reshape(3, *transfer.shape[1:]).tolist(),
                fitted_omega_rad_s=w.tolist(), maximum_residual=float(np.max(abs(prediction-target))),
                relative_maximum_residual=float(np.max(abs(prediction-target))/scale))


def analyze(omega, transfer, heldout_omega, heldout_transfer, *, dense):
    report = dict(total_surface_ready=False, dense_direct_frequency_data=bool(dense),
        convention='exp(-iwt); Hv=i*eta_per_unit_displacement/omega',
        assumptions=['Diagnostic zero low-frequency endpoint; not yet an accepted endpoint policy',
                     'Finite-band cosine inversion imposes causality; complex recovery is measured independently',
                     'Sparse interpolation does not qualify frequency quadrature or poles'])
    report['endpoint_values'] = dict(omega_rad_s=omega.tolist(),
        max_abs_transfer_over_points_and_modes=np.max(abs(transfer), axis=(1, 2)).tolist(),
        max_abs_real=np.max(abs(transfer.real), axis=(1, 2)).tolist(),
        max_abs_imaginary=np.max(abs(transfer.imag), axis=(1, 2)).tolist(),
        max_abs_omega_times_transfer=np.max(abs(omega[:, None, None]*transfer), axis=(1, 2)).tolist())
    report['high_frequency_fits'] = [endpoint_fit(omega, transfer, np.arange(len(omega)-n, len(omega))) for n in (3, 4)]
    report['low_frequency_values'] = dict(real=transfer[:4].real.tolist(), imaginary=transfer[:4].imag.tolist())
    lag = np.arange(3001)*.02
    # Same data interpolant in diagnostic mode; in dense mode every .05-grid
    # node below24 is a direct source solve. No amplitudes/phases are fitted.
    f = interp1d(np.r_[0., omega], np.concatenate([np.zeros_like(transfer[:1]), transfer]), axis=0,
                 kind='linear', bounds_error=True)
    direct = np.zeros(transfer.shape[1:])
    candidates, fields = {}, {}
    for cutoff, spacing in [(12., .1), (12., .05), (24., .1), (24., .05)]:
        w = np.arange(round(cutoff/spacing)+1)*spacing
        candidate = cosine_kernel(w, f(w), lag, direct_velocity=direct)
        recovered = recover_transfer(candidate, heldout_omega)
        error = recovered-heldout_transfer
        per_mode_scale = np.max(abs(heldout_transfer), axis=(0, 1))
        relative = np.max(abs(error), axis=(0, 1))/np.maximum(per_mode_scale, 1e-30)
        label = f'cutoff{cutoff:g}_step{spacing:g}'
        candidates[label] = dict(max_abs_complex_error=float(np.max(abs(error))),
            max_abs_real_error=float(np.max(abs(error.real))),
            max_abs_imaginary_error=float(np.max(abs(error.imag))),
            relative_error_per_mode=relative.tolist(),
            max_abs_kernel=float(np.max(abs(candidate.values))),
            kernel_value_at_zero=candidate.values[0].tolist(),
            heldout_omega_rad_s=heldout_omega.tolist(),
            recovered_real=recovered.real.tolist(), recovered_imaginary=recovered.imag.tolist())
        fields[label] = candidate.values
    report['finite_band_complex_recovery'] = candidates
    report['kernel_sensitivities'] = {
        'cutoff12_to24_step005_max_difference': float(np.max(abs(fields['cutoff12_step0.05']-fields['cutoff24_step0.05']))),
        'step01_to005_cutoff24_max_difference': float(np.max(abs(fields['cutoff24_step0.1']-fields['cutoff24_step0.05']))),
    }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--refine', action='store_true')
    parser.add_argument('--physical-case', type=Path)
    parser.add_argument('--native-mat', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    report = dict(complete=False, total_surface_ready=False)
    try:
        from threadpoolctl import threadpool_limits
        with np.load(args.input, allow_pickle=False) as data:
            omega = data['omega'].copy(); xy = data['xy'].copy()
            transfer = data['velocity_transfer'].copy()
        if args.refine:
            if args.physical_case is None or args.native_mat is None:
                raise ValueError('--refine requires --physical-case and --native-mat')
            from scipy.io import loadmat
            physical = json.loads(args.physical_case.read_text(encoding='utf-8-sig'))
            native = loadmat(args.native_mat)

            def budget():
                if time.monotonic()-started > 300:
                    raise TimeoutError('300 second source-refinement budget exceeded')

            np.testing.assert_array_equal(native['M'], physical['dry_mass_matrix_world'])
            triangles = np.asarray(physical['minimum_wetted_triangles_world_m'])
            if triangles.shape != (17, 3, 3):
                raise ValueError('The bounded refinement requires unchanged17 wet triangles')
            hydro = build_hydrodynamics(triangles, native['cg'].reshape(3), native['M'], args.output/'native',
                hydrostatic_stiffness=native['K'], cb=native['cb'].reshape(3), volume=float(native['volume'].item()),
                density=float(native['rho'].item()), gravity=float(native['g'].item()),
                radiation_omega=np.array([.1, 1.]), excitation_omega=np.array([.1, 1.]),
                headings=np.array([0., np.pi/4, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]),
                checkpoint=budget, green_cache_dir=args.output/'green-cache')
            omega = np.unique(np.r_[.0125, .025, np.arange(1, 481)/20, 32., 48., 64., 96.])
            heldout_omega = np.array([.137, .687, 1.337, 2.637, 5.137, 9.937])
            with threadpool_limits(1):
                transfer = radiation_velocity_transfer(hydro, xy, omega, checkpoint=budget)
                heldout_transfer = radiation_velocity_transfer(hydro, xy, heldout_omega, checkpoint=budget)
            np.savez(args.output/'refined-transfer.npz', omega=omega, xy=xy, velocity_transfer=transfer,
                     heldout_omega=heldout_omega, heldout_velocity_transfer=heldout_transfer)
        else:
            # Leave every second original frequency out of the interpolant.
            # This exposes interpolation error; it is not a new direct solve.
            heldout_indices = np.arange(1, len(omega)-1, 2)
            heldout_omega = omega[heldout_indices]
            heldout_transfer = transfer[heldout_indices]
            retained = np.ones(len(omega), dtype=bool); retained[heldout_indices] = False
            omega = omega[retained]; transfer = transfer[retained]
        with threadpool_limits(1):
            report.update(analyze(omega, transfer, heldout_omega, heldout_transfer, dense=args.refine))
        report['complete'] = True
    except BaseException as error:
        report['error'] = dict(type=type(error).__name__, message=str(error))
    report['wall_seconds'] = time.monotonic()-started
    (args.output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: report[key] for key in ('complete', 'total_surface_ready', 'wall_seconds')}))
    return 0 if report['complete'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
