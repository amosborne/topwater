"""Selected-case transient surface diagnostic; coordinator allocation required.

One fixed accepted motion/train, three exterior near/far points. About5800 tiny
radiation solves plus one combined diffraction solve per actual train sideband;
no body integration, load modification, rendering or high band beyond96rad/s.
The600s wall budget is checked between solves. This checkpoint measures actual
surface sensitivity, not a universal mesh/convergence qualification.
"""
import argparse
import json
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np
from scipy.io import loadmat

from topwater._capytaine_surface import (
    build_hydrodynamics, combined_diffraction, convolve_radiation, cosine_kernel,
    radiation_acceleration_direct, radiation_velocity_transfer, recover_transfer,
    synthesize_spectrum, exterior_trace_points, incident_spectrum,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--physical-case', type=Path, required=True)
    parser.add_argument('--native-mat', type=Path, required=True)
    parser.add_argument('--motion', type=Path, required=True)
    parser.add_argument('--train', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    report = dict(complete=False, total_surface_ready=False,
        scope='selected accepted trajectory near/far transient diagnostic',
        raw_velocity='saved velocity_world, native modelworld CG linear velocity and world angular omega; not Euler rates')
    try:
        from threadpoolctl import threadpool_limits
        physical = json.loads(args.physical_case.read_text(encoding='utf-8-sig'))
        native = loadmat(args.native_mat)
        motion = loadmat(args.motion)
        stamp = np.asarray(motion['time']).reshape(-1)
        velocity = np.asarray(motion['velocity_world'])
        acceleration = np.asarray(motion['acceleration_world'])
        prepared = np.asarray(motion['excitation_prepared'])
        if velocity.shape != (len(stamp), 6) or acceleration.shape != velocity.shape:
            raise ValueError('Complete accepted six modal velocities and accelerations are required')
        if np.any(velocity[0] != 0) or stamp[0] != 0:
            raise ValueError('Selected checkpoint requires exact zero velocity and zero-prehistory origin')
        with np.load(args.train, allow_pickle=False) as data:
            train = SimpleNamespace(omega=data['omega_rad_s'].copy(), k=data['k_abs_m_inv'].copy(),
                                    theta=data['directions_body_rad'].copy(), C=data['C'].copy())

        def budget():
            if time.monotonic()-started > 600:
                raise TimeoutError('600 second selected-case surface-response checkpoint exceeded')

        triangles = np.asarray(physical['minimum_wetted_triangles_world_m'])
        if triangles.shape != (17, 3, 3):
            raise ValueError('This checkpoint requires unchanged17 wet triangles')
        np.testing.assert_array_equal(native['M'], physical['dry_mass_matrix_world'])
        hydro = build_hydrodynamics(triangles, native['cg'].reshape(3), native['M'], args.output/'native',
            hydrostatic_stiffness=native['K'], cb=native['cb'].reshape(3), volume=float(native['volume'].item()),
            density=float(native['rho'].item()), gravity=float(native['g'].item()),
            radiation_omega=np.array([.1, 1.]), excitation_omega=np.array([.1, 1.]),
            headings=np.array([0., np.pi/4, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]),
            checkpoint=budget, green_cache_dir=args.output/'green-cache')
        direct_xy = np.array([[6., 0.], [0., 1.25], [3.5, 0.]])
        strip_xy = np.array([[0., 1.0], [3.2, 0.]])
        trace_a, _ = exterior_trace_points(hydro, strip_xy, .025)
        trace_b, _ = exterior_trace_points(hydro, strip_xy, .05)
        xy = np.r_[direct_xy, trace_a, trace_b]
        target_xy = np.r_[direct_xy, strip_xy, strip_xy]
        omega = np.arange(1, 961)/40
        with threadpool_limits(1):
            transfer = radiation_velocity_transfer(hydro, xy, omega, checkpoint=budget, shared_green=True)
            direct_acceleration = radiation_acceleration_direct(hydro, xy)
            high_omega = np.array([32., 48., 64., 96.])
            high_transfer = radiation_velocity_transfer(hydro, xy, high_omega, checkpoint=budget, shared_green=True)
            report['infinity_candidate_vs_high_frequency'] = dict(
                omega_rad_s=high_omega.tolist(), real=high_transfer.real.tolist(),
                imaginary=high_transfer.imag.tolist(),
                imaginary_divided_by_omega=(high_transfer.imag/high_omega[:, None, None]).tolist(),
                exact_acceleration_term=direct_acceleration.tolist(),
                meaning='Derived exact endpoint candidate only; not a fitted term and not automatically used')
            heldout = np.array([.137, .687, 1.337, 2.637, 5.137, 9.937])
            heldout_transfer = radiation_velocity_transfer(hydro, xy, heldout, checkpoint=budget, shared_green=True)
            # Independent public evaluator equivalence: a real API comparison,
            # not merely comparison of algebraic expressions in this module.
            public = radiation_velocity_transfer(hydro, xy, np.array([.1, 1., 12.]), checkpoint=budget)
            shared = radiation_velocity_transfer(hydro, xy, np.array([.1, 1., 12.]), checkpoint=budget, shared_green=True)
            np.testing.assert_allclose(shared, public, rtol=1e-12, atol=1e-13)
            report['shared_green_max_error'] = float(np.max(abs(shared-public)))
            spectra = combined_diffraction(hydro, train, xy, return_forces=True, checkpoint=budget)
            incident = synthesize_spectrum(train.omega, incident_spectrum(train, target_xy), stamp)
            scattered = synthesize_spectrum(train.omega, spectra[1], stamp)
            direct_force = synthesize_spectrum(train.omega, spectra[2], stamp)
            response = {}
            recovered = {}
            for cutoff, step in [(12., .05), (12., .025), (24., .05), (24., .025)]:
                indices = np.flatnonzero((omega <= cutoff+1e-12) &
                    (abs(omega/step-np.round(omega/step)) < 1e-10))
                w = np.r_[0., omega[indices]]
                H = np.concatenate([np.zeros_like(transfer[:1]), transfer[indices]])
                candidate = cosine_kernel(w, H, stamp, direct_velocity=np.zeros((len(xy), 6)))
                # Diagnostic computation only. This local acceptance permits a
                # sensitivity measurement and does NOT open production readiness.
                candidate.accepted = True
                label = f'cutoff{cutoff:g}_step{step:g}'
                response[label] = convolve_radiation(candidate, stamp, velocity)
                recovered[label] = recover_transfer(candidate, heldout)
            baseline = response['cutoff24_step0.05']
            reference = response['cutoff24_step0.025']
            extended = baseline + acceleration @ direct_acceleration.T
            report['response'] = dict(
                max_abs_radiation_per_point_m=np.max(abs(reference), axis=0).tolist(),
                cutoff12_to24_step005_max_difference_per_point_m=np.max(abs(response['cutoff12_step0.05']-baseline), axis=0).tolist(),
                step005_to0025_cutoff24_max_difference_per_point_m=np.max(abs(baseline-reference), axis=0).tolist(),
                exact_infinity_acceleration_candidate=direct_acceleration.tolist(),
                candidate_acceleration_contribution_max_m=float(np.max(abs(extended-baseline))),
                baseline_initial_radiation_m=baseline[0].tolist(),
                incident_initial_m=incident[0].tolist(), scattered_initial_m=scattered[0].tolist(),
                total_initial_zero_acceleration_term_m=(incident[0]+scattered[0]+baseline[0]).tolist(),
                total_initial_infinity_candidate_m=(incident[0]+scattered[0]+extended[0]).tolist())
            report['field_vs_prepared_loads'] = dict(
                max_abs_difference_per_axis=np.max(abs(direct_force-prepared), axis=0).tolist(),
                max_abs_prepared_per_axis=np.max(abs(prepared), axis=0).tolist(),
                initial_difference=(direct_force[0]-prepared[0]).tolist(),
                meaning='New exact-sideband field loads versus unchanged accepted interpolated-force operator; no loads were modified')
            report['exterior_trace_offset_sensitivity'] = dict(
                target_xy_m=strip_xy.tolist(), offsets_m=[.025, .05],
                max_abs_total_difference_per_target_m=np.max(abs(
                    (scattered[:, 3:5]+reference[:, 3:5])-(scattered[:, 5:7]+reference[:, 5:7])), axis=0).tolist(),
                meaning='Same exact incident at targets; disturbance continued from two explicit exterior traces; no nonlinear-contact error bound')
            report['complex_recovery'] = {}
            for label, H in recovered.items():
                error = H-heldout_transfer
                report['complex_recovery'][label] = dict(
                    max_abs_error=float(np.max(abs(error))),
                    relative_error_per_mode=(np.max(abs(error), axis=(0, 1))/np.maximum(1e-30, np.max(abs(heldout_transfer), axis=(0, 1)))).tolist())
            np.savez(args.output/'response.npz', omega=omega, xy=xy, target_xy=target_xy, velocity_transfer=transfer,
                     heldout_omega=heldout, heldout_velocity_transfer=heldout_transfer,
                     direct_acceleration=direct_acceleration, time_s=stamp, incident=incident,
                     scattered=scattered, direct_force=direct_force,
                     **{'radiation_'+key: value for key, value in response.items()})
        report.update(complete=True, field_points_m=xy.tolist(),
            physical_limit='Selected fixed17-panel linear mean-hull reference; no nonlinear instantaneous-contact or universal error-bound claim')
    except BaseException as error:
        report['error'] = dict(type=type(error).__name__, message=str(error))
    report['wall_seconds'] = time.monotonic()-started
    (args.output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: report[key] for key in ('complete', 'total_surface_ready', 'wall_seconds')}))
    return 0 if report['complete'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
