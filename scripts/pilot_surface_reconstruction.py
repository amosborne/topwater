"""Bounded retained-field endpoint/linearity pilot; never a total-field success.

SOURCE ONLY until the coordinator freezes this producer and allocates the sole
numerical lease.  No body integration, MATLAB, motion modification or render.
The sparse frequency set diagnoses endpoint structure; it cannot qualify a
quadrature kernel.  Expected 127 coefficient and 10 linearity BEM solves on the
supplied unchanged wet mesh, three field points.  Hard wall budget 180 seconds.
"""
import argparse
import json
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from topwater._capytaine_surface import (
    MODES, build_hydrodynamics, combined_diffraction, radiation_velocity_transfer,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--physical-case', type=Path, required=True)
    parser.add_argument('--native-mat', type=Path, required=True,
                        help='Accepted exact M/K/cg/cb/volume/rho/g; no regenerated hydrostatics')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    physical = json.loads(args.physical_case.read_text(encoding='utf-8-sig'))
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    report = dict(complete=False, total_surface_ready=False,
                  scope='retained-source endpoint/linearity pilot only',
                  physical_case=str(args.physical_case.resolve()),
                  endpoint_and_tail_policy='UNQUALIFIED; sparse probes only')
    try:
        import capytaine as cpt
        from capytaine.bem.airy_waves import airy_waves_free_surface_elevation, froude_krylov_force
        from scipy.io import loadmat

        def budget():
            if time.monotonic()-started > 180:
                raise TimeoutError('180 second retained-field pilot budget exceeded')

        native = loadmat(args.native_mat)
        triangles = np.asarray(physical['minimum_wetted_triangles_world_m'])
        if triangles.shape != (17, 3, 3):
            raise ValueError('This bounded pilot requires the accepted unchanged 17 wet triangles')
        np.testing.assert_array_equal(native['M'], physical['dry_mass_matrix_world'])
        np.testing.assert_allclose(native['cg'].reshape(3), physical['world_cg_m'], rtol=0, atol=1e-12)
        rad_w = np.array([.025, .05, .1, .2, .4, .8, 1.6, 3.2, 6.4, 12., 24., 48., 96.])
        angles = np.array([0., np.pi/4, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
        hydro = build_hydrodynamics(triangles, native['cg'].reshape(3), native['M'], args.output/'native',
            hydrostatic_stiffness=native['K'], cb=native['cb'].reshape(3), volume=float(native['volume'].item()),
            density=float(native['rho'].item()), gravity=float(native['g'].item()),
            radiation_omega=rad_w, excitation_omega=np.array([.4, .8]), headings=angles, checkpoint=budget)
        xy = np.array([[6., 0.], [0., 6.], [7., 7.]])
        transfer = radiation_velocity_transfer(hydro, xy, checkpoint=budget)
        # Two directions, both signed-kappa branches, two active frequencies.
        C = np.array([[[0, .03+.01j, .02-.005j], [0, .002-.003j, -.004+.001j]],
                      [[0, -.01+.02j, .005+.01j], [0, .001+.002j, .003-.002j]]])
        train = SimpleNamespace(omega=np.array([0., .4, .8]), theta=np.array([.2, 1.1]), C=C)
        incident, scattered, force = combined_diffraction(hydro, train, xy, return_forces=True, checkpoint=budget)
        independent_eta = np.zeros_like(scattered)
        independent_force = np.zeros_like(force)
        independent_incident = np.zeros_like(incident)
        for iw in (1, 2):
            for j in range(2):
                for branch in range(2):
                    budget()
                    problem = cpt.DiffractionProblem(body=hydro.body, omega=float(train.omega[iw]),
                        rho=hydro.density, g=hydro.gravity, water_depth=hydro.depth,
                        wave_direction=float((train.theta[j]+branch*np.pi) % (2*np.pi)))
                    result = hydro.solver.solve(problem, keep_details=True, n_threads=1)
                    c = C[j, branch, iw]
                    independent_eta[iw] += c*hydro.solver.compute_free_surface_elevation(xy, result)
                    independent_incident[iw] += c*airy_waves_free_surface_elevation(xy, problem)
                    fk = froude_krylov_force(problem)
                    independent_force[iw] += c*np.array([result.forces[m]+fk[m] for m in MODES])
        eta_error = float(np.max(abs(scattered-independent_eta)))
        force_error = float(np.max(abs(force-independent_force)))
        np.testing.assert_allclose(scattered, independent_eta, rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(force, independent_force, rtol=1e-10, atol=1e-9)
        np.testing.assert_allclose(incident, independent_incident, rtol=0, atol=1e-14)
        if not np.isfinite(transfer).all():
            raise ValueError('Nonfinite endpoint field transfer')
        np.savez(args.output/'pilot.npz', omega=rad_w, xy=xy, velocity_transfer=transfer,
                 incident=incident, scattered=scattered, direct_force=force,
                 omega_times_velocity_transfer=rad_w[:, None, None]*transfer)
        report.update(complete=True, linearity_eta_max_error_m=eta_error,
                      linearity_force_max_error=force_error, total_surface_ready=False,
                      angular_frequency_rad_s=rad_w.tolist(), field_points_m=xy.tolist(),
                      interpretation='Endpoint structure saved for independent review; no kernel acceptance or body-motion accuracy claim')
    except BaseException as error:
        report['error'] = dict(type=type(error).__name__, message=str(error))
    report['wall_seconds'] = time.monotonic()-started
    (args.output/'pilot-report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))
    return 0 if report['complete'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
