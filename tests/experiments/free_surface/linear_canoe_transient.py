"""Same finite-tank linear control for a prescribed nonlinear canoe trajectory.

This retains the rest operator, consistent weak sponge, constant mode and zero
initial fluid state. It predicts startup as well as the harmonic tail; it does
not establish nonlinear or outgoing-wave qualification.
"""
import argparse
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from topwater.hydrostatics import hydrostatics
from .canoe_radiation import MODES, ramped_motion, sponge
from .moving_canoe import make_domain, operator_system
from .surface_damping import matrices


@dataclass(frozen=True)
class LinearTank:
    mass: object
    damping: object
    schur: np.ndarray
    coupling: np.ndarray
    added_mass: np.ndarray
    restoring: np.ndarray
    buoyancy: np.ndarray


def prepare(domain, *, sponge_strength=4., sponge_buffer=.5):
    """Eliminate only interior potentials from the identical rest operator.

    E phi = B V on interior rows gives FS flux S psi - C V and
    body potential momentum A_inf V + rho C.T psi. No mean is removed.
    """
    points, _, _, _, energy, load = operator_system(domain, domain.rest.position,
        domain.rest.orientation, np.zeros(len(domain.free_nodes)))
    free = domain.free_nodes
    active = np.flatnonzero(energy.diagonal() > 0)
    interior = np.setdiff1d(active, free)
    factor = splu(energy[interior][:, interior].tocsc())
    extension = factor.solve(energy[interior][:, free].toarray())
    response = factor.solve(load[interior])
    schur = energy[free][:, free].toarray()-energy[free][:, interior]@extension
    coupling = load[free]-energy[free][:, interior]@response
    added_mass = 1000.*load[interior].T@response
    sigma, _ = sponge(points[free], sponge_strength, sponge_buffer)
    mass, damping = matrices(points, domain.free_faces, free, sigma)
    hydro = hydrostatics(domain.body, domain.rest)
    return LinearTank(mass, damping, schur, coupling, added_mass,
                      hydro.stiffness, np.r_[hydro.force, hydro.moment])


def integrate(model, times, motion, *, rtol=1e-10, atol=1e-12):
    """Predict eta/psi and all six loads at requested times, starting at rest.

    motion(t) returns six world-axis displacement, velocity, acceleration
    vectors. Rotation displacement is the infinitesimal world rotation at CG.
    """
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0
            or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0)):
        raise ValueError('strictly increasing finite times starting at zero required')
    mass = splu(model.mass.tocsc())
    size = model.mass.shape[0]

    def rate(time, state):
        eta, psi = state[:size], state[size:]
        _, velocity, _ = motion(time)
        return np.r_[mass.solve(model.schur@psi-model.coupling@velocity-model.damping@eta),
                     -9.81*eta-mass.solve(model.damping@psi)]

    solved = solve_ivp(rate, (0., times[-1]), np.zeros(2*size), t_eval=times,
                       method='DOP853', rtol=rtol, atol=atol)
    if not solved.success:
        raise RuntimeError(solved.message)
    states = solved.y.T
    motions = np.asarray([motion(time) for time in times])
    psi_rate = np.asarray([rate(t, state)[size:] for t, state in zip(times, states)])
    # Fluid pressure load is -d/dt(A_inf V + rho C.T psi), plus
    # buoyancy and its hydrostatic displacement derivative. Rigid mass is absent.
    loads = (model.buoyancy-motions[:, 0]@model.restoring.T
             -motions[:, 2]@model.added_mass.T-1000.*psi_rate@model.coupling)
    return dict(times=times, eta=states[:, :size], psi=states[:, size:],
                hydrodynamic_load=loads, prescribed_motion=motions)


def compare(report_path, *, gpu=True):
    """Compare an existing variational NPZ; never run another nonlinear case."""
    report_path = Path(report_path)
    report = json.loads(report_path.read_text())
    if (report['formulation'] != 'variational'
            or report['integrated_fluid_variables'] != 'eta,Pi'):
        raise ValueError('comparison requires the variational prescribed runner')
    saved = np.load(report_path.with_suffix('.npz'))
    physical = saved['states']
    times = saved['times']
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0
            or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0)
            or physical.ndim != 2 or physical.shape[0] != len(times)
            or physical.shape[1] < 15 or (physical.shape[1]-13) % 2
            or not np.all(np.isfinite(physical))
            or saved['hydrodynamic_load'].shape != (len(times), 6)
            or not np.all(np.isfinite(saved['hydrodynamic_load']))):
        raise ValueError('saved times, physical eta/psi state and loads have invalid shape or values')
    domain = make_domain(report['n'], gpu, report['wet_subdivisions'], report['exterior_scale'])
    model = prepare(domain, sponge_strength=report['sponge_strength'],
                    sponge_buffer=report['sponge_buffer'])
    mode = MODES.index(report['mode'])

    def motion(time):
        result = np.zeros((3, 6))
        result[:, mode] = ramped_motion(time, report['amplitude'], report['omega'], report['ramp_periods'])
        return result

    linear = integrate(model, times, motion)
    if not np.allclose(linear['prescribed_motion'][:, :, mode], saved['prescribed_motion'],
                       rtol=1e-13, atol=1e-15):
        raise ValueError('saved prescribed motion does not match run metadata')
    size = len(domain.free_nodes)
    if physical.shape != (len(times), 13+2*size):
        raise ValueError('saved physical surface state does not match reconstructed domain')
    difference = saved['hydrodynamic_load']-linear['hydrodynamic_load']
    dynamic = linear['hydrodynamic_load']-model.buoyancy
    duration = times[-1]-times[0]
    summary = dict(qualified=False, source=str(report_path),
        scope='Same finite-tank linear startup control; not spatial/nonlinear/Capytaine qualification.',
        maximum_load_difference_by_row=np.max(abs(difference), axis=0).tolist(),
        load_difference_time_rms_by_row=np.sqrt(trapezoid(difference**2, times, axis=0)/duration).tolist(),
        linear_dynamic_load_time_rms_by_row=np.sqrt(trapezoid(dynamic**2, times, axis=0)/duration).tolist(),
        load_units=['N']*3+['N m']*3,
        surface_norm='Trapezoidal time RMS of the constant rest flat projected consistent-mass L2 norm.')
    for name, actual in (('eta', physical[:, 13:13+size]), ('psi', physical[:, 13+size:])):
        difference = actual-linear[name]
        error_squared = np.sum(difference*(model.mass@difference.T).T, axis=1)
        signal_squared = np.sum(linear[name]*(model.mass@linear[name].T).T, axis=1)
        error_integral = float(trapezoid(error_squared, times))
        signal_integral = float(trapezoid(signal_squared, times))
        summary[name+'_maximum_absolute_difference'] = float(np.max(abs(difference)))
        summary[name+'_mass_l2_time_rms_difference'] = float(np.sqrt(error_integral/duration))
        summary[name+'_mass_l2_time_rms_reference'] = float(np.sqrt(signal_integral/duration))
        summary[name+'_relative_mass_time_l2_difference'] = (
            float(np.sqrt(error_integral/signal_integral)) if signal_integral > 0 else None)
    return summary, linear


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('--cpu', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    with threadpool_limits(1):
        summary, arrays = compare(args.report, gpu=not args.cpu)
    output = args.output or args.report.with_name(args.report.stem+'-linear-control.json')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2)+'\n')
    np.savez(output.with_suffix('.npz'), **arrays)
    print(json.dumps(summary), flush=True)
