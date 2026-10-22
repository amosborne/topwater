"""Prescribed-heave radiation check through the nonlinear moving-canoe stage.

This is an unqualified finite absorbing-tank experiment, not a radiation model
used by production dynamics. The short default only exercises the runner.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import RK23, trapezoid
from scipy.optimize import brentq

from topwater.hydrostatics import hydrostatics
from .moving_canoe import make_domain, initial_state, rhs, field


def heave(time, amplitude, omega, ramp_periods=2.):
    """Displacement, velocity and acceleration, with a C2 quintic startup."""
    if amplitude <= 0 or omega <= 0 or ramp_periods <= 0:
        raise ValueError('positive amplitude, frequency and ramp periods required')
    duration = ramp_periods*2*np.pi/omega
    s = np.clip(time/duration, 0., 1.)
    ramp = s**3*(10+s*(-15+6*s))
    first = 30*s*s*(1-s)**2/duration
    second = 60*s*(1-s)*(1-2*s)/duration**2
    c, s_phase = np.cos(omega*time), np.sin(omega*time)
    return amplitude*np.array([ramp*c, first*c-omega*ramp*s_phase,
        second*c-2*omega*first*s_phase-omega**2*ramp*c])


def sponge(points, strength, buffer=.5):
    """Smooth bounded outer damping, zero on the local pocket plus buffer."""
    if strength < 0 or buffer < 0:
        raise ValueError('nonnegative sponge strength and buffer required')
    xy = np.abs(points[:, :2])
    outer = xy.max(axis=0)
    inner = np.array([2*np.pi/3, np.pi/3])+buffer
    if np.any(inner >= outer):
        raise ValueError('no exterior sponge width remains after the buffer')
    s = np.clip((xy-inner)/(outer-inner), 0., 1.)
    ramp = s**3*(10+s*(-15+6*s))
    return strength*(ramp[:, 0]+ramp[:, 1]-ramp[:, 0]*ramp[:, 1]), outer-inner


def harmonic_fit(times, forces, omega, amplitude, restoring):
    """Time-weighted DC/cos/sin fit, using the exp(-i omega t) convention.

    Total heave force = DC + a*(omega^2*A-C)*cos + a*omega*B*sin.
    Trapezoidal time weights avoid bias from adaptive sampling density.
    """
    t, f = np.asarray(times), np.asarray(forces)
    if len(t) < 4 or np.any(np.diff(t) <= 0):
        raise ValueError('at least four strictly increasing sample times required')
    weights = np.empty(len(t))
    weights[0], weights[-1] = (t[1]-t[0])/2, (t[-1]-t[-2])/2
    weights[1:-1] = (t[2:]-t[:-2])/2
    design = np.column_stack((np.ones(len(t)), np.cos(omega*t), np.sin(omega*t)))
    coefficients, _, rank, _ = np.linalg.lstsq(
        design*np.sqrt(weights[:, None]), f*np.sqrt(weights), rcond=None)
    if rank != 3:
        raise ValueError('harmonic fit has deficient rank')
    residual = f-design@coefficients
    rms = np.sqrt(weights@(residual**2)/weights.sum())
    velocity = -amplitude*omega*np.sin(omega*t)
    # Integrate the large constant buoyancy term analytically. Over complete
    # periods its work is exactly zero, even with irregular accepted samples.
    work = (trapezoid((f-coefficients[0])*velocity, t)+
            coefficients[0]*amplitude*(np.cos(omega*t[-1])-np.cos(omega*t[0])))
    return dict(dc=float(coefficients[0]), cosine=float(coefficients[1]),
        sine=float(coefficients[2]), added_mass=float((coefficients[1]/amplitude+restoring)/omega**2),
        radiation_damping=float(coefficients[2]/(amplitude*omega)),
        harmonic_residual_rms=float(rms),
        mean_fluid_work=float(work/(t[-1]-t[0])))


def _window(times, force, start, stop):
    """Include exact period endpoints by interpolating accepted samples."""
    inside = (times > start) & (times < stop)
    t = np.r_[start, times[inside], stop]
    return t, np.interp(t, times, force)


def run(n=12, *, periods=.05, omega=4., amplitude=1e-4, ramp_periods=2.,
        fit_periods=3, gpu=True, material_subdivisions=2, exterior_scale=1.,
        sponge_strength=4., sponge_buffer=.5, rtol=1e-4, derivative_step=1e-4,
        max_step=.01, capytaine=False):
    if not all(np.isfinite(v) and v > 0 for v in
               (periods, omega, amplitude, ramp_periods, rtol, derivative_step, max_step)):
        raise ValueError('positive finite run settings required')
    if rtol >= 1 or fit_periods < 1 or int(fit_periods) != fit_periods:
        raise ValueError('rtol < 1 and positive integer fit periods required')
    domain = make_domain(n, gpu, material_subdivisions=material_subdivisions,
                         exterior_scale=exterior_scale)
    size = len(domain.free_nodes)
    base = initial_state(domain, False)
    damping, widths = sponge(domain.reference_points[domain.free_nodes],
                             sponge_strength, sponge_buffer)
    period = 2*np.pi/omega
    calls = 0
    cache = {}

    def evaluate(time, fluid):
        nonlocal calls
        motion = heave(time, amplitude, omega, ramp_periods)
        state = base.copy()
        state[2] += motion[0]
        state[9] = motion[1]
        state[13:] = fluid
        diagnostics = {}
        rate = rhs(domain, state, derivative_step, diagnostics=diagnostics, damping=damping)
        acceleration = np.zeros(6); acceleration[2] = motion[2]
        force = (np.asarray(diagnostics['known_pressure_load'])-
                 np.asarray(diagnostics['added_mass'])@acceleration)
        cache.update(time=time, fluid=fluid.copy(), force=force)
        calls += 1
        if calls%20==0:
            print(json.dumps(dict(progress=True,simulated_time=time,rhs_evaluations=calls)),flush=True)
        return rate[13:]

    fluid = np.zeros(2*size)
    evaluate(0., fluid)
    times = [0.]; states = [fluid.copy()]; forces = [cache['force'].copy()]
    solver = RK23(evaluate, 0., fluid, periods*period, rtol=rtol, atol=1e-9,
                  max_step=min(max_step, period/40))
    while solver.status == 'running':
        message = solver.step()
        if solver.status == 'failed':
            raise RuntimeError(message)
        # RK23's final FSAL stage normally supplies the accepted endpoint.
        if cache['time'] != solver.t or not np.array_equal(cache['fluid'], solver.y):
            evaluate(solver.t, solver.y)
        times.append(solver.t); states.append(solver.y.copy()); forces.append(cache['force'].copy())
    times, states, forces = np.asarray(times), np.asarray(states), np.asarray(forces)
    restoring = float(hydrostatics(domain.body, domain.rest).stiffness[2, 2])
    fit = None; period_fits = []
    if periods >= ramp_periods+fit_periods:
        start = (periods-fit_periods)*period
        t, f = _window(times, forces[:, 2], start, periods*period)
        fit = harmonic_fit(t, f, omega, amplitude, restoring)
        for j in range(fit_periods):
            t, f = _window(times, forces[:, 2], start+j*period, start+(j+1)*period)
            period_fits.append(harmonic_fit(t, f, omega, amplitude, restoring))
        transfer = np.array([[row['cosine'], row['sine']] for row in period_fits])
        changes = np.linalg.norm(np.diff(transfer, axis=0), axis=1)
        fit['period_to_period_force_change'] = changes.tolist()
        harmonic_amplitude = np.hypot(fit['cosine'], fit['sine'])
        fit['period_to_period_force_change_relative'] = (
            (changes/harmonic_amplitude).tolist() if harmonic_amplitude else None)
        fit['harmonic_residual_relative'] = (
            fit['harmonic_residual_rms']/harmonic_amplitude if harmonic_amplitude else None)
    k = brentq(lambda value: 9.81*value*np.tanh(4*value)-omega**2, 1e-10, omega**2/9.81+10)
    triangles = domain.reference_points[domain.free_faces]
    maximum_edge = float(np.linalg.norm(triangles-np.roll(triangles, 1, axis=1), axis=2).max())
    report = dict(qualified=False, case='prescribed-heave nonlinear radiation experiment',
        n=n, gpu=gpu, material_subdivisions=material_subdivisions, exterior_scale=exterior_scale,
        periods=periods, ramp_periods=ramp_periods, fit_periods=fit_periods, omega=omega,
        amplitude=amplitude, rtol=rtol, derivative_step=derivative_step, max_step=max_step,
        water_depth=4., wavelength=2*np.pi/k, maximum_surface_edge=maximum_edge,
        wavelength_over_maximum_edge=2*np.pi/k/maximum_edge,
        domain_half_widths=np.max(np.abs(domain.reference_points[:, :2]), axis=0).tolist(),
        sponge_widths=widths.tolist(), sponge_strength=sponge_strength, sponge_buffer=sponge_buffer,
        restoring_heave=restoring, rhs_evaluations=calls, accepted_steps=len(times)-1,
        harmonic_fit=fit, individual_period_fits=period_fits,
        scope='Finite absorbing tank; fixed contact topology. No outgoing-boundary, spatial, '
              'amplitude or transient convergence established. Short smoke runs have no harmonic fit. '
              'No performance or finite-frequency qualification claim.')
    if capytaine:
        from tests.capytaine_reference import frequency_response
        current = field(domain, base[:3], base[3:7], np.zeros(size), np.zeros(size), np.zeros(6))
        wet = current['points'][domain.local_mapping][domain.contact.faces[domain.contact.body_faces]]
        reference = frequency_response(wet[:, ::-1], base[:3], domain.body.mass,
                                       [omega], [], water_depth=4.)
        report['capytaine'] = dict(added_mass_heave=float(reference['added_mass'][0, 2, 2]),
            damping_heave=float(reference['radiation_damping'][0, 2, 2]),
            scope='Same rest wet panels and depth; unbounded horizontal reference.')
        if fit is not None:
            report['capytaine']['added_mass_relative_error']=(
                fit['added_mass']/report['capytaine']['added_mass_heave']-1)
            report['capytaine']['damping_relative_error']=(
                fit['radiation_damping']/report['capytaine']['damping_heave']-1)
    motions = np.array([heave(t, amplitude, omega, ramp_periods) for t in times])
    full_states = np.broadcast_to(base, (len(times), len(base))).copy()
    full_states[:, 2] += motions[:, 0]
    full_states[:, 9] = motions[:, 1]
    full_states[:, 13:] = states
    return report, dict(times=times, states=full_states, heave=motions, hydrodynamic_load=forces)


if __name__ == '__main__':
    import argparse
    from threadpoolctl import threadpool_limits
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=12)
    parser.add_argument('--periods', type=float, default=.05)
    parser.add_argument('--omega', type=float, default=4.)
    parser.add_argument('--amplitude', type=float, default=1e-4)
    parser.add_argument('--ramp-periods', type=float, default=2.)
    parser.add_argument('--fit-periods', type=int, default=3)
    parser.add_argument('--material-subdivisions', type=int, default=2)
    parser.add_argument('--exterior-scale', type=float, default=1.)
    parser.add_argument('--sponge-strength', type=float, default=4.)
    parser.add_argument('--sponge-buffer', type=float, default=.5)
    parser.add_argument('--rtol', type=float, default=1e-4)
    parser.add_argument('--derivative-step', type=float, default=1e-4)
    parser.add_argument('--max-step', type=float, default=.01)
    parser.add_argument('--cpu', action='store_true')
    parser.add_argument('--capytaine', action='store_true')
    args = vars(parser.parse_args()); args['gpu'] = not args.pop('cpu')
    with threadpool_limits(1):
        report, states = run(**args)
    path = Path('artifacts'); path.mkdir(exist_ok=True)
    suffix = (f"n{args['n']}-material{args['material_subdivisions']}-periods{args['periods']:g}"
        f"-omega{args['omega']:g}-amplitude{args['amplitude']:g}-extent{args['exterior_scale']:g}"
        f"-sponge{args['sponge_strength']:g}-rtol{args['rtol']:g}-step{args['max_step']:g}")
    np.savez(path/f'canoe-radiation-{suffix}.npz', **states)
    (path/f'canoe-radiation-{suffix}.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
