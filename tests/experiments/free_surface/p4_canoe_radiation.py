"""Closed-tank full nonlinear P4 prescribed-heave startup transient.

Integrates eta/Pi; every reaction uses the current nonlinear action and three
sequential fresh fluid factorizations. This is not an outgoing/Capytaine test.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from .canoe_radiation import ramped_motion
from .p4_whole_geometry import P4WholeCanoeGeometry
from .p4_prescribed_canoe import reaction
from .p4_first_variation import _hash, _progress, _revision
from .geometry_rk23 import GeometryRK23


def run(*, amplitude=1e-4, omega=4., periods=1., ramp_periods=.25,
        max_step=.02, rtol=1e-4, atol=1e-10, derivative_step=1e-4):
    """Zero-start fluid, C2 heave startup, all six accepted reaction rows."""
    settings = (amplitude, omega, periods, ramp_periods, max_step, rtol, atol, derivative_step)
    if not all(np.isfinite(value) and value > 0 for value in settings):
        raise ValueError('positive finite motion and integration controls required')
    owner = P4WholeCanoeGeometry.make()
    _progress('actual P4 prescribed-heave geometry built')
    count = len(owner.free_nodes)
    duration = periods*2*np.pi/omega
    base = np.r_[owner.domain.rest.position, owner.domain.rest.orientation, np.zeros(6)]
    cache = {}
    calls = 0

    def evaluate(time, fluid):
        nonlocal calls
        motion = ramped_motion(time, amplitude, omega, ramp_periods)
        body = base.copy(); body[2] += motion[0]; body[9] = motion[1]
        acceleration = np.zeros(6); acceleration[2] = motion[2]
        calls += 1
        _progress(f'prescribed heave RHS{calls}, t={time:.12g}, begin')
        current = reaction(owner, body, fluid[:count], fluid[count:], acceleration,
                           step=derivative_step)
        cache.update(time=time, fluid=fluid.copy(), result=current, motion=motion.copy())
        _progress(f'prescribed heave RHS{calls}, t={time:.12g}, complete')
        return np.r_[current['eta_rate'], current['pi_rate']]

    # Only physical eta/psi history is retained; no duplicate accepted eta/Pi
    # or full canonical history. The integrator owns its current eta/Pi state.
    times, states, loads, controllers, motions = [], [], [], [], []
    energies, volumes, clearances, pushforward, normal, work = [], [], [], [], [], []

    def record(time, fluid):
        if cache.get('time') != time or not np.array_equal(cache.get('fluid'), fluid):
            evaluate(time, fluid)
        current = cache['result']
        times.append(time); states.append(current['physical_state'].copy())
        loads.append(current['hydrodynamic_load'].copy())
        controllers.append(current['controller_load'].copy()); motions.append(cache['motion'].copy())
        energies.append(current['hamiltonian']); volumes.append(current['fluid_volume'])
        clearances.append(current['minimum_rim_clearance_lower_bound'])
        pushforward.append(float(np.max(abs(current['fluid_momentum_pushforward_error']))))
        normal.append(current['full_weak_normal_velocity_relative_residual'])
        work.append(current['fluid_work_identity_absolute_residual'])

    solver = GeometryRK23(evaluate, 0., np.zeros(2*count), duration, rtol=rtol, atol=atol,
                  max_step=max_step, first_step=min(duration, max_step))
    record(0., solver.y)
    while solver.status == 'running':
        message = solver.step()
        if solver.status == 'failed':
            raise RuntimeError(message)
        record(solver.t, solver.y)
    arrays = dict(times=np.asarray(times), states=np.asarray(states),
        hydrodynamic_load=np.asarray(loads), controller_load=np.asarray(controllers),
        prescribed_motion=np.asarray(motions), hamiltonian=np.asarray(energies),
        fluid_volume=np.asarray(volumes), rim_clearance_lower_bound=np.asarray(clearances),
        fluid_momentum_pushforward_maximum_error=np.asarray(pushforward),
        weak_normal_relative_residual=np.asarray(normal, dtype=float),
        fluid_work_identity_absolute_residual=np.asarray(work),
        reference_points=owner.mesh.points, mesh_elements=owner.mesh.elements,
        free_nodes=owner.free_nodes)
    files = ['p4_canoe_radiation.py', 'p4_prescribed_canoe.py', 'p4_canoe_field.py',
        'p4_variational_canoe.py', 'p4_whole_geometry.py', 'p4_volume_map.py',
        'p4_canoe_contact.py', 'p4_geometry.py', 'p4_surface.py', 'p4_boundary_action.py',
        'p4_tetra.py', 'canoe_radiation.py', 'canonical_canoe.py', 'canoe_contact.py',
        'p4_canoe_release.py', 'geometry_rk23.py']
    report = dict(qualified=False, case='closed-tank nonlinear P4 prescribed heave',
        formulation='P4 canonical prescribed motion', integrated_fluid_variables='eta,Pi',
        initial_case='equilibrium hull, eta=psi=0, zero body velocity, C2 prescribed startup',
        revision=_revision(), source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                            for name in files},
        master_hashes=dict(reference_points_sha256=_hash(owner.mesh.points, '<f8'),
            elements_sha256=_hash(owner.mesh.elements, '<i8'), free_nodes_sha256=_hash(owner.free_nodes, '<i8')),
        rest_position=owner.domain.rest.position.tolist(),
        rest_orientation=owner.domain.rest.orientation.tolist(), eta_modes=count,
        mode='heave', amplitude=amplitude, omega=omega, periods=periods, ramp_periods=ramp_periods,
        duration=duration, max_step=max_step, rtol=rtol, atol=atol, derivative_step=derivative_step,
        boundary='closed finite tank, no sponge or PML', density=1000., gravity=9.81,
        rhs_evaluations=calls, accepted_steps=len(times)-1,
        invalid_geometry_trial_rejections=solver.invalid_geometry_rejections,
        fresh_current_fluid_factorizations_per_rhs=3,
        minimum_rim_clearance_lower_bound_m=float(min(clearances)),
        load_units=['N']*3+['N m']*3,
        scope='Actual nonlinear moving-hull heave with every P4 surface mode and all six raw '
            'controller/hydrodynamic load rows. Accepted diagnostics use exact matching states. '
            'Fixed-chart/cell and conservative rim guards remain active. Forced-system energy '
            'is not claimed conserved. No six-body-axis trajectory, temporal/spatial, outgoing, '
            'Capytaine or throughput qualification.')
    return report, arrays


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--amplitude', type=float, default=1e-4)
    parser.add_argument('--omega', type=float, default=4.)
    parser.add_argument('--periods', type=float, default=1.)
    parser.add_argument('--ramp-periods', type=float, default=.25)
    parser.add_argument('--max-step', type=float, default=.02)
    parser.add_argument('--rtol', type=float, default=1e-4)
    parser.add_argument('--atol', type=float, default=1e-10)
    parser.add_argument('--derivative-step', type=float, default=1e-4)
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-prescribed-heave.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        report, arrays = run(amplitude=args.amplitude, omega=args.omega, periods=args.periods,
            ramp_periods=args.ramp_periods, max_step=args.max_step, rtol=args.rtol,
            atol=args.atol, derivative_step=args.derivative_step)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output.with_suffix('.npz'), **arrays)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
