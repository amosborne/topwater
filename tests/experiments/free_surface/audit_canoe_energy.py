"""Directional conservation audit of the actual nonlinear canoe stage.

This diagnoses spatial work consistency independently of time integration.
The reference is an undamped closed tank: physical energy and fluid volume
must have zero derivative. A nonzero resolved derivative is not an RK error.
"""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
from threadpoolctl import threadpool_limits

from topwater.frames import rotation_matrix
from topwater.hydrostatics import hydrostatics
from .moving_canoe import make_domain, initial_state, rhs, field


def equilibrium_energy(domain):
    """Independent physical hydrostatics and exact rectangular tank moment."""
    wet = hydrostatics(domain.body, domain.rest)
    tank = domain.reference_points[:domain.volume_nodes]
    lower, upper = tank.min(axis=0), tank.max(axis=0)
    volume = float(np.prod(upper-lower))
    tank_first_z = volume*(lower[2]+upper[2])/2
    fluid_first_z = tank_first_z-wet.displaced_volume*wet.center_of_buoyancy[2]
    energy = 9810*fluid_first_z+domain.body.mass*9.81*domain.rest.position[2]
    return float(energy), volume-wet.displaced_volume


def invariants(domain, state):
    """Energy and volume evaluated at a state, without another RHS evaluation."""
    size = len(domain.free_nodes)
    orientation = state[3:7]/np.linalg.norm(state[3:7])
    solution = field(domain, state[:3], orientation, state[13:13+size],
        state[13+size:], state[7:13])
    inertia = rotation_matrix(orientation)@domain.body.inertia_cg@rotation_matrix(orientation).T
    energy = (solution['fluid_energy']+domain.body.mass*9.81*state[2]
        +.5*domain.body.mass*state[7:10]@state[7:10]
        +.5*state[10:13]@inertia@state[10:13])
    return float(energy), solution['fluid_volume']


def run(n=12, wet_subdivisions=2, initial_scale=10., gpu=True):
    domain = make_domain(n, gpu, wet_subdivisions)
    state = initial_state(domain, scale=initial_scale)
    diagnostics = {}
    direction = rhs(domain, state, diagnostics=diagnostics)
    rest_energy, rest_volume = equilibrium_energy(domain)
    energy = diagnostics['mechanical_energy']
    excitation = energy-rest_energy
    if excitation <= 0:
        raise RuntimeError('initial excitation energy must be positive for normalization')
    rows = []
    for h in (2e-4, 1e-4, 5e-5):
        plus = invariants(domain, state+h*direction)
        minus = invariants(domain, state-h*direction)
        energy_rate = (plus[0]-minus[0])/(2*h)
        volume_rate = (plus[1]-minus[1])/(2*h)
        rows.append(dict(directional_step_seconds=h,
            energy_derivative_watts=energy_rate,
            energy_derivative_relative_per_second=energy_rate/excitation,
            energy_change_over_one_second_relative=energy_rate*1./excitation,
            volume_derivative_m3_per_second=volume_rate,
            volume_derivative_displacement_relative_per_second=(
                volume_rate/(domain.body.mass/1000)),
            energy_plus_joules=plus[0], energy_minus_joules=minus[0],
            volume_plus_m3=plus[1], volume_minus_m3=minus[1]))
    revision = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
        text=True, check=True).stdout.strip()
    return dict(qualified=False, source_revision=revision, n=n,
        wet_subdivisions=wet_subdivisions, gpu=gpu, initial_scale=initial_scale,
        physical_time_scale_seconds=1., density_kg_per_m3=1000., gravity_m_per_s2=9.81,
        equilibrium_mechanical_energy_joules=rest_energy,
        initial_mechanical_energy_joules=energy, initial_excitation_energy_joules=excitation,
        equilibrium_fluid_volume_m3=rest_volume,
        initial_fluid_volume_m3=diagnostics['fluid_volume'], results=rows,
        scope='Actual full nonlinear six-DOF RHS at the stated initial perturbation; '
            'undamped reflecting tank, current physical-facet ALE mesh and pressure. '
            'Central state-direction derivatives test energy/work and volume consistency '
            'without integrating time. Energy normalization uses independently integrated '
            'flat-water equilibrium, not the large constant tank gravity baseline. '
            'Finite-difference cancellation and spatial refinement must be assessed; '
            'this is a diagnostic, not a physical qualification or timing benchmark.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=12)
    parser.add_argument('--wet-subdivisions', type=int, default=2)
    parser.add_argument('--initial-scale', type=float, default=10.)
    parser.add_argument('--cpu', action='store_true')
    args = parser.parse_args()
    with threadpool_limits(1):
        report = run(args.n, args.wet_subdivisions, args.initial_scale, not args.cpu)
    destination = Path('artifacts')
    destination.mkdir(exist_ok=True)
    suffix = f'n{args.n}-wet{args.wet_subdivisions}-initial{args.initial_scale:g}'
    (destination/f'canoe-energy-audit-{suffix}.json').write_text(
        json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report), flush=True)
