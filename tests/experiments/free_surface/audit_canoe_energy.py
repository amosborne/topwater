"""Resolve the nonlinear canoe's energy derivative into physical work terms.

Individual contributions generally do not vanish: the undamped closed system's
sum must vanish. Split state directions localize exchanges without labeling a
nonzero body, surface, kinetic, or gravity contribution as an error by itself.
"""
import argparse
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import numpy as np
from threadpoolctl import threadpool_limits

from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from .canoe import body_load
from .moving_canoe import make_domain, initial_state, rhs, field, equilibrium_invariants


TERMS = ('fluid_kinetic', 'fluid_gravity', 'body_kinetic', 'body_gravity')


def kinetic_and_load(domain, solution, position, velocity):
    """Mixed identity E phi = B V + scatter_FS(flux), including rim rows."""
    mapping = domain.local_mapping
    points = solution['points'][mapping]
    panels = TrianglePanels(points[domain.contact.faces], 1)
    local = body_load(SimpleNamespace(vertices=points, faces=domain.contact.faces),
        panels, domain.contact.body_faces, position)
    load = np.zeros(len(solution['points']))
    load[mapping] = local@velocity
    load[domain.free_nodes] += solution['flux']
    return float(500*(solution['potential']@load)), load


def gravity_moments(domain, points):
    """Exact element contributions to the fluid's vertical first moment."""
    triangles = points[domain.local_mapping][domain.contact.faces]
    cross_z = np.cross(triangles[:, 1]-triangles[:, 0],
                       triangles[:, 2]-triangles[:, 0])[:, 2]
    z = triangles[:, :, 2]
    local = cross_z*(z.sum(axis=1)**2+(z*z).sum(axis=1))/48
    tets = points[domain.exterior_tets]
    volume = np.abs(np.linalg.det(tets[:, 1:]-tets[:, :1]))/6
    return local, volume*tets[:, :, 2].mean(axis=1)


def energy_terms(domain, state, reference_moments):
    size = len(domain.free_nodes)
    orientation = state[3:7]/np.linalg.norm(state[3:7])
    velocity = state[7:13]
    psi = state[13+size:]
    solution = field(domain, state[:3], orientation, state[13:13+size], psi, velocity)
    fluid_kinetic, _ = kinetic_and_load(domain, solution, state[:3], velocity)
    moments = gravity_moments(domain, solution['points'])
    fluid_gravity = 9810*sum(np.sum(current-reference)
        for current, reference in zip(moments, reference_moments))
    rotation = rotation_matrix(orientation)
    inertia = rotation@domain.body.inertia_cg@rotation.T
    body_kinetic = (.5*domain.body.mass*velocity[:3]@velocity[:3]
                    +.5*velocity[3:]@inertia@velocity[3:])
    body_gravity = domain.body.mass*9.81*(state[2]-domain.rest.position[2])
    reconstruction = (fluid_kinetic+9810*sum(np.sum(moment) for moment in moments)
                      -solution['fluid_energy'])
    return np.array([fluid_kinetic, fluid_gravity, body_kinetic, body_gravity]), (
        solution['fluid_volume']), float(reconstruction)


def run(n=12, wet_subdivisions=2, initial_scale=10., gpu=True):
    domain = make_domain(n, gpu, wet_subdivisions)
    state = initial_state(domain, scale=initial_scale)
    diagnostics = {}
    rate = rhs(domain, state, diagnostics=diagnostics)
    reference = gravity_moments(domain, domain.reference_points)
    initial, volume, reconstruction = energy_terms(domain, state, reference)
    rest_energy, _ = equilibrium_invariants(domain)
    excitation = diagnostics['mechanical_energy']-rest_energy
    if excitation <= 0:
        raise RuntimeError('initial excitation must be positive for normalization')
    size = len(domain.free_nodes)
    directions = {'full': rate}
    for name, block in [('body_pose', slice(0, 7)), ('body_velocity', slice(7, 13)),
                        ('fluid_eta', slice(13, 13+size)), ('fluid_psi', slice(13+size, None))]:
        direction = np.zeros_like(rate)
        direction[block] = rate[block]
        directions[name] = direction
    results = []
    for h in (2e-4, 1e-4):
        rows = {}
        for name, direction in directions.items():
            plus = energy_terms(domain, state+h*direction, reference)
            minus = energy_terms(domain, state-h*direction, reference)
            power = (plus[0]-minus[0])/(2*h)
            rows[name] = dict(component_power_watts=dict(zip(TERMS, power.tolist())),
                total_power_watts=float(power.sum()),
                total_power_excitation_relative_per_second=float(power.sum()/excitation),
                fluid_volume_rate_m3_per_second=(plus[1]-minus[1])/(2*h))
            print(json.dumps(dict(progress=True, step=h, direction=name,
                                 total_power_watts=float(power.sum()))), flush=True)
        summed = np.sum([list(rows[name]['component_power_watts'].values())
                         for name in rows if name != 'full'], axis=0)
        full = np.array(list(rows['full']['component_power_watts'].values()))
        results.append(dict(directional_step_seconds=h, directions=rows,
            sum_split_minus_full_component_power_watts=dict(zip(TERMS, (summed-full).tolist())),
            sum_split_minus_full_total_power_watts=float((summed-full).sum())))
    revision = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
        text=True, check=True).stdout.strip()
    return dict(qualified=False, source_revision=revision, n=n, gpu=gpu,
        wet_subdivisions=wet_subdivisions, initial_scale=initial_scale,
        initial_component_energy_relative_to_rest_geometry_joules=dict(zip(TERMS, initial.tolist())),
        initial_total_excitation_energy_joules=excitation,
        initial_component_sum_minus_independent_excitation_joules=float(initial.sum()-excitation),
        initial_fluid_volume_m3=volume,
        fluid_energy_reconstruction_error_joules=reconstruction,
        initial_hydrodynamic_body_power_watts=float(state[7:13]@diagnostics['hydrodynamic_load']),
        results=results,
        scope='Actual undamped nonlinear six-DOF stage. Fluid/body kinetic and gravity '
            'work are separately differentiated along the full RHS and four disjoint '
            'state blocks. Gravity differences are accumulated per element before '
            'summing, reducing cancellation from the constant tank baseline. '
            'Full and summed split derivatives must agree within finite-difference '
            'error; individual nonzero exchanges do not identify a defect. '
            'This localizes the unresolved total work defect; no timing claim.')


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
