"""Slow canonical reduced-Luke reference for the actual nonlinear canoe.

This independently differentiates the discrete physical Hamiltonian at fixed
canonical momenta. It is a correctness reference, not a production integrator
or an energy projection. Conservation by construction does not establish the
spatial accuracy of its shared BEM/FEM approximation.
"""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from topwater.frames import rotation_matrix
from .audit_canoe_energy import gravity_moments
from .moving_canoe import make_domain, initial_state, geometry, field, rhs

RHO = 1000.


def body_map(rotvec):
    """Map local exponential-coordinate rates to world translation/angular rate."""
    x, y, z = rotvec
    skew = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    squared = float(np.dot(rotvec, rotvec))
    if squared < 1e-8:
        a = .5-squared/24+squared*squared/720
        b = 1/6-squared/120+squared*squared/5040
    else:
        angle = np.sqrt(squared)
        a = (1-np.cos(angle))/squared
        b = (angle-np.sin(angle))/(squared*angle)
    mapping = np.eye(6)
    mapping[3:, 3:] += a*skew+b*skew@skew
    return mapping


def configuration(domain, state):
    rotation = rotation_matrix(state[3:7])@rotation_matrix(domain.rest.orientation).T
    size = len(domain.free_nodes)
    return np.r_[state[:3], Rotation.from_matrix(rotation).as_rotvec(), state[13:13+size]]


def _surface_maps(domain, points, height_jacobians, body_jacobian):
    lookup = np.full(len(points), -1, dtype=int)
    lookup[domain.free_nodes] = np.arange(len(domain.free_nodes))
    local = lookup[domain.free_faces]
    triangles = points[domain.free_faces]
    normal_area = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    area = normal_area[:, 2]/2
    slope = -normal_area[:, :2]/normal_area[:, 2, None]
    element_mass = area[:, None, None]*(np.ones((3, 3))+np.eye(3))/12
    indices = (np.repeat(local, 3, axis=1).ravel(), np.tile(local, (1, 3)).ravel())
    shape = (len(domain.free_nodes), len(domain.free_nodes))
    mass = coo_matrix((element_mass.ravel(), indices), shape=shape).tocsr()
    eta_map = mass.copy()
    body_transport = np.zeros((shape[0], 6))
    for axis in (0, 1):
        transport = coo_matrix(((element_mass*slope[:, axis, None, None]).ravel(),
            indices), shape=shape).tocsr()
        eta_map -= transport@height_jacobians[axis][domain.free_nodes]
        body_transport -= transport@body_jacobian[domain.free_nodes, axis]
    return eta_map, body_transport


def _prepare(domain, q, *, psi=None, pi=None):
    orientation = (Rotation.from_rotvec(q[3:6])*
        Rotation.from_quat(domain.rest.orientation, scalar_first=True)).as_quat(scalar_first=True)
    points, height_maps, body_jacobian = geometry(domain, q[:3], orientation, q[6:])
    eta_map, body_transport = _surface_maps(domain, points, height_maps, body_jacobian)
    factor = splu(eta_map.tocsc())
    if pi is not None:
        psi = factor.solve(pi/RHO, trans='T')
    solution = field(domain, q[:3], orientation, q[6:], psi, np.zeros(6))
    rotation = rotation_matrix(orientation)
    rigid = np.zeros((6, 6))
    rigid[:3, :3] = domain.body.mass*np.eye(3)
    rigid[3:, 3:] = rotation@domain.body.inertia_cg@rotation.T
    effective_mass = rigid+solution['added_mass']
    coupling = -solution['body_surface_flux']+body_transport
    moments = gravity_moments(domain, solution['points'])
    reference = gravity_moments(domain, domain.reference_points)
    gravity = (9810*sum(np.sum(value-rest) for value, rest in zip(moments, reference))
        +domain.body.mass*9.81*(q[2]-domain.rest.position[2]))
    return dict(orientation=orientation, psi=psi, eta_map=eta_map, factor=factor,
        body_map=body_map(q[3:6]), effective_mass=effective_mass,
        coupling=coupling, solution=solution, gravity=float(gravity))


def encode(domain, state):
    """Physical pose/world twist/eta/psi to canonical q and conjugate momenta."""
    q = configuration(domain, state)
    data = _prepare(domain, q, psi=state[13+len(domain.free_nodes):])
    psi, velocity = data['psi'], state[7:13]
    pi = RHO*data['eta_map'].T@psi
    py = data['body_map'].T@(data['effective_mass']@velocity
                            +RHO*data['coupling'].T@psi)
    return np.r_[q, py, pi]


def decode(domain, canonical, *, diagnostics=None):
    """Return physical state, Hamiltonian and configuration rate dH/dp."""
    size = 6+len(domain.free_nodes)
    q, momentum = canonical[:size], canonical[size:]
    data = _prepare(domain, q, pi=momentum[6:])
    psi = data['psi']
    velocity = np.linalg.solve(data['effective_mass'],
        np.linalg.solve(data['body_map'].T, momentum[:6])-RHO*data['coupling'].T@psi)
    qdot = np.r_[np.linalg.solve(data['body_map'], velocity),
        data['factor'].solve(data['solution']['flux']-data['coupling']@velocity)]
    hamiltonian = (.5*velocity@data['effective_mass']@velocity
        +data['solution']['kinetic_energy']+data['gravity'])
    state = np.r_[q[:3], data['orientation'], velocity, q[6:], psi]
    if diagnostics is not None:
        solved=data['solution']
        diagnostics.update(points=solved['points'],
            potential=solved['potential']+solved['acceleration_potentials']@velocity,
            height_jacobians=solved['height_jacobians'],
            body_jacobian=solved['body_jacobian'],fluid_volume=solved['fluid_volume'])
    return state, float(hamiltonian), qdot


def canonical_rhs(domain, canonical, config_step=1e-5, *, progress=False):
    """Hamilton's equations; configuration derivative holds all momenta fixed."""
    if not np.isfinite(config_step) or config_step<=0:
        raise ValueError('configuration difference step must be finite and positive')
    size = len(canonical)//2
    _, _, qdot = decode(domain, canonical)
    gradient = np.empty(size)
    for i in range(size):
        direction = np.zeros_like(canonical)
        direction[i] = config_step
        plus = decode(domain, canonical+direction)[1]
        minus = decode(domain, canonical-direction)[1]
        gradient[i] = (plus-minus)/(2*config_step)
        if progress and (i+1) % 20 == 0:
            print(json.dumps(dict(progress=True, config_step=config_step,
                differentiated_coordinates=i+1, coordinates=size)), flush=True)
    return np.r_[qdot, -gradient]


def physical_rate(domain, canonical, canonical_rate, step=1e-5):
    """Differentiate the inverse canonical map along an actual canonical RHS."""
    if not np.isfinite(step) or step<=0:
        raise ValueError('inverse-map difference step must be finite and positive')
    base = decode(domain, canonical)[0]
    plus = decode(domain, canonical+step*canonical_rate)[0]
    minus = decode(domain, canonical-step*canonical_rate)[0]
    for state in (plus, minus):
        if state[3:7]@base[3:7] < 0:
            state[3:7] *= -1
    return (plus-minus)/(2*step)


def run(n=12, wet_subdivisions=2, initial_scale=10., gpu=True):
    domain = make_domain(n, gpu, wet_subdivisions)
    initial = initial_state(domain, scale=initial_scale)
    canonical = encode(domain, initial)
    recovered, energy, _ = decode(domain, canonical)
    actual = rhs(domain, initial)
    rows = []
    rates = []
    for step in (2e-5, 1e-5):
        derivative = canonical_rhs(domain, canonical, step, progress=True)
        recovered_rate = physical_rate(domain, canonical, derivative)
        rates.append(recovered_rate)
        size = 6+len(domain.free_nodes)
        # This cancellation is Hamiltonian algebra, not independent validation.
        gradient = -derivative[size:]
        canonical_work = gradient@derivative[:size]+derivative[:size]@derivative[size:]
        rows.append(dict(configuration_step=step,
            canonical_rate=derivative.tolist(),physical_rate=recovered_rate.tolist(),
            canonical_work_cancellation_watts=float(canonical_work),
            world_acceleration=recovered_rate[7:13].tolist(),
            current_rhs_world_acceleration=actual[7:13].tolist(),
            acceleration_difference_norm=float(np.linalg.norm(recovered_rate[7:13]-actual[7:13])),
            surface_potential_rate_difference_norm=float(np.linalg.norm(
                recovered_rate[13+len(domain.free_nodes):]-actual[13+len(domain.free_nodes):])),
            surface_height_rate_difference_maximum=float(np.max(np.abs(
                recovered_rate[13:13+len(domain.free_nodes)]-actual[13:13+len(domain.free_nodes)])))))
    revision = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
        text=True, check=True).stdout.strip()
    return dict(qualified=False, source_revision=revision, n=n,
        initial_canonical_state=canonical.tolist(),initial_physical_state=initial.tolist(),
        wet_subdivisions=wet_subdivisions, initial_scale=initial_scale, gpu=gpu,
        initial_hamiltonian_joules=energy,
        physical_round_trip_maximum_error=float(np.max(np.abs(recovered-initial))),
        decoded_rate_refinement_acceleration_difference=float(np.linalg.norm(rates[1][7:13]-rates[0][7:13])),
        decoded_rate_refinement_potential_difference=float(np.linalg.norm(
            rates[1][13+len(domain.free_nodes):]-rates[0][13+len(domain.free_nodes):])),
        results=rows,
        scope='Independent reduced-Luke canonical reference sharing current spatial field solver. '
            'Full configuration Hamiltonian derivative by finite differences at fixed momenta; '
            'physical rates recovered through the inverse momentum map. No trajectory, '
            'production performance, nonlinear physical accuracy, or qualification claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=12)
    parser.add_argument('--wet-subdivisions', type=int, default=2)
    parser.add_argument('--initial-scale', type=float, default=10.)
    parser.add_argument('--cpu', action='store_true')
    args = parser.parse_args()
    with threadpool_limits(1):
        report = run(args.n, args.wet_subdivisions, args.initial_scale, not args.cpu)
    suffix = f'n{args.n}-wet{args.wet_subdivisions}-initial{args.initial_scale:g}'
    Path('artifacts').mkdir(exist_ok=True)
    Path(f'artifacts/canonical-canoe-{suffix}.json').write_text(
        json.dumps(report, indent=2)+'\n', encoding='utf-8')
    summary={key:value for key,value in report.items()
        if key not in ('initial_canonical_state','initial_physical_state','results')}
    summary['results']=[{key:value for key,value in row.items()
        if key not in ('canonical_rate','physical_rate')} for row in report['results']]
    print(json.dumps(summary), flush=True)
