"""Independent physical energy/volume work along a SAVED actual P4 stage RHS.

Build the geometry once; every perturbed canonical state gets a fresh current
fluid operator and solve. This does not recompute or certify the saved RHS.
New seeds carry exact master identities; legacy seeds require deterministic
meshing and an independent fresh decoded physical-state match.
No trajectory, Capytaine, spatial-accuracy or performance qualification follows.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np
from threadpoolctl import threadpool_limits

from topwater.frames import rotation_matrix
from .p4_whole_geometry import P4WholeCanoeGeometry
from .p4_canoe_field import prepare
from .p4_variational_canoe import decode
from .p4_canoe_stage import verify_master_identity
from .p4_geometry import quadrature as cell_quadrature
from .p4_surface import quadrature as face_quadrature
from .p4_tetra import face_nodes


def _physical_state(owner, canonical, label):
    """Release the current stiffness/factor after retaining physical fields."""
    size = len(canonical)//2
    data = prepare(owner, canonical[:size])
    print(json.dumps(dict(progress=True, stage='current operator prepared', state=label)), flush=True)
    decoded = decode(data, canonical[size:])
    return dict(points=data['state'].points, potential=decoded['potential'],
                physical=decoded['physical_state'], hamiltonian=decoded['hamiltonian'])


def _kinetic_difference(owner, plus, minus=None, *, order=10):
    """Integrate physical velocity squared, subtracting quadrature densities.

    Reference potential derivatives are transformed directly by the current
    geometry inverse; no assembled stiffness or reduced-energy formula is used.
    Small batches avoid all-cell quadrature tensors. None denotes zero KE.
    """
    _, derivative, weights = cell_quadrature(order)
    result = 0.
    for start in range(0, len(owner.mesh.elements), 16):
        nodes = owner.mesh.elements[start:start+16]
        densities = []
        for state in (plus, minus):
            if state is None:
                densities.append(np.zeros((len(nodes), len(weights))))
                continue
            coordinates = state['points'][nodes]
            jacobian = np.einsum('qaj,tai->tqji', derivative, coordinates)
            determinant = np.linalg.det(jacobian)
            if np.any(determinant <= 0) or not np.all(np.isfinite(determinant)):
                raise ValueError('nonpositive sampled P4 audit geometry Jacobian')
            reference_gradient = np.einsum('qaj,ta->tqj', derivative, state['potential'][nodes])
            velocity = np.einsum('tqij,tqj->tqi', np.linalg.inv(jacobian), reference_gradient)
            densities.append(500*determinant*np.einsum('tqi,tqi->tq', velocity, velocity))
        result += np.einsum('q,tq->', weights, densities[0]-densities[1])
    return float(result)


def _boundary_difference(owner, plus, minus, *, order=12):
    """Independent closed-boundary V/G; difference integrands before totals."""
    values, derivative, weights = face_quadrature(order)
    nodes = face_nodes(owner.mesh, owner.faces)
    volume = gravity = 0.
    for start in range(0, len(nodes), 32):
        terms = []
        for points in (plus, minus):
            coordinates = points[nodes[start:start+32]]
            tangent = np.einsum('qaj,fai->fqji', derivative, coordinates)
            normal_z = np.cross(tangent[:, :, 0], tangent[:, :, 1])[:, :, 2]
            height = np.einsum('qa,fa->fq', values, coordinates[:, :, 2])
            terms.append((height*normal_z, 4905*height*height*normal_z))
        volume += np.einsum('q,fq->', weights, terms[0][0]-terms[1][0])
        gravity += np.einsum('q,fq->', weights, terms[0][1]-terms[1][1])
    return float(volume), float(gravity)


def _rigid_kinetic(owner, physical):
    rotation = rotation_matrix(physical[3:7])
    velocity, omega = physical[7:10], physical[10:13]
    inertia = rotation@owner.domain.body.inertia_cg@rotation.T
    return float(.5*owner.domain.body.mass*(velocity@velocity)+.5*omega@inertia@omega)


def _source_revision():
    try:
        return subprocess.check_output(['git', 'rev-parse', 'HEAD'],
            cwd=Path(__file__).resolve().parents[3], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def audit(stage_path):
    """Run the five fresh PDEs; always return an unqualified diagnostic report."""
    stage_path = Path(stage_path)
    with np.load(stage_path) as saved:
        canonical = np.array(saved['canonical'], dtype=float)
        rate = np.array(saved['canonical_rate'], dtype=float)
        physical = np.array(saved['physical'], dtype=float)
        q = np.array(saved['q'], dtype=float)
    if (canonical.ndim != 1 or len(canonical) % 2 or rate.shape != canonical.shape
            or q.shape != (len(canonical)//2,) or not np.array_equal(q, canonical[:len(q)])
            or physical.shape != (13+2*(len(q)-6),)
            or not all(np.all(np.isfinite(value)) for value in (canonical, rate, q, physical))):
        raise ValueError('incompatible finite saved canonical state/rate and physical stage state')
    metadata_path = stage_path.with_suffix('.json')
    source_metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else None
    owner = P4WholeCanoeGeometry.make()
    master_identity = verify_master_identity(owner, stage_path)
    size = 6+len(owner.free_nodes)
    if canonical.shape != (2*size,):
        raise ValueError('saved canonical dimensions differ from rebuilt master space')
    print(json.dumps(dict(progress=True, stage='single geometry map ready',
                          active_nodes=len(np.unique(owner.mesh.elements)))), flush=True)
    base = _physical_state(owner, canonical, 'base')
    base_errors = {}
    for label, section in (('position', slice(0, 3)), ('quaternion', slice(3, 7)),
                           ('world_velocity', slice(7, 13)),
                           ('eta', slice(13, 13+len(owner.free_nodes))),
                           ('psi', slice(13+len(owner.free_nodes), None))):
        base_errors[label] = float(np.max(abs(base['physical'][section]-physical[section])))
    if max(base_errors.values()) > 1e-9:
        raise ValueError(f'fresh physical base does not match saved master ordering/state: {base_errors}')
    fluid_kinetic = _kinetic_difference(owner, base)
    _, gravity_relative = _boundary_difference(owner, base['points'], owner.mesh.points)
    rigid_kinetic = _rigid_kinetic(owner, base['physical'])
    body_gravity = float(owner.domain.body.mass*9.81*(base['physical'][2]-owner.domain.rest.position[2]))
    base_energy = float(fluid_kinetic+rigid_kinetic+gravity_relative+body_gravity)
    energy_tolerance = max(1e-7, 2e-6*abs(base_energy)/1.)  # 2ppm excitation per 1s, with W floor.
    volume_tolerance = 1e-8
    results = []
    for step in (2e-4, 1e-4):
        plus = _physical_state(owner, canonical+step*rate, f'plus h={step:g}')
        minus = _physical_state(owner, canonical-step*rate, f'minus h={step:g}')
        fluid_difference = _kinetic_difference(owner, plus, minus)
        rigid_difference = _rigid_kinetic(owner, plus['physical'])-_rigid_kinetic(owner, minus['physical'])
        volume_difference, gravity_difference = _boundary_difference(owner, plus['points'], minus['points'])
        body_difference = float(owner.domain.body.mass*9.81*(plus['physical'][2]-minus['physical'][2]))
        energy_rate = (fluid_difference+rigid_difference+gravity_difference+body_difference)/(2*step)
        primal_rate = (plus['hamiltonian']-minus['hamiltonian'])/(2*step)
        volume_rate = volume_difference/(2*step)
        result = dict(step=step, physical_energy_rate_watts=energy_rate,
            absolute_energy_work_error_watts=abs(energy_rate),
            physical_energy_rate_over_excitation_per_second=energy_rate/abs(base_energy) if base_energy else None,
            primal_hamiltonian_rate_watts=primal_rate,
            physical_minus_primal_energy_rate_watts=energy_rate-primal_rate,
            fluid_volume_rate_m3_per_second=volume_rate,
            absolute_fluid_volume_rate_error=abs(volume_rate),
            energy_within_diagnostic_tolerance=abs(energy_rate) <= energy_tolerance,
            volume_within_diagnostic_tolerance=abs(volume_rate) <= volume_tolerance,
            component_energy_rates_watts=dict(fluid_kinetic=fluid_difference/(2*step),
                rigid_kinetic=rigid_difference/(2*step), fluid_gravity=gravity_difference/(2*step),
                body_gravity=body_difference/(2*step)))
        results.append(result)
        print(json.dumps(dict(progress=True, stage='independent physical pair complete', **result)), flush=True)
        del plus, minus
    convergence = {}
    for key in ('absolute_energy_work_error_watts', 'absolute_fluid_volume_rate_error'):
        convergence[key+'_fine_over_coarse'] = results[1][key]/results[0][key] if results[0][key] else None
    return dict(qualified=False, case='saved actual nonlinear P4 canoe stage physical work audit',
        source_npz=str(stage_path), source_npz_sha256=sha256(stage_path.read_bytes()).hexdigest(),
        source_stage_metadata=source_metadata, audit_revision=_source_revision(),
        audit_script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        saved_rhs_recomputed=False,
        master_identity_assumption=master_identity,
        regenerated_elements_sha256=sha256(np.asarray(owner.mesh.elements, dtype='<i8').tobytes()).hexdigest(),
        regenerated_free_nodes_sha256=sha256(np.asarray(owner.free_nodes, dtype='<i8').tobytes()).hexdigest(),
        fresh_physical_base_maximum_absolute_errors=base_errors,
        base_physical_excitation_energy_joules=float(base_energy),
        base_primal_hamiltonian_joules=float(base['hamiltonian']),
        base_physical_minus_primal_energy_joules=float(base_energy-base['hamiltonian']),
        base_physical_energy_components_joules=dict(fluid_kinetic=fluid_kinetic,
            rigid_kinetic=rigid_kinetic, fluid_gravity_relative=gravity_relative, body_gravity_relative=body_gravity),
        independent_kinetic_quadrature_order=10, independent_boundary_quadrature_order=12,
        expected_undamped_energy_rate_watts=0., expected_closed_fluid_volume_rate=0.,
        diagnostic_tolerances=dict(energy_watts=energy_tolerance, excitation_fraction_per_1s=2e-6,
            energy_absolute_floor_watts=1e-7, fluid_volume_m3_per_second=volume_tolerance),
        differences=results, convergence=convergence,
        scope='One saved canonical RHS direction with fresh CURRENT perturbed solves. Physical kinetic '
            'energy uses direct quadrature velocities, gravity/volume use closed physical boundary '
            'quadrature, with integrands subtracted before aggregation. No projection, operator lag, '
            'RHS recomputation, trajectory, temporal/spatial/Capytaine or performance qualification.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', type=Path)
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-stage-physical-work.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        report = audit(args.stage)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
