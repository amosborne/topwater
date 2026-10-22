"""Reproducible actual nonlinear P4 canoe seed and scalar-action falsifier.

Default initialization and action differences are the verified original stage
experiment; --calm selects a separate initially flat-surface physical case.
This is an experimental case driver, not a public solver API
or trajectory/Capytaine qualification.
"""
import argparse
from hashlib import sha256
from importlib.metadata import version
import json
import subprocess
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from topwater.frames import rotation_matrix
from .canonical_canoe import body_map
from .p4_whole_geometry import P4WholeCanoeGeometry
from .p4_canoe_field import prepare, boundary_moments
from .p4_variational_canoe import encode, evaluate
from .p4_geometry import integrals as cell_integrals
from .p4_boundary_action import integrals as face_integrals
from .p4_tetra import face_nodes


def progress(stage, **values):
    print(json.dumps(dict(stage=stage, **values)), flush=True)


def _array_hash(array, dtype):
    return sha256(np.asarray(array, dtype=dtype).tobytes()).hexdigest()


def verify_master_identity(geometry, stage_path):
    """Check new seed IDs exactly; explicitly retain legacy ordering assumptions."""
    names = ('reference_points', 'mesh_elements', 'free_nodes')
    with np.load(stage_path) as saved:
        present = [name in saved for name in names]
        if any(present) and not all(present):
            raise ValueError('stage master identity arrays must be present all-or-none')
        if all(present):
            actual = (geometry.mesh.points, geometry.mesh.elements, geometry.free_nodes)
            for name, array in zip(names, actual):
                if not np.array_equal(saved[name], array):
                    raise ValueError('saved/regenerated master identity mismatch: '+name)
            return ('Saved reference_points, mesh_elements and ordered free_nodes exactly '
                    'match the regenerated geometry; physical-state checks are additional.')
    return ('Legacy stage NPZ lacks explicit master identity arrays. Same source and '
            'deterministic mesher/master ordering are assumed; regenerated hashes and '
            'the independent fresh physical-state match are additional checks, not '
            'proof of identity with absent saved IDs.')


def initial_state(geometry):
    """Verified finite-pose/all-six-velocity seed, including high contact modes."""
    xy = geometry.mesh.points[geometry.free_nodes, :2]
    eta = 8e-5*np.cos(1.3*xy[:, 0])+3e-5*np.sin(2*xy[:, 1])
    high = geometry.chart.high_contact_nodes
    columns = geometry.chart_eta_columns[np.searchsorted(geometry.chart.free_nodes, high)]
    eta[columns] += 2e-6*np.sin(np.arange(len(columns)))
    q = np.r_[geometry.domain.rest.position+[.003, -.002, .005], [.008, -.006, .004], eta]
    psi = 2e-4*np.cos(.9*xy[:, 0]+.6*xy[:, 1])+4e-5*np.sin(1.3*xy[:, 0]-.7*xy[:, 1])
    psi[columns] += 1e-6*np.cos(np.arange(len(columns)))
    velocity = np.array([.01, .015, -.007, .02, -.03, .01])
    return q, psi, velocity


def calm_state(geometry):
    """Same finite pose/world twist, initially flat surface and zero trace psi.

    The interior potential still solves the moving-body normal flux; this is
    not initially motionless water. Nonflat water develops under the full RHS.
    The original high-contact-mode case remains the default derivative audit.
    """
    q, psi, velocity = initial_state(geometry)
    q, psi = q.copy(), psi.copy()
    q[6:] = 0.
    psi[:] = 0.
    return q, psi, velocity


def fixed_action(geometry, phi, qdot):
    """Physical scalar L(q; fixed nodal phi, canonical qdot), orders8/9.

    Only geometry is reevaluated: no perturbed fluid solve or force derivative.
    This is the same primal action used by the original verified stage driver.
    """
    nodes_on_faces = face_nodes(geometry.mesh, geometry.faces)
    _, reference_gravity = boundary_moments(geometry.mesh, geometry.faces, geometry.mesh.points)

    def scalar_action(configuration):
        orientation = (Rotation.from_rotvec(configuration[3:6])*Rotation.from_quat(
            geometry.domain.rest.orientation, scalar_first=True)).as_quat(scalar_first=True)
        state = geometry.evaluate(configuration[:3], orientation, configuration[6:])
        twist = body_map(configuration[3:6])@qdot[:6]
        motion = state.velocity(qdot[6:], twist)
        boundary = sum(face_integrals(state.points[nodes], phi[nodes], motion[nodes])['action']
                       for nodes in nodes_on_faces)
        kinetic = sum(cell_integrals(state.points[nodes], phi[nodes])['kinetic_energy']
                      for nodes in geometry.mesh.elements)
        _, fluid_gravity = boundary_moments(geometry.mesh, geometry.faces, state.points)
        gravity = (np.sum(fluid_gravity-reference_gravity)
                   +geometry.domain.body.mass*9.81*(configuration[2]-geometry.domain.rest.position[2]))
        omega_body = rotation_matrix(orientation).T@twist[3:]
        rigid = (.5*geometry.domain.body.mass*np.dot(twist[:3], twist[:3])
                 +.5*omega_body@geometry.domain.body.inertia_cg@omega_body)
        return float(rigid+1000*boundary-kinetic-gravity)

    return scalar_action


def run(*, calm=False):
    """One owner/current solve and the verified fixed-phi action differences."""
    start = perf_counter()
    progress('build actual whole-canoe geometry')
    geometry = P4WholeCanoeGeometry.make()
    progress('geometry ready', active_nodes=len(np.unique(geometry.mesh.elements)),
             surface_nodes=len(geometry.free_nodes), geometry_preflight=geometry.preflight)
    q, psi, velocity = calm_state(geometry) if calm else initial_state(geometry)
    eta = q[6:]
    progress('prepare fresh current curved fluid operator')
    data = prepare(geometry, q)
    progress('current field ready', fluid_volume=data['fluid_volume'],
             added_mass_diagonal=np.diag(data['added_mass']).tolist())
    canonical = encode(data, psi, velocity)
    progress('evaluate full nonlinear canonical stage')
    result = evaluate(geometry, data, canonical[len(q):])
    progress('canonical stage evaluated', hamiltonian=result['hamiltonian'],
             normal_relative_residual=result['full_weak_normal_velocity_relative_residual'],
             fluid_work_absolute_residual=result['fluid_work_identity_absolute_residual'])
    phi, qdot = result['potential'], result['qdot']
    scalar_action = fixed_action(geometry, phi, qdot)

    direction = np.r_[[.07, -.04, .03, .02, -.03, .015],
                      3e-4*np.cos(np.arange(len(eta))*.9)]
    expected = float(result['rate'][len(q):]@direction)
    differences = []
    for step in (2e-4, 1e-4):
        progress('independent fixed-phi/qdot scalar-action differences', step=step)
        actual = (scalar_action(q+step*direction)-scalar_action(q-step*direction))/(2*step)
        differences.append(dict(step=step, finite_difference=actual, expected=expected,
                                absolute_error=abs(actual-expected)))
        progress('action difference completed', **differences[-1])
    sources = ['p4_canoe_stage.py', 'p4_canoe_field.py', 'p4_variational_canoe.py',
        'p4_whole_geometry.py', 'p4_volume_map.py', 'p4_canoe_contact.py',
        'p4_geometry.py', 'p4_surface.py', 'p4_boundary_action.py', 'p4_tetra.py',
        'canonical_canoe.py']
    file_hashes = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                   for name in sources}
    mesh_hashes = dict(reference_points_sha256=_array_hash(geometry.mesh.points, '<f8'),
        elements_sha256=_array_hash(geometry.mesh.elements, '<i8'),
        free_nodes_sha256=_array_hash(geometry.free_nodes, '<i8'),
        active_nodes_sha256=_array_hash(np.unique(geometry.mesh.elements), '<i8'))
    report = dict(qualified=False, case='actual nonlinear P4 canonical-canoe stage',
        initial_case='flat-surface-release' if calm else 'nonflat-high-contact-audit',
        initial_case_description=('Original finite pose and all six world velocities; eta=psi=0. '
            'Compatible interior potential follows the moving-body flux, not globally still water. '
            'All nonlinear dynamics remain active; no deliberate high-contact-mode excitation.'
            if calm else 'Original verified nonflat eta/psi including deliberate independent '
            'high-contact modes, finite pose and all six world velocities.'),
        source_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'],
            cwd=Path(__file__).resolve().parents[3], text=True).strip(),
        source_sha256=file_hashes, master_hashes=mesh_hashes,
        package_versions={name: version(name) for name in ('meshpy', 'numpy', 'scipy')},
        master_identity='NPZ contains reference_points, ordered free_nodes, sorted active_nodes and mesh_elements; '
            'all use the shared global node IDs of the deterministic whole geometry. '
            'Reference coordinates are identified by the recorded float64 hash.',
        active_p4_nodes=len(np.unique(geometry.mesh.elements)), eta_modes=len(eta),
        body_world_velocity=velocity.tolist(), pose=q[:6].tolist(),
        hamiltonian=result['hamiltonian'], fluid_volume=result['fluid_volume'],
        full_weak_normal_velocity_residual=result['full_weak_normal_velocity_residual'],
        full_weak_normal_velocity_relative_residual=result['full_weak_normal_velocity_relative_residual'],
        fluid_work_identity_absolute_residual=result['fluid_work_identity_absolute_residual'],
        action_directional_differences=differences, diagnostic_seconds=perf_counter()-start,
        scope='One current nonlinear stage with all six velocities and finite pose. '
              'The initial_case and its description identify the surface excitation. Scalar action '
              'differences hold phi and coordinate qdot fixed; no perturbed fluid solve. '
              'Not independent physical Hamiltonian work along the full RHS, temporal/spatial '
              'accuracy, P4 rim/failure coverage, trajectory, Capytaine or throughput qualification.')
    arrays = dict(q=q, canonical=canonical, canonical_rate=result['rate'],
        physical=result['physical_state'], free_nodes=geometry.free_nodes,
        active_nodes=np.unique(geometry.mesh.elements), mesh_elements=geometry.mesh.elements,
        reference_points=geometry.mesh.points,
        action_direction=direction, source_revision=np.asarray(report['source_revision']),
        initial_case=np.asarray(report['initial_case']),
        source_sha256_json=np.asarray(json.dumps(file_hashes, sort_keys=True)),
        master_hashes_json=np.asarray(json.dumps(mesh_hashes, sort_keys=True)))
    return report, arrays


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calm', action='store_true',
                        help='same finite pose/all-six twist with initial eta=psi=0')
    parser.add_argument('--output', type=Path,
                        help='default p4-calm-stage.json with --calm, p4-current-stage.json otherwise')
    args = parser.parse_args()
    output = args.output or Path('artifacts/p4-calm-stage.json' if args.calm
                                else 'artifacts/p4-current-stage.json')
    with threadpool_limits(1):
        report, arrays = run(calm=args.calm)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output.with_suffix('.npz'), **arrays)
    output.write_text(json.dumps(report, indent=2)+'\n')
    progress('completed', **report)
