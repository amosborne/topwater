"""Current curved P4 field for the actual canoe geometry candidate.

This is a validation reference. Interior elimination is applied through solves;
the volume-to-surface extension and dense surface Schur matrix are never formed.
No nonlinear stage or open-water qualification follows from this helper alone.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.spatial.transform import Rotation

from topwater.frames import rotation_matrix
from .canonical_canoe import body_map
from .p4_tetra import face_nodes
from .p4_surface import quadrature, integrals as surface_integrals
from .p4_geometry import stiffness


def boundary_moments(mesh, faces, points, *, order=9):
    """Physical closed-boundary volume/gravity contributions, face by face."""
    values, derivative, weights = quadrature(order)
    coordinates = points[face_nodes(mesh, faces)]
    tangent = np.einsum('qaj,fai->fqji', derivative, coordinates)
    normal = np.cross(tangent[:, :, 0], tangent[:, :, 1])
    height = np.einsum('qa,fa->fq', values, coordinates[:, :, 2])
    return (np.einsum('q,fq,fq->f', weights, height, normal[:, :, 2]),
            4905*np.einsum('q,fq,fq->f', weights, height**2, normal[:, :, 2]))


def body_load(mesh, faces, points, center, *, order=9):
    """Curved physical hull flux tested against every P4 boundary basis."""
    values, derivative, weights = quadrature(order)
    nodes = face_nodes(mesh, faces)
    coordinates = points[nodes]
    tangent = np.einsum('qaj,fai->fqji', derivative, coordinates)
    normal = np.cross(tangent[:, :, 0], tangent[:, :, 1])
    samples = np.einsum('qa,fai->fqi', values, coordinates)
    modes = np.concatenate((normal, np.cross(samples-center, normal)), axis=2)
    local = np.einsum('q,qa,fqk->fak', weights, values, modes)
    result = np.zeros((len(points), 6))
    np.add.at(result, nodes.ravel(), local.reshape(-1, 6))
    return result


def surface_maps(geometry, state, *, order=9):
    """Sparse Ceta/mass and six-column Cv, using actual geometry Jacobians."""
    nodes = face_nodes(geometry.mesh, geometry.free_faces)
    free_lookup = np.full(len(state.points), -1, dtype=int)
    free_lookup[geometry.free_nodes] = np.arange(len(geometry.free_nodes))
    boundary_lookup = np.full(len(state.points), -1, dtype=int)
    boundary_lookup[geometry.volume_map.boundary] = np.arange(len(geometry.volume_map.boundary))
    eta_rows, eta_columns, eta_values = [], [], []
    mass_rows, mass_columns, mass_values = [], [], []
    transport = np.zeros((len(geometry.free_nodes), 6))
    for face in nodes:
        rows = free_lookup[face]
        boundary_rows = boundary_lookup[face]
        if np.any(rows < 0) or np.any(boundary_rows < 0):
            raise ValueError('physical surface trace is missing its shared geometry nodes')
        height = [matrix[boundary_rows] for matrix in state.boundary_height_jacobians]
        columns = np.unique(np.concatenate([matrix.indices for matrix in height]))
        local = surface_integrals(state.points[face], tuple(matrix[:, columns] for matrix in height),
            state.boundary_body_jacobian[boundary_rows], order=order)
        eta_rows.append(np.repeat(rows, len(columns)))
        eta_columns.append(np.tile(columns, 15))
        eta_values.append(local['eta_map'].ravel())
        mass_rows.append(np.repeat(rows, 15))
        mass_columns.append(np.tile(rows, 15))
        mass_values.append(local['mass'].ravel())
        np.add.at(transport, rows, local['body_transport'])
    shape = (len(geometry.free_nodes),)*2
    eta_map = coo_matrix((np.concatenate(eta_values),
        (np.concatenate(eta_rows), np.concatenate(eta_columns))), shape=shape).tocsc()
    mass = coo_matrix((np.concatenate(mass_values),
        (np.concatenate(mass_rows), np.concatenate(mass_columns))), shape=shape).tocsr()
    return eta_map, transport, mass


def prepare(geometry, q, *, volume_order=8, boundary_order=9):
    """Assemble the current physical system once for canonical decoding.

    q=(CG position, local rotation vector from rest, all surface heights).
    Density/gravity follow the existing canonical-canoe reference.
    """
    q = np.array(q, copy=True)
    if q.shape != (6+len(geometry.free_nodes),) or not np.all(np.isfinite(q)):
        raise ValueError('finite body configuration and all P4 surface heights required')
    if np.linalg.norm(q[3:6]) >= np.pi:
        raise ValueError('local rotation chart boundary reached; recharting required')
    domain = geometry.domain
    orientation = (Rotation.from_rotvec(q[3:6])*Rotation.from_quat(
        domain.rest.orientation, scalar_first=True)).as_quat(scalar_first=True)
    state = geometry.evaluate(q[:3], orientation, q[6:])
    # Surface orientation is checked at its own quadrature, independently of
    # sampled volume positivity and the local contact chart's guards.
    physical_moments = geometry.boundary_moments(state, order=boundary_order)
    energy = stiffness(geometry.mesh, state.points, order=volume_order)
    load = body_load(geometry.mesh, geometry.body_faces, state.points, q[:3], order=boundary_order)
    active = np.unique(geometry.mesh.elements)
    free = geometry.free_nodes
    inside = np.setdiff1d(active, free)
    factor = splu(energy[inside][:, inside].tocsc())
    response = factor.solve(load[inside])
    added = 1000*load[inside].T@response
    eta_map, body_transport, mass = surface_maps(geometry, state, order=boundary_order)
    eta_factor = splu(eta_map)
    coupling = load[free]-energy[free][:, inside]@response+body_transport
    rotation = rotation_matrix(orientation)
    rigid = np.zeros((6, 6))
    rigid[:3, :3] = domain.body.mass*np.eye(3)
    rigid[3:, 3:] = rotation@domain.body.inertia_cg@rotation.T
    volume, gravity = boundary_moments(geometry.mesh, geometry.faces, state.points, order=boundary_order)
    _, reference_gravity = boundary_moments(geometry.mesh, geometry.faces,
                                           geometry.mesh.points, order=boundary_order)
    gravity_relative = float(np.sum(gravity-reference_gravity)
                             +domain.body.mass*9.81*(q[2]-domain.rest.position[2]))
    return dict(q=q, orientation=orientation, state=state, energy=energy, load=load,
        free=free, inside=inside, factor=factor, response=response,
        added_mass=added, rigid_mass=rigid, effective_mass=rigid+added,
        eta_map=eta_map, eta_factor=eta_factor, body_transport=body_transport,
        coupling=coupling, mass=mass, gravity=gravity_relative,
        fluid_volume=float(volume.sum()), body_map=body_map(q[3:6]),
        volume_order=volume_order, boundary_order=boundary_order,
        physical_boundary_moments=physical_moments)


def solve(data, psi, velocity):
    """Current Dirichlet/body-flux field; every body surface row is retained."""
    psi, velocity = np.asarray(psi), np.asarray(velocity)
    if (psi.shape != (len(data['free']),) or velocity.shape != (6,)
            or not np.all(np.isfinite(psi)) or not np.all(np.isfinite(velocity))):
        raise ValueError('finite P4 surface potential and world body twist required')
    inside, free = data['inside'], data['free']
    surface_response = data['factor'].solve(-data['energy'][inside][:, free]@psi)
    potential = np.zeros(data['energy'].shape[0])
    potential[inside] = surface_response+data['response']@velocity
    potential[free] = psi
    flux = (data['energy']@potential-data['load']@velocity)[free]
    eta_rate = data['eta_factor'].solve(flux-data['body_transport']@velocity)
    surface_only = potential.copy()
    surface_only[inside] = surface_response
    kinetic = float(500*surface_only@(data['energy']@surface_only))
    hamiltonian = float(.5*velocity@data['effective_mass']@velocity+kinetic+data['gravity'])
    return dict(potential=potential, flux=flux, eta_rate=eta_rate,
                hamiltonian=hamiltonian, surface_kinetic_energy=kinetic)
