"""Whole current canoe P4 geometry chart; no fluid operator or dynamics.

Only physical boundary coordinates are prescribed. Artificial pocket walls
are interior to the merged mesh and move through the fixed reference ALE map.
Sampled positive Jacobians are necessary, not a validity certificate.

The reference edge-spring graph alone failed the unchanged wet2 finite-pose
control: even a corner-only cell inverted (determinant ratio -0.34677).
The reference P1 FEM corner extension kept that same control positive; the
hierarchical map therefore uses the graph only for genuine P4 boundary modes.
This records a rejected coordinate map, not fluid-accuracy qualification.
"""
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, diags

from .moving_canoe import make_domain
from .pocket import boundary
from .p4_whole_canoe import make_volume, face_keys
from .p4_canoe_contact import P4CanoeContact
from .p4_volume_map import make_hierarchical_map
from .p4_geometry import quadrature as volume_quadrature
from .p4_surface import quadrature as surface_quadrature
from . import p4_tetra


@dataclass(frozen=True)
class GeometryState:
    owner: object
    points: np.ndarray
    boundary_height_jacobians: tuple
    boundary_body_jacobian: np.ndarray

    def velocity(self, eta_rate, world_twist):
        eta_rate, world_twist = np.asarray(eta_rate), np.asarray(world_twist)
        if (eta_rate.shape != (len(self.owner.free_nodes),) or world_twist.shape != (6,)
                or not np.all(np.isfinite(eta_rate)) or not np.all(np.isfinite(world_twist))):
            raise ValueError('finite global eta rate and world six-twist required')
        value = np.einsum('nij,j->ni', self.boundary_body_jacobian, world_twist)
        value += np.column_stack([matrix@eta_rate for matrix in self.boundary_height_jacobians])
        return self.owner.volume_map.lift(value)

    def pullback(self, nodal_force):
        """Return world-body and global-height forces: exact transpose chart."""
        force = self.owner.volume_map.pullback(nodal_force)
        body = np.einsum('nij,ni->j', self.boundary_body_jacobian, force)
        height = sum(matrix.T@force[:, axis]
                     for axis, matrix in enumerate(self.boundary_height_jacobians))
        return body, height


@dataclass(frozen=True)
class P4WholeCanoeGeometry:
    domain: object
    mesh: object
    faces: np.ndarray
    free_faces: np.ndarray
    body_faces: np.ndarray
    volume_map: object
    chart: object
    chart_to_global: np.ndarray
    chart_eta_columns: np.ndarray
    override_chart_nodes: np.ndarray
    override_boundary_rows: np.ndarray
    free_nodes: np.ndarray
    free_boundary_rows: np.ndarray
    preflight: dict

    @classmethod
    def make(cls, domain=None):
        domain = make_domain(12, False, 2) if domain is None else domain
        original, old_faces, _, body_faces, preflight = make_volume(domain)
        # The static builder accepts either orientation. Nonlinear calculus
        # requires every positive orientation BEFORE constructing its P4 IDs.
        points = original.points[:original.vertex_count].copy()
        tets = original.tets.copy()
        determinant = np.linalg.det(points[tets[:, 1:]]-points[tets[:, :1]])
        negative = determinant < 0
        tets[negative, 0], tets[negative, 1] = tets[negative, 1].copy(), tets[negative, 0].copy()
        mesh = p4_tetra.make_mesh(points, tets)
        faces = boundary(points, tets)
        if face_keys(faces) != face_keys(old_faces):
            raise RuntimeError('positive orientation changed the physical boundary')
        free_faces = faces[np.all(points[faces, 2] == 0, axis=1)]
        chart = P4CanoeContact.make(domain.contact)
        chart_to_global = np.full(len(chart.reference_points), -1, dtype=int)
        for face_id, face in enumerate(domain.contact.faces):
            for a, alpha in enumerate(p4_tetra.indices(3)):
                node = chart.faces[face_id, a]
                target = mesh.node_by_key[p4_tetra.key(domain.local_mapping[face], alpha)]
                if chart_to_global[node] not in (-1, target):
                    raise RuntimeError('inconsistent shared chart topology')
                chart_to_global[node] = target
        if np.any(chart_to_global < 0) or not np.allclose(
                mesh.points[chart_to_global], chart.reference_points, atol=2e-15, rtol=0):
            raise RuntimeError('chart/reference topology or coordinates differ')
        free = np.unique(p4_tetra.face_nodes(mesh, free_faces))
        local_free = chart_to_global[chart.free_nodes]
        eta_columns = np.searchsorted(free, local_free)
        if np.any(eta_columns >= len(free)) or not np.array_equal(free[eta_columns], local_free):
            raise RuntimeError('local free heights do not belong to global surface masters')
        prescribed = np.unique(p4_tetra.face_nodes(mesh, faces))
        active = np.unique(mesh.elements)
        # Count the actual graph input before its single reference factor. This
        # bounds input storage only; exposed factor storage is checked after it.
        alpha = p4_tetra.indices(4)
        adjacent = np.sum(abs(alpha[:, None]-alpha[None, :]), axis=2) == 2
        local_edges = np.column_stack(np.where(np.triu(adjacent, 1)))
        edges = np.unique(np.sort(mesh.elements[:, local_edges].reshape(-1, 2), axis=1), axis=0)
        graph_entries = 2*len(edges)+len(active)
        graph_bytes = 12*graph_entries+4*(len(mesh.points)+1)
        if len(active) > 100000 or len(edges) > 1000000 or graph_bytes > 64000000:
            raise RuntimeError('geometry graph exceeded its input preflight caps')
        preflight = dict(preflight, active_geometry_nodes=len(active),
            physical_boundary_nodes=len(prescribed), graph_edges=len(edges),
            graph_input_entries_upper_bound=graph_entries,
            graph_input_bytes_upper_bound=graph_bytes, positively_reoriented_tets=int(negative.sum()))
        del edges, original
        volume_map = make_hierarchical_map(mesh)
        # Report both actual factors; these are exposed arrays, not peak-memory
        # bounds. Input CSR storage does not bound factorization workspace.
        total_bytes = 0
        for name, mapping in (('p4_residual_graph', volume_map.graph),
                              ('p1_fem_corners', volume_map.corners)):
            lower, upper = mapping.factor.L, mapping.factor.U
            factor_bytes = sum(value.data.nbytes+value.indices.nbytes+value.indptr.nbytes
                               for value in (lower, upper))
            preflight.update({name+'_lower_nonzeros': lower.nnz,
                              name+'_upper_nonzeros': upper.nnz,
                              name+'_exposed_factor_array_bytes': int(factor_bytes)})
            total_bytes += factor_bytes
            del lower, upper
        if total_bytes > 2000000000:
            raise RuntimeError('geometry factors exceeded exposed-array storage cap')
        preflight['total_exposed_factor_array_bytes'] = int(total_bytes)
        row = np.full(len(mesh.points), -1, dtype=int)
        row[volume_map.boundary] = np.arange(len(volume_map.boundary))
        override = np.unique(chart.faces[np.r_[domain.contact.body_faces, domain.contact.free_faces]])
        override_rows = row[chart_to_global[override]]
        if np.any(override_rows < 0) or len(np.unique(override_rows)) != len(override_rows):
            raise RuntimeError('physical chart rows must have unique boundary ownership')
        return cls(domain, mesh, faces, free_faces, body_faces, volume_map, chart,
                   chart_to_global, eta_columns, override, override_rows, free, row[free], preflight)

    def evaluate(self, position, orientation, eta):
        eta = np.asarray(eta)
        if eta.shape != (len(self.free_nodes),) or not np.all(np.isfinite(eta)):
            raise ValueError('finite eta at every global P4 free-surface master required')
        chart = self.chart.evaluate(position, orientation, eta[self.chart_eta_columns])
        boundary_points = self.mesh.points[self.volume_map.boundary].copy()
        boundary_points[self.free_boundary_rows, 2] = eta
        boundary_points[self.override_boundary_rows] = chart.points[self.override_chart_nodes]
        shape = (len(boundary_points), len(eta))
        default_z = coo_matrix((np.ones(len(eta)), (self.free_boundary_rows, np.arange(len(eta)))),
                               shape=shape).tocsr()
        keep = np.ones(len(boundary_points)); keep[self.override_boundary_rows] = 0
        height = []
        for axis in range(3):
            local = chart.height_jacobians[axis][self.override_chart_nodes].tocoo()
            matrix = coo_matrix((local.data, (self.override_boundary_rows[local.row],
                                             self.chart_eta_columns[local.col])), shape=shape).tocsr()
            if axis == 2:
                matrix += diags(keep)@default_z
            height.append(matrix)
        body = np.zeros((len(boundary_points), 3, 6))
        body[self.override_boundary_rows] = chart.body_jacobian[self.override_chart_nodes]
        displacement = boundary_points-self.mesh.points[self.volume_map.boundary]
        points = self.mesh.points+self.volume_map.lift(displacement)
        return GeometryState(self, points, tuple(height), body)

    def boundary_moments(self, state, *, order=9):
        """Closed physical volume/gravity and independent normal-work gradients."""
        values, derivative, weights = surface_quadrature(order)
        nodes = p4_tetra.face_nodes(self.mesh, self.faces)
        coordinates = state.points[nodes]
        tangent = np.einsum('qaj,fai->fqji', derivative, coordinates)
        normal = np.cross(tangent[:, :, 0], tangent[:, :, 1])
        reference = self.mesh.points[self.faces]
        reference_normal = np.cross(reference[:, 1]-reference[:, 0], reference[:, 2]-reference[:, 0])
        body_keys = face_keys(self.body_faces)
        body = np.array([tuple(sorted(face)) in body_keys for face in self.faces])
        # Body orientation is checked by the actual chart against CURRENT
        # physical facet normals, not against normals before its rigid rotation.
        if (not np.all(np.isfinite(normal)) or np.any(np.linalg.norm(normal, axis=2) == 0)
                or np.any(np.einsum('fqi,fi->fq', normal[~body], reference_normal[~body]) <= 0)):
            raise ValueError('nonpositive sampled physical boundary patch orientation')
        z = np.einsum('qa,fa->fq', values, coordinates[:, :, 2])
        volume_gradient = np.zeros_like(state.points)
        gravity_gradient = np.zeros_like(state.points)
        local_volume = np.einsum('q,qa,fqi->fai', weights, values, normal)
        local_gravity = 9810*np.einsum('q,qa,fq,fqi->fai', weights, values, z, normal)
        np.add.at(volume_gradient, nodes, local_volume)
        np.add.at(gravity_gradient, nodes, local_gravity)
        return dict(volume=float(np.einsum('q,fq,fq->', weights, z, normal[:, :, 2])),
            gravity=float(4905*np.einsum('q,fq,fq->', weights, z*z, normal[:, :, 2])),
            volume_gradient=volume_gradient, gravity_gradient=gravity_gradient)

    def volume_moments(self, state, *, order=8):
        """Physical cell moments with explicit sampled positive-Jacobian guard.

        This computes geometry only; no inverse, stiffness or fluid solve.
        Evaluation alone does not certify the interior map. The fluid assembler
        must retain its own current sampled guards at its integration points.
        """
        values, derivative, weights = volume_quadrature(order)
        volume = gravity = 0.; minimum = np.inf
        for start in range(0, len(self.mesh.elements), 32):
            coordinates = state.points[self.mesh.elements[start:start+32]]
            jacobian = np.einsum('qaj,tai->tqji', derivative, coordinates)
            determinant = np.linalg.det(jacobian)
            if not np.all(np.isfinite(determinant)) or np.any(determinant <= 0):
                raise ValueError('nonpositive sampled whole-canoe P4 volume Jacobian')
            z = np.einsum('qa,ta->tq', values, coordinates[:, :, 2])
            volume += np.einsum('q,tq->', weights, determinant)
            gravity += 9810*np.einsum('q,tq,tq->', weights, determinant, z)
            minimum = min(minimum, determinant.min())
        return dict(volume=float(volume), gravity=float(gravity),
                    minimum_sampled_jacobian=float(minimum))
