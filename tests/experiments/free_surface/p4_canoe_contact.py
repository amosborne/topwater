"""Actual canonical canoe's enriched P4 boundary chart; no volume solver.

Reuse the established physical-facet P1 chart at its old nodes. Additional
waterline heights slide in their physical facet, producing curved waterlines
without curving the physical hull. Other new wet nodes remain barycentric
lifts. Sampled patch guards delimit this small-motion chart, not certification
or general remeshing/flooding support.
"""
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags

from topwater.frames import rotation_matrix
from .canoe_contact import CanoeContact, ContactGeometry, ContactTopologyError, _rigid_jacobian
from .p4_tetra import indices, key, basis
from .p4_surface import quadrature


@dataclass(frozen=True)
class P4CanoeContact:
    linear: CanoeContact
    reference_points: np.ndarray
    faces: np.ndarray
    lift: csr_matrix
    free_nodes: np.ndarray
    contact_nodes: np.ndarray
    high_contact_nodes: np.ndarray
    contact_tangents: np.ndarray
    old_eta_indices: np.ndarray

    @classmethod
    def make(cls, contact):
        """Enrich fixed boundary connectivity; retain every original vertex ID."""
        count = len(contact.reference_points)
        lookup = {((i, 4),): i for i in range(count)}
        rows, columns, weights = list(range(count)), list(range(count)), [1.]*count
        faces = np.empty((len(contact.faces), 15), dtype=int)
        for face_id, face in enumerate(contact.faces):
            for j, alpha in enumerate(indices(3)):
                identity = key(face, alpha)
                if identity not in lookup:
                    node = len(lookup)
                    lookup[identity] = node
                    for vertex, numerator in identity:
                        rows.append(node); columns.append(vertex); weights.append(numerator/4)
                faces[face_id, j] = lookup[identity]
        lift = coo_matrix((weights, (rows, columns)), shape=(len(lookup), count)).tocsr()
        free = np.unique(faces[contact.free_faces])
        body = np.unique(faces[contact.body_faces])
        contact_nodes = np.intersect1d(free, body)
        high = contact_nodes[contact_nodes >= count]
        tangents = []
        for node in high:
            parents = np.unique(contact.body_source_faces[
                np.any(faces[contact.body_faces] == node, axis=1)])
            support = lift.indices[lift.indptr[node]:lift.indptr[node+1]]
            if len(parents) != 1 or len(support) != 2 or not np.all(np.isin(support, contact.contact_nodes)):
                raise ContactTopologyError('new P4 contact requires one physical facet and an old waterline edge')
            triangle = contact.physical_vertices[contact.physical_faces[parents[0]]]
            normal = np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0])
            normal /= np.linalg.norm(normal)
            tangents.append(np.array([0., 0., 1.])-normal[2]*normal)
        return cls(contact, lift@contact.reference_points, faces, lift, free,
                   contact_nodes, high, np.asarray(tangents).reshape(-1, 3),
                   np.searchsorted(free, contact.free_nodes))

    def evaluate(self, position, orientation, eta):
        eta = np.asarray(eta)
        if eta.shape != (len(self.free_nodes),) or not np.all(np.isfinite(eta)):
            raise ValueError('finite eta indexed by enriched free_nodes required')
        old = self.linear.evaluate(position, orientation, eta[self.old_eta_indices])
        selector = csr_matrix((np.ones(len(self.old_eta_indices)),
            (np.arange(len(self.old_eta_indices)), self.old_eta_indices)),
            shape=(len(self.old_eta_indices), len(self.free_nodes)))
        points = self.lift@old.points
        height = [self.lift@matrix@selector for matrix in old.height_jacobians]
        body = (self.lift@old.body_jacobian.reshape(len(old.points), 18)).reshape(-1, 3, 6)
        high = self.high_contact_nodes
        tangent = self.contact_tangents@rotation_matrix(orientation).T
        if np.any(np.abs(tangent[:, 2]) <= np.finfo(float).eps):
            raise ContactTopologyError('P4 facet contact tangent is horizontal')
        tangent /= tangent[:, 2, None]
        own = np.searchsorted(self.free_nodes, high)
        mismatch = eta[own]-points[high, 2]
        points[high] += mismatch[:, None]*tangent
        derivative = _rigid_jacobian(tangent)
        derivative[:, :, :3] = 0.
        derivative -= tangent[:, :, None]*derivative[:, 2:3]
        body[high] -= tangent[:, :, None]*body[high, 2:3]
        body[high] += mismatch[:, None, None]*derivative
        extraction = csr_matrix((np.ones(len(high)), (np.arange(len(high)), high)),
                                shape=(len(high), len(points)))
        own_height = csr_matrix((np.ones(len(high)), (np.arange(len(high)), own)),
                               shape=(len(high), len(self.free_nodes)))
        correction = own_height-extraction@height[2]
        height = [matrix+extraction.T@diags(tangent[:, axis])@correction
                  for axis, matrix in enumerate(height)]
        # All enriched FS heights are independent absolute z, including the
        # ordinary high nodes that do not belong to the hull.
        points[self.free_nodes, 2] = eta
        body[self.free_nodes, 2] = 0.
        keep = np.ones(len(points)); keep[self.free_nodes] = 0.
        height[2] = diags(keep)@height[2]+csr_matrix((np.ones(len(eta)),
            (self.free_nodes, np.arange(len(eta)))), shape=(len(points), len(eta)))
        self._check(points, position, orientation)
        return ContactGeometry(points, tuple(matrix.tocsr() for matrix in height), body)

    def _check(self, points, position, orientation):
        """Current P4 sampled geometry, including interpolation nodes; no certificate."""
        values, derivative, _ = quadrature(9)
        node_values, node_derivative = basis(indices(3)/4)
        values = np.vstack((values, node_values))
        derivative = np.concatenate((derivative, node_derivative[:, :, 1:]-node_derivative[:, :, :1]))
        body_faces = self.faces[self.linear.body_faces]
        samples = np.einsum('qa,fai->fqi', values, points[body_faces])
        tangent = np.einsum('qaj,fai->fqji', derivative, points[body_faces])
        area = np.cross(tangent[:, :, 0], tangent[:, :, 1])
        world = (self.linear.physical_vertices-self.linear.cg)@rotation_matrix(orientation).T+position
        parent = world[self.linear.physical_faces[self.linear.body_source_faces]]
        edges = parent[:, 1:]-parent[:, :1]
        normal = np.cross(edges[:, 0], edges[:, 1])
        tolerance = 4096*np.finfo(float).eps
        distance = np.einsum('fqi,fi->fq', samples-parent[:, :1], normal)
        if np.any(abs(distance) > tolerance*np.linalg.norm(normal, axis=1)[:, None]):
            raise ContactTopologyError('sampled P4 wet patch left its physical plane')
        if np.any(np.einsum('fqi,fi->fq', area, normal) >= 0):
            raise ContactTopologyError('sampled P4 wet patch folded')
        gram = np.einsum('fji,fki->fjk', edges, edges)
        uv = np.einsum('fqi,fji->fqj', samples-parent[:, :1], edges)
        uv = np.einsum('fqj,fjk->fqk', uv, np.linalg.inv(gram))
        if np.any(uv < -tolerance) or np.any(uv.sum(axis=2) > 1+tolerance):
            raise ContactTopologyError('sampled P4 wet patch left its physical facet')
        tangent = np.einsum('qaj,fai->fqji', derivative,
                            points[self.faces[self.linear.free_faces]])
        if np.any(np.cross(tangent[:, :, 0], tangent[:, :, 1])[:, :, 2] <= 0):
            raise ContactTopologyError('sampled P4 free-surface graph folded')
