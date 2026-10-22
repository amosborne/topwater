"""Physical-facet ALE mesh for the canonical canoe's moving waterline.

Fine contact heights remain independent. Contact on physical edges slides on
those edges; contact inside a facet follows a body-fixed tangent. The wet mesh
redistributes within each planar facet by a precomputed harmonic extension of
boundary displacement. Numerical subdivision vertices do not pin the chart.
"""
from dataclasses import dataclass

import numpy as np
from scipy.sparse import csr_matrix

from topwater._clipping import clip_flat_water
from topwater.frames import rotation_matrix
from topwater.hydrostatics import FloodingError
from tests.experiments.capytaine_mesh_cost import subdivide
from .canoe import join_wetted_hull


class ContactTopologyError(ValueError):
    """A physical facet transition or folded mesh requires a new contact chart."""


class RimCoverageError(ContactTopologyError):
    """The represented free surface does not cover a projected opening rim."""


def _rigid_jacobian(offset):
    out = np.zeros((len(offset), 3, 6))
    out[:, :, :3] = np.eye(3)
    out[:, :, 3:] = np.cross(np.eye(3)[None], offset[:, None]).transpose(0, 2, 1)
    return out


@dataclass(frozen=True)
class ContactGeometry:
    points: np.ndarray
    height_jacobians: tuple[csr_matrix, csr_matrix, csr_matrix]
    body_jacobian: np.ndarray

    def velocity(self, eta_rate, twist):
        """Nodal ALE velocity, excluding separately owned exterior wall motion."""
        rate = np.einsum('nij,j->ni', self.body_jacobian, twist)
        return rate+np.column_stack([matrix@eta_rate for matrix in self.height_jacobians])


@dataclass(frozen=True)
class CanoeContact:
    reference_points: np.ndarray
    faces: np.ndarray
    body_faces: np.ndarray
    free_faces: np.ndarray
    wall_faces: np.ndarray
    free_nodes: np.ndarray
    contact_nodes: np.ndarray
    wall_nodes: np.ndarray
    physical_vertices: np.ndarray
    physical_faces: np.ndarray
    body_source_faces: np.ndarray
    body_nodes: np.ndarray
    body_reference: np.ndarray
    coarse_nodes: np.ndarray
    contact_edges: np.ndarray
    contact_blend: np.ndarray
    contact_tangents: np.ndarray
    displacement_weights: csr_matrix
    cg: np.ndarray
    rim_vertices: np.ndarray

    @classmethod
    def make(cls, body, state, walls, surface_grid, *, subdivisions=2, constrained=False):
        """Clip physical facets at flat water, then subdivide the wet mesh.

        ``walls`` wind out of the fluid. ``eta`` is absolute world height indexed
        by sorted ``free_nodes``. Wall deformation belongs to the exterior mesh.
        Subdivision changes resolution, not the chart's original physical edges.
        """
        if isinstance(subdivisions, (bool, np.bool_)) or not isinstance(
                subdivisions, (int, np.integer)) or subdivisions < 0:
            raise ValueError('subdivisions must be a nonnegative integer')
        rotation = rotation_matrix(state.orientation)
        world = (body.mesh.vertices-body.cg)@rotation.T+state.position
        if np.any(world[:, 2] == 0):
            raise ContactTopologyError('initial waterline passes through a physical vertex')
        wet, source, segments, *_ = clip_flat_water(world, body.mesh.faces,
            body.mesh.physical_faces, 0., True)
        if not len(wet) or not len(segments):
            raise ContactTopologyError('initial canoe requires a closed waterline')
        coarse_line = np.unique(segments.reshape(-1, 3), axis=0)
        for _ in range(subdivisions):
            wet = subdivide(wet)
            source = np.tile(source, 4)
        fine_edges = np.stack((wet, np.roll(wet, -1, axis=1)), axis=2).reshape(-1, 2, 3)
        segments = fine_edges[np.all(fine_edges[:, :, 2] == 0, axis=1)]
        mesh, body_faces, free_faces = join_wetted_hull(wet, segments, walls, surface_grid,
                                                       constrained=constrained)
        lookup = {tuple(p): i for i, p in enumerate(mesh.vertices)}
        body_nodes = np.unique(mesh.faces[body_faces])
        free_nodes = np.unique(mesh.faces[free_faces])
        contact_nodes = np.intersect1d(body_nodes, free_nodes)
        coarse_nodes = np.array([lookup[tuple(p)] for p in coarse_line])
        physical = body.mesh.faces[body.mesh.physical_faces]
        edges = np.unique(np.sort(np.concatenate((physical[:, [0, 1]],
            physical[:, [1, 2]], physical[:, [2, 0]])), axis=1), axis=0)
        edge_lookup = {}
        for i, j in edges:
            if (world[i, 2] < 0) == (world[j, 2] < 0):
                continue
            a, b = (i, j) if world[i, 2] < 0 else (j, i)
            p = world[a]-world[a, 2]/(world[b, 2]-world[a, 2])*(world[b]-world[a])
            p[2] = 0.
            edge_lookup[tuple(p)] = (a, b)
        contact_edges = np.array([edge_lookup[tuple(p)] for p in coarse_line])
        blend = np.zeros((len(contact_nodes), len(coarse_nodes)))
        tangents = np.zeros((len(contact_nodes), 3))
        coarse_index = {node: i for i, node in enumerate(coarse_nodes)}
        contact_index = {node: i for i, node in enumerate(contact_nodes)}
        normals = np.cross(world[body.mesh.faces[:, 1]]-world[body.mesh.faces[:, 0]],
                           world[body.mesh.faces[:, 2]]-world[body.mesh.faces[:, 0]])
        normals /= np.linalg.norm(normals, axis=1)[:, None]
        for i, node in enumerate(contact_nodes):
            if node in coarse_index:
                blend[i, coarse_index[node]] = 1.
                continue
            parent = np.unique(source[np.any(mesh.faces[body_faces] == node, axis=1)])
            if len(parent) != 1:
                raise ValueError('fine contact must belong to one physical facet')
            face = body.mesh.faces[parent[0]]
            controls = [j for j, edge in enumerate(contact_edges) if np.all(np.isin(edge, face))]
            if len(controls) != 2:
                raise ValueError('a cut physical facet must have two edge contacts')
            a, b = coarse_line[controls]
            delta = b-a
            fraction = (mesh.vertices[node]-a)@delta/(delta@delta)
            blend[i, controls] = [1-fraction, fraction]
            n = normals[parent[0]]
            tangent = np.array([0., 0., 1.])-n[2]*n
            tangents[i] = tangent@rotation

        weights = np.zeros((len(mesh.vertices), len(contact_nodes)))
        weights[contact_nodes] = np.eye(len(contact_nodes))
        coarse_contact_indices = np.array([contact_index[node] for node in coarse_nodes])
        for parent in np.unique(source):
            faces = mesh.faces[body_faces[source == parent]]
            patch_nodes = np.unique(faces)
            patch_edges, counts = np.unique(np.sort(np.concatenate((faces[:, [0, 1]],
                faces[:, [1, 2]], faces[:, [2, 0]])), axis=1), axis=0, return_counts=True)
            boundary = np.unique(patch_edges[counts == 1])
            corners = body.mesh.faces[parent]
            triangle = world[corners]
            edge_vectors = triangle[1:]-triangle[:1]
            inverse_gram = np.linalg.inv(edge_vectors@edge_vectors.T)
            for node in np.setdiff1d(boundary, contact_nodes):
                uv = (mesh.vertices[node]-triangle[0])@edge_vectors.T@inverse_gram
                bary = np.r_[1-uv.sum(), uv]
                opposite = int(np.argmin(np.abs(bary)))
                edge = corners[[j for j in range(3) if j != opposite]]
                if (world[edge[0], 2] < 0) != (world[edge[1], 2] < 0):
                    a, b = edge if world[edge[0], 2] < 0 else edge[::-1]
                    control = np.flatnonzero(np.all(contact_edges == [a, b], axis=1))[0]
                    delta = coarse_line[control]-world[a]
                    fraction = (mesh.vertices[node]-world[a])@delta/(delta@delta)
                    weights[node, coarse_contact_indices[control]] = fraction
            interior = np.setdiff1d(patch_nodes, boundary)
            if len(interior):
                numbering = {node: i for i, node in enumerate(patch_nodes)}
                laplacian = np.zeros((len(patch_nodes), len(patch_nodes)))
                for a, b in patch_edges:
                    i, j = numbering[a], numbering[b]
                    laplacian[i, i] += 1
                    laplacian[j, j] += 1
                    laplacian[i, j] -= 1
                    laplacian[j, i] -= 1
                ii = np.array([numbering[node] for node in interior])
                bb = np.array([numbering[node] for node in boundary])
                weights[interior] = np.linalg.solve(laplacian[np.ix_(ii, ii)],
                    -laplacian[np.ix_(ii, bb)]@weights[boundary])
        wall_faces = np.arange(len(walls))
        body_reference = (mesh.vertices[body_nodes]-state.position)@rotation+body.cg
        return cls(mesh.vertices, mesh.faces, body_faces, free_faces, wall_faces,
            free_nodes, contact_nodes, np.unique(mesh.faces[wall_faces]),
            np.array(body.mesh.vertices), np.array(body.mesh.faces), source,
            body_nodes, body_reference, coarse_nodes, contact_edges, blend, tangents,
            csr_matrix(weights), np.array(body.cg),
            np.array(body.mesh.vertices[body.mesh.opening_rim]))

    def evaluate(self, position, orientation, eta):
        eta = np.asarray(eta)
        if eta.shape != (len(self.free_nodes),) or not np.all(np.isfinite(eta)):
            raise ValueError('eta must be finite and indexed by free_nodes')
        rotation = rotation_matrix(orientation)
        world = (self.physical_vertices-self.cg)@rotation.T+position
        coarse_eta = eta[np.searchsorted(self.free_nodes, self.coarse_nodes)]
        a, b = world[self.contact_edges[:, 0]], world[self.contact_edges[:, 1]]
        edge = b-a
        if np.any(edge[:, 2] <= 0):
            raise ContactTopologyError('physical contact edge no longer crosses upward')
        fraction = (coarse_eta-a[:, 2])/edge[:, 2]
        if np.any((fraction <= 0) | (fraction >= 1)):
            i = np.flatnonzero((fraction <= 0) | (fraction >= 1))[0]
            endpoint = self.contact_edges[i, 0 if fraction[i] <= 0 else 1]
            raise ContactTopologyError(f'waterline reached physical vertex {endpoint} '
                f'on contact edge {tuple(self.contact_edges[i])}; remeshing required')
        coarse = a+fraction[:, None]*edge
        coarse[:, 2] = coarse_eta
        slope = edge/edge[:, 2, None]
        coarse_body = _rigid_jacobian(coarse-position)
        coarse_body -= slope[:, :, None]*coarse_body[:, 2:3, :]
        coarse_height = np.zeros((len(coarse), 3, len(self.free_nodes)))
        for i, column in enumerate(np.searchsorted(self.free_nodes, self.coarse_nodes)):
            coarse_height[i, :, column] = slope[i]
        contact = self.contact_blend@coarse
        contact_body = np.einsum('ic,cjk->ijk', self.contact_blend, coarse_body)
        contact_height = np.einsum('ic,cjk->ijk', self.contact_blend, coarse_height)
        fine = ~np.isin(self.contact_nodes, self.coarse_nodes)
        tangent = self.contact_tangents[fine]@rotation.T
        if np.any(np.abs(tangent[:, 2]) <= np.finfo(float).eps):
            raise ContactTopologyError('facet contact tangent is horizontal')
        tangent /= tangent[:, 2, None]
        heights = eta[np.searchsorted(self.free_nodes, self.contact_nodes)]
        mismatch = heights[fine]-contact[fine, 2]
        contact[fine] += mismatch[:, None]*tangent
        derivative = _rigid_jacobian(tangent)
        derivative[:, :, :3] = 0.
        derivative -= tangent[:, :, None]*derivative[:, 2:3, :]
        contact_body[fine] -= tangent[:, :, None]*contact_body[fine, 2:3, :]
        contact_body[fine] += mismatch[:, None, None]*derivative
        contact_height[fine] -= tangent[:, :, None]*contact_height[fine, 2:3, :]
        for row, v in zip(np.flatnonzero(fine), tangent):
            contact_height[row, :, np.searchsorted(self.free_nodes, self.contact_nodes[row])] += v
        contact[:, 2] = heights
        contact_body[:, 2] = 0.
        contact_height[:, 2] = 0.
        contact_height[np.arange(len(contact)), 2,
            np.searchsorted(self.free_nodes, self.contact_nodes)] = 1.
        points = self.reference_points.copy()
        rigid = (self.body_reference-self.cg)@rotation.T+position
        points[self.body_nodes] = rigid
        rigid_contact = points[self.contact_nodes].copy()
        rigid_body = np.zeros((len(points), 3, 6))
        rigid_body[self.body_nodes] = _rigid_jacobian(rigid-position)
        points += self.displacement_weights@(contact-rigid_contact)
        correction = self.displacement_weights@(
            contact_body-rigid_body[self.contact_nodes]).reshape(len(contact), 18)
        body_jacobian = rigid_body+correction.reshape(len(points), 3, 6)
        height_jacobians = [self.displacement_weights@csr_matrix(contact_height[:, axis])
                            for axis in range(3)]
        ordinary = np.setdiff1d(self.free_nodes, self.contact_nodes)
        height_jacobians[2] += csr_matrix((np.ones(len(ordinary)),
            (ordinary, np.searchsorted(self.free_nodes, ordinary))),
            shape=(len(points), len(self.free_nodes)))
        points[self.free_nodes, 2] = eta
        self._check_patch(points, world)
        return ContactGeometry(points, tuple(height_jacobians), body_jacobian)

    def _check_patch(self, points, world):
        triangles = points[self.faces[self.body_faces]]
        parent = world[self.physical_faces[self.body_source_faces]]
        edges = parent[:, 1:]-parent[:, :1]
        normal = np.cross(edges[:, 0], edges[:, 1])
        area = np.einsum('fi,fi->f', normal,
            np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]))
        if np.any(area >= 0):  # Boundary winds into body, opposite physical hull.
            raise ContactTopologyError('wet physical-facet ALE patch folded')
        gram = np.einsum('fji,fki->fjk', edges, edges)
        uv = np.einsum('fci,fji->fcj', triangles-parent[:, :1], edges)
        uv = np.einsum('fcj,fjk->fck', uv, np.linalg.inv(gram))
        tolerance = 2048*np.finfo(float).eps
        if np.any(uv < -tolerance) or np.any(uv.sum(axis=2) > 1+tolerance):
            raise ContactTopologyError('waterline or wet mesh left its physical facet')
        distance = np.einsum('fci,fi->fc', triangles-parent[:, :1], normal)
        if np.any(np.abs(distance) > tolerance*np.linalg.norm(normal, axis=1)[:, None]):
            raise ContactTopologyError('wet mesh left its physical facet plane')
        p = points[self.faces[self.free_faces]]
        if np.any(np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0])[:, 2] <= 0):
            raise ContactTopologyError('free-surface projected triangle inverted')

    def rim_clearance(self, points, position, orientation):
        """Exact minimum clearance over rim segments on the current P1 graph.

        Clip each projected rim segment against existing surface triangles.
        Clearance is affine on each resulting interval, so its minimum is at
        an interval endpoint. No surface is extrapolated into uncovered areas.
        Coverage is checked before classifying immersion as physical flooding.
        """
        if not len(self.rim_vertices):
            return None
        rim = (self.rim_vertices-self.cg)@rotation_matrix(orientation).T+position
        triangles = np.asarray(points)[self.faces[self.free_faces]]
        lower = triangles[:, :, :2].min(axis=1)
        upper = triangles[:, :, :2].max(axis=1)
        minimum = np.inf
        tolerance = 512*np.finfo(float).eps  # Parameter-space roundoff only.
        for a, b in zip(rim, np.roll(rim, -1, axis=0)):
            candidates = np.all(upper >= np.minimum(a[:2], b[:2]), axis=1)&np.all(
                lower <= np.maximum(a[:2], b[:2]), axis=1)
            p = triangles[candidates]
            if not len(p):
                raise RimCoverageError('opening rim lies outside represented free surface')
            edges = p[:, 1:, :2]-p[:, :1, :2]
            inverse = np.linalg.inv(edges)
            at_a = np.einsum('fi,fij->fj', a[:2]-p[:, 0, :2], inverse)
            change = np.einsum('i,fij->fj', b[:2]-a[:2], inverse)
            at_a = np.column_stack((1-at_a.sum(axis=1), at_a))
            change = np.column_stack((-change.sum(axis=1), change))
            intervals = []
            for triangle, bary, delta in zip(p, at_a, change):
                lo, hi = 0., 1.
                for value, rate in zip(bary, delta):
                    if rate > 0:
                        lo = max(lo, -value/rate)
                    elif rate < 0:
                        hi = min(hi, -value/rate)
                    elif value < -tolerance:
                        hi = -1.
                if hi < lo:
                    continue
                intervals.append((lo, hi))
                for t in (lo, hi):
                    clearance = a[2]+t*(b[2]-a[2])-(bary+t*delta)@triangle[:, 2]
                    minimum = min(minimum, clearance)
            covered = 0.
            for lo, hi in sorted(intervals):
                if lo > covered+tolerance:
                    raise RimCoverageError('opening rim crosses an unrepresented surface region')
                covered = max(covered, hi)
            if covered < 1-tolerance:
                raise RimCoverageError('opening rim crosses an unrepresented surface region')
        if minimum <= 0:
            raise FloodingError(f'opening rim reaches the water surface (clearance {minimum:g} m)')
        return float(minimum)
