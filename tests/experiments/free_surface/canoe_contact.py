"""Body-fitted ALE waterline for the canonical canoe, with fixed connectivity.

This is the geometry of a moving case, not a remesher. Contact nodes slide on
their original material edges. A material vertex crossing the waterline ends
the admissible topology; the caller must reject the step or remesh explicitly.
Interface/wall deformation belongs to the coupled fluid mesh and is not inferred
here. All twists and the body Jacobian use world translation/rotation axes.
"""
from dataclasses import dataclass

import numpy as np

from topwater.frames import rotation_matrix
from topwater.hydrostatics import FloodingError
from .canoe import join_wetted_hull
from .material_layout import material_mesh


class ContactTopologyError(ValueError):
    """The waterline has reached a material vertex or inverted a surface face."""


class RimCoverageError(ContactTopologyError):
    """The represented free surface does not cover a projected opening rim."""


def _rigid_jacobian(offset):
    out = np.zeros((len(offset), 3, 6))
    out[:, :, :3] = np.eye(3)
    # Each column is e_axis cross offset.
    out[:, :, 3:] = np.cross(np.eye(3)[None], offset[:, None]).transpose(0, 2, 1)
    return out


@dataclass(frozen=True)
class ContactGeometry:
    points: np.ndarray
    eta_jacobian: np.ndarray
    body_jacobian: np.ndarray

    def velocity(self, free_nodes, eta_rate, twist):
        """Nodal ALE velocity at this geometry, excluding exterior wall motion."""
        rate = np.einsum('nij,j->ni', self.body_jacobian, twist)
        rate[free_nodes] += self.eta_jacobian[free_nodes]*np.asarray(eta_rate)[:, None]
        return rate


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
    material_vertices: np.ndarray
    material_nodes: np.ndarray
    material_ids: np.ndarray
    contact_edges: np.ndarray
    cg: np.ndarray
    rim_vertices: np.ndarray

    @classmethod
    def make(cls, body, state, walls, surface_grid, *, transverse=False, subdivisions=2):
        """Initialize at flat z=0 using a rigid canonical hull inside the pocket.

        ``walls`` are fluid-outward triangles, ``surface_grid`` supplies the top
        of that box. The returned face winding is fluid-outward everywhere.
        ``eta`` in evaluate is indexed by ``free_nodes``, including contact.
        """
        vertices, faces = material_mesh(body, transverse, subdivisions=subdivisions)
        world = (vertices-body.cg)@rotation_matrix(state.orientation).T+state.position
        if np.any(world[:, 2] == 0):
            raise ContactTopologyError('initial waterline passes through a material vertex')
        coordinates = {}
        records = {}
        triangles = []
        segments = []

        def vertex(i):
            key = ('v', int(i))
            coordinates[key] = world[i]
            records[key] = (int(i),)
            return key

        def crossing(i, j):
            # Same wet-to-dry arithmetic as physical clipping; incident faces
            # must share the identical coordinate, not tolerance-welded copies.
            wet, dry = (i, j) if world[i, 2] < 0 else (j, i)
            key = ('e', int(wet), int(dry))
            t = -world[wet, 2]/(world[dry, 2]-world[wet, 2])
            p = world[wet]+t*(world[dry]-world[wet])
            p[2] = 0.
            coordinates[key] = p
            records[key] = (int(wet), int(dry))
            return key

        for face in faces:
            polygon = []
            previous = face[-1]
            for current in face:
                inside = world[current, 2] < 0
                if inside != (world[previous, 2] < 0):
                    polygon.append(crossing(previous, current))
                if inside:
                    polygon.append(vertex(current))
                previous = current
            for j in range(1, len(polygon)-1):
                triangles.append([coordinates[k] for k in (polygon[0], polygon[j], polygon[j+1])])
            for a, b in zip(polygon, polygon[1:]+polygon[:1]):
                if a[0] == b[0] == 'e':
                    segments.append([coordinates[a], coordinates[b]])
        if not triangles or not segments:
            raise ContactTopologyError('initial canoe requires a closed waterline')
        mesh, body_faces, free_faces = join_wetted_hull(
            np.asarray(triangles), np.asarray(segments), walls, surface_grid)
        lookup = {tuple(p): i for i, p in enumerate(mesh.vertices)}
        material_nodes = []
        material_ids = []
        contact_nodes = []
        contact_edges = []
        for key, point in coordinates.items():
            if key[0] == 'v':
                material_nodes.append(lookup[tuple(point)])
                material_ids.append(records[key][0])
            else:
                contact_nodes.append(lookup[tuple(point)])
                contact_edges.append(records[key])
        wall_faces = np.arange(len(walls))
        return cls(mesh.vertices, mesh.faces, body_faces, free_faces, wall_faces,
            np.unique(mesh.faces[free_faces]), np.asarray(contact_nodes),
            np.unique(mesh.faces[wall_faces]), vertices, np.asarray(material_nodes),
            np.asarray(material_ids), np.asarray(contact_edges), np.array(body.cg),
            np.array(body.mesh.vertices[body.mesh.opening_rim]))

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

    def evaluate(self, position, orientation, eta):
        eta = np.asarray(eta)
        if eta.shape != (len(self.free_nodes),) or not np.all(np.isfinite(eta)):
            raise ValueError('eta must be finite and indexed by free_nodes')
        world = (self.material_vertices-self.cg)@rotation_matrix(orientation).T+position
        points = self.reference_points.copy()
        points[self.material_nodes] = world[self.material_ids]
        points[self.free_nodes, 2] = eta
        eta_jacobian = np.zeros_like(points)
        eta_jacobian[self.free_nodes, 2] = 1.
        body_jacobian = np.zeros((len(points), 3, 6))
        body_jacobian[self.material_nodes] = _rigid_jacobian(
            points[self.material_nodes]-position)
        a, b = world[self.contact_edges[:, 0]], world[self.contact_edges[:, 1]]
        d = b-a
        heights = points[self.contact_nodes, 2]
        if np.any(d[:, 2] <= 0):
            raise ContactTopologyError('contact material edge no longer crosses upward')
        fraction = (heights-a[:, 2])/d[:, 2]
        if np.any((fraction <= 0) | (fraction >= 1)):
            i = np.flatnonzero((fraction <= 0) | (fraction >= 1))[0]
            endpoint = self.contact_edges[i, 0 if fraction[i] <= 0 else 1]
            raise ContactTopologyError(f'waterline reached material vertex {endpoint} '
                f'on contact edge {tuple(self.contact_edges[i])}; remeshing required')
        points[self.contact_nodes] = a+fraction[:, None]*d
        points[self.contact_nodes, 2] = heights
        slope = d/d[:, 2, None]
        eta_jacobian[self.contact_nodes] = slope
        rigid = _rigid_jacobian(points[self.contact_nodes]-position)
        body_jacobian[self.contact_nodes] = rigid-slope[:, :, None]*rigid[:, 2:3, :]
        p = points[self.faces[self.free_faces]]
        if np.any(np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0])[:, 2] <= 0):
            raise ContactTopologyError('free-surface projected triangle inverted')
        return ContactGeometry(points, eta_jacobian, body_jacobian)
