"""Trim exported surface geometry against a posed convex hull, without physics.

This is a bounded convex-halfspace operation, not a general mesh Boolean or a
fluid-state remesher. Input triangles retain their plane and winding. The
exported contact is shared by water triangles and agrees geometrically with
the hull; the rigid hull's stored connectivity is not changed.
"""

import itertools

import numpy as np
from scipy.spatial import cKDTree

from ._validation import array, indices
from .inputs import HullMesh


def _mesh_arrays(vertices_m, faces):
    vertices = array(vertices_m, (None, 3), 'surface vertices_m')
    raw = np.asarray(faces)
    if (raw.ndim != 2 or raw.shape[1:] != (3,) or raw.dtype.kind not in 'iu'
            or not len(raw) or len(vertices) < 3
            or np.any(raw < 0) or np.any(raw >= len(vertices))):
        raise ValueError('surface faces must be nonempty integer triangles with valid indices')
    faces = np.array(raw, dtype=np.int64, copy=True)
    if len(np.unique(np.sort(faces, axis=1), axis=0)) != len(faces):
        raise ValueError('surface has duplicate faces')
    if len(np.unique(faces)) != len(vertices):
        raise ValueError('surface has unused vertices')
    triangles = vertices[faces]
    a, b = triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]
    cross = np.cross(a, b)
    area = np.linalg.norm(cross, axis=1)
    if np.any(area <= 64*np.finfo(float).eps*np.linalg.norm(a, axis=1)*np.linalg.norm(b, axis=1)):
        raise ValueError('surface has degenerate triangles')
    return vertices, faces


def validate_surface_topology(vertices_m, faces):
    """Check a nonempty oriented 2-manifold mesh, with boundaries permitted.

    Edge incidence, shared orientation and vertex links are checked. Disjoint
    components, multiple boundaries and non-graph surfaces are permitted. This
    does not detect geometric self intersections or certify physical contact.
    """
    vertices, faces = _mesh_arrays(vertices_m, faces)
    edges, incident = {}, [[] for _ in vertices]
    for i, triangle in enumerate(faces):
        for v in triangle:
            incident[v].append(i)
        for a, b in zip(triangle, np.roll(triangle, -1)):
            edges.setdefault(tuple(sorted((int(a), int(b)))), []).append((int(a), int(b)))
    for owners in edges.values():
        if len(owners) > 2:
            raise ValueError('surface has a non-manifold edge')
        if len(owners) == 2 and owners[0] != owners[1][::-1]:
            raise ValueError('surface has inconsistent shared-edge winding')
    for vertex, owners in enumerate(incident):
        link = {}
        for i in owners:
            a, b = [int(v) for v in faces[i] if v != vertex]
            link.setdefault(a, set()).add(b); link.setdefault(b, set()).add(a)
        degrees = [len(neighbors) for neighbors in link.values()]
        if any(d not in (1, 2) for d in degrees) or degrees.count(1) not in (0, 2):
            raise ValueError('surface has a non-manifold vertex link')
        seen, todo = set(), [next(iter(link))]
        while todo:
            v = todo.pop()
            if v not in seen:
                seen.add(v); todo.extend(link[v]-seen)
        if len(seen) != len(link):
            raise ValueError('surface has a pinched vertex with disconnected links')


def _split_polygon(polygon, normal, offset, eps, vertex):
    """Return inside/outside polygon IDs against one outward convex plane."""
    inside, outside = [], []
    has_positive = False
    previous = polygon[-1]
    dp = vertex.points[previous] @ normal+offset
    if abs(dp) <= eps:
        dp = 0.
    for current in polygon:
        dc = vertex.points[current] @ normal+offset
        if abs(dc) <= eps:
            dc = 0.
        has_positive |= dc > 0
        if dp*dc < 0:
            # Always interpolate in index order for common original/cut edges.
            a, b = sorted((previous, current))
            pa, pb = vertex.points[a], vertex.points[b]
            da, db = pa @ normal+offset, pb @ normal+offset
            cut = vertex.add(pa+(pb-pa)*(da/(da-db)))
            inside.append(cut); outside.append(cut)
        if dc <= 0:
            inside.append(current)
        if dc >= 0:
            outside.append(current)
        previous, dp = current, dc
    return _clean(inside), _clean(outside) if has_positive else []


def _clean(polygon):
    clean = []
    for i in polygon:
        if not clean or clean[-1] != i:
            clean.append(i)
    if len(clean) > 1 and clean[0] == clean[-1]:
        clean.pop()
    return clean


def _area_vector(polygon, points):
    if len(polygon) < 3:
        return np.zeros(3)
    q = np.asarray([points[i] for i in polygon])
    q = q-q[0]
    return np.cross(q, np.roll(q, -1, axis=0)).sum(axis=0)/2


def _valid_triangle(points):
    a, b = points[None, 1]-points[None, 0], points[None, 2]-points[None, 0]
    cross = np.cross(a, b)
    threshold = 64*np.finfo(float).eps*np.linalg.norm(a, axis=1)*np.linalg.norm(b, axis=1)
    return bool(np.linalg.norm(cross, axis=1)[0] > threshold[0])


def _validate_reconstructed_partition(polygon, triangles, world_points, normal):
    """Reject a guessed quad fan unless it covers its oriented parent once."""
    parent = world_points[np.asarray(polygon)]
    relative = parent-parent[0]
    parent_vector = np.cross(relative, np.roll(relative, -1, axis=0)).sum(axis=0)
    child = world_points[np.asarray(triangles)]
    cross = np.cross(child[:, 1]-child[:, 0], child[:, 2]-child[:, 0])
    parent_area = float(np.linalg.norm(parent_vector))
    child_area = float(np.linalg.norm(cross, axis=1).sum())
    roundoff = 256*np.finfo(float).eps*max(parent_area, child_area)
    if (parent_vector @ normal <= 0 or np.any(cross @ normal <= 0)
            or abs(parent_area-child_area) > roundoff
            or np.linalg.norm(parent_vector-cross.sum(axis=0)) > roundoff):
        raise ValueError('reconstructed export triangulation does not preserve oriented parent coverage')


def _export_ears(polygon, points, normal, world_points):
    """Retain constrained vertices if projected ear roundoff leaves a sliver.

    The shared triangulator is unchanged. Its usual output is retained exactly.
    A triangle failing the ordinary FINAL-WORLD degeneracy predicate can only
    be repaired by flipping an internal diagonal with its neighbor: all polygon
    vertices and directed boundary edges stay fixed, both new triangles must
    pass that same predicate, and oriented geometry/area must be preserved.
    No triangles or contact constraints are dropped and no vertices are moved.
    """
    from .fluid_boundary import _ears

    reconstructed = False
    try:
        triangles = [list(face) for face in _ears(polygon, points, normal)]
    except ValueError as error:
        # Under rigid transforms the same four-vertex contact degeneracy can
        # fall just below the shared routine's projected final-turn tolerance.
        # Reconstruct its two-diagonal quad partition for the same WORLD tests
        # below. Larger/other triangulation failures remain fatal.
        if len(polygon) != 4 or str(error) != 'degenerate final wet-facet triangle':
            raise
        triangles = [[polygon[0], polygon[1], polygon[2]],
                     [polygon[2], polygon[3], polygon[0]]]
        reconstructed = True
    boundary = {tuple(sorted((int(a), int(b))))
                for a, b in zip(polygon, np.roll(polygon, -1))}
    while True:
        bad = [i for i, face in enumerate(triangles) if not _valid_triangle(world_points[face])]
        if not bad:
            if reconstructed:
                _validate_reconstructed_partition(polygon, triangles, world_points, normal)
            return triangles
        repaired = False
        for index in bad:
            face = triangles[index]
            for start, end in zip(face, np.roll(face, -1)):
                edge = tuple(sorted((int(start), int(end))))
                if edge in boundary:
                    continue
                neighbors = [j for j, other in enumerate(triangles)
                             if j != index and start in other and end in other]
                if len(neighbors) != 1:
                    continue
                neighbor = neighbors[0]
                other = triangles[neighbor]
                if (end, start) not in list(zip(other, np.roll(other, -1))):
                    continue
                opposite = next(v for v in face if v not in edge)
                opposite_other = next(v for v in other if v not in edge)
                if opposite == opposite_other:
                    continue
                replacements = [[opposite, start, opposite_other], [opposite, opposite_other, end]]
                if not all(_valid_triangle(world_points[t]) for t in replacements):
                    continue
                old_points = world_points[np.asarray([face, other])]
                new_points = world_points[np.asarray(replacements)]
                old_cross = np.cross(old_points[:, 1]-old_points[:, 0], old_points[:, 2]-old_points[:, 0])
                new_cross = np.cross(new_points[:, 1]-new_points[:, 0], new_points[:, 2]-new_points[:, 0])
                if np.any(new_cross @ normal <= 0):
                    continue
                old_area = float(np.linalg.norm(old_cross, axis=1).sum())
                new_area = float(np.linalg.norm(new_cross, axis=1).sum())
                roundoff = 256*np.finfo(float).eps*max(old_area, new_area)
                if (abs(old_area-new_area) > roundoff
                        or np.linalg.norm(old_cross.sum(axis=0)-new_cross.sum(axis=0)) > roundoff):
                    continue
                # Replacing a bad pair by two valid faces strictly reduces the
                # count of bad triangles, so this bounded fallback terminates.
                triangles[index], triangles[neighbor] = replacements
                repaired = True
                break
            if repaired:
                break
        if not repaired:
            raise ValueError('export polygon cannot retain all contact constraints with nondegenerate triangles')


class _Vertices:
    def __init__(self, points, eps):
        self.points = list(points.copy())
        self.original_tree = cKDTree(points)
        self.eps = eps
        self.cells = {}
        self.cuts = set()

    def add(self, point):
        near = self.original_tree.query_ball_point(point, self.eps)
        if len(near) > 1:
            raise ValueError('surface vertices are indistinguishable at contact tolerance; no implicit welding')
        if near:
            self.cuts.add(near[0])
            return near[0]
        key = tuple(np.floor(point/self.eps).astype(np.int64))
        found = []
        for delta in itertools.product((-1, 0, 1), repeat=3):
            cell = tuple(a+b for a, b in zip(key, delta))
            found.extend(i for i in self.cells.get(cell, ())
                         if np.linalg.norm(self.points[i]-point) <= self.eps)
        if len(found) > 1:
            raise ValueError('ambiguous cut-vertex contact at export tolerance')
        if found:
            self.cuts.add(found[0])
            return found[0]
        index = len(self.points)
        self.points.append(point)
        self.cells.setdefault(key, []).append(index)
        self.cuts.add(index)
        return index


def _opening_clearance(vertices, faces, hull, rim, eps):
    """Check complete rim edges against the represented piecewise-linear field."""
    if not len(rim):
        return
    triangles = vertices[faces]
    xy = triangles[:, :, :2]
    a, b = xy[:, 1]-xy[:, 0], xy[:, 2]-xy[:, 0]
    determinant = a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]
    lower, upper = xy.min(axis=1), xy.max(axis=1)
    for point, end in zip(hull[rim], hull[np.roll(rim, -1)]):
        candidates = np.flatnonzero(np.all(np.maximum(point[:2], end[:2]) >= lower-eps, axis=1)
                                    & np.all(np.minimum(point[:2], end[:2]) <= upper+eps, axis=1)
                                    & (abs(determinant) > eps**2))
        d = point[:2]-xy[candidates, 0]
        de = end[:2]-xy[candidates, 0]
        aa, bb = a[candidates], b[candidates]
        s = (d[:, 0]*bb[:, 1]-d[:, 1]*bb[:, 0])/determinant[candidates]
        t = (aa[:, 0]*d[:, 1]-aa[:, 1]*d[:, 0])/determinant[candidates]
        se = (de[:, 0]*bb[:, 1]-de[:, 1]*bb[:, 0])/determinant[candidates]
        te = (aa[:, 0]*de[:, 1]-aa[:, 1]*de[:, 0])/determinant[candidates]
        intervals = []
        for row, i in enumerate(candidates):
            start_weights = np.array([1-s[row]-t[row], s[row], t[row]])
            changes = np.array([s[row]+t[row]-se[row]-te[row], se[row]-s[row], te[row]-t[row]])
            lo, hi = 0., 1.
            for weight, change in zip(start_weights, changes):
                if abs(change) <= eps:
                    if weight < -eps:
                        lo, hi = 1., 0.; break
                elif change > 0:
                    lo = max(lo, (-eps-weight)/change)
                else:
                    hi = min(hi, (-eps-weight)/change)
            if hi < lo:
                continue
            z0 = start_weights @ triangles[i, :, 2]
            dz = changes @ triangles[i, :, 2]
            for parameter in (lo, hi):
                if z0+parameter*dz >= point[2]+parameter*(end[2]-point[2])-eps:
                    raise ValueError('surface reaches physical opening rim; flooding is unsupported')
            intervals.append((lo, hi, z0, dz))
        intervals.sort()
        covered = 0.
        for j, (lo, hi, z0, dz) in enumerate(intervals):
            if lo > covered+8*eps:
                raise ValueError('total surface does not provide complete opening-rim coverage')
            covered = max(covered, hi)
            for old_lo, old_hi, old_z, old_dz in intervals[:j]:
                overlap_lo, overlap_hi = max(lo, old_lo), min(hi, old_hi)
                if overlap_hi >= overlap_lo and max(abs(z0-old_z+x*(dz-old_dz))
                                                     for x in (overlap_lo, overlap_hi)) > 8*eps:
                    raise ValueError('multiple surface elevations at opening rim; local compatibility is ambiguous')
        if covered < 1.-8*eps:
            raise ValueError('total surface does not provide complete opening-rim coverage')


def _cap_contact(vertices, faces, hull, cap_faces, eps):
    from .boat_geometry import _triangle_contacts

    water = vertices[faces]
    lower, upper = water.min(axis=1), water.max(axis=1)
    for triangle in hull[cap_faces]:
        normal = np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0])
        normal /= np.linalg.norm(normal)
        candidates = np.flatnonzero(np.all(upper >= triangle.min(axis=0)-eps, axis=1)
                                    & np.all(lower <= triangle.max(axis=0)+eps, axis=1))
        for i in candidates:
            other = water[i]
            other_normal = np.cross(other[1]-other[0], other[2]-other[0])
            other_normal /= np.linalg.norm(other_normal)
            if _triangle_contacts(triangle, other, normal, other_normal, eps):
                raise ValueError('surface reaches nonphysical opening closure; flooding is unsupported')


def export_surface(vertices_m, faces, hull_vertices_m, hull_faces, *, opening_rim=(),
                   artificial_faces=(), tolerance_m=None):
    """Return total-field triangles outside a posed, closed convex hull.

    All inputs are WORLD metres. Hull faces wind outward from the body and
    include any internal nonphysical closure. opening_rim refers to that hull's
    vertices and artificial_faces identifies its internal cap triangles. Supply
    both together for an open boat, or neither for a sealed hull. The input
    total surface must cover the complete rim before trimming; local
    multiple elevations or rim immersion are rejected. This local compatibility
    check does not impose a graph or single-component constraint on the generic
    topology validator.

    Point/edge tangencies at physical hull panels remove no area and are
    retained. Contact with an open rim or nonphysical cap, including exact
    tangency, is rejected as unsupported flooding onset. Positive-area
    coplanar hull contact is ambiguous and rejected. New intersection vertices
    share indices within the water mesh, using an explicit tolerance defaulting
    to 1e-10 of hull span; ambiguous tolerance collisions fail without repairing
    physical geometry. No hull or fluid-state arrays are mutated.
    """
    vertices, faces = _mesh_arrays(vertices_m, faces)
    validate_surface_topology(vertices, faces)
    hull = HullMesh(hull_vertices_m, hull_faces)
    origin = hull.vertices.mean(axis=0)
    scale = float(np.max(np.ptp(hull.vertices, axis=0)))
    eps = 1e-10 if tolerance_m is None else float(tolerance_m)/scale
    if not np.isfinite(eps) or eps <= 0:
        raise ValueError('surface export tolerance_m must be positive and finite')
    hp = (hull.vertices-origin)/scale
    vp = (vertices-origin)/scale
    rim = indices(opening_rim, len(hp), 'opening_rim')
    artificial = indices(artificial_faces, len(hull.faces), 'artificial_faces')
    if bool(len(rim)) != bool(len(artificial)) or (len(rim) and len(rim) < 3):
        raise ValueError('open hull export needs opening_rim and artificial_faces together')
    if len(rim):
        HullMesh(hull.vertices, hull.faces, opening_rim=rim, artificial_faces=artificial)
    _opening_clearance(vp, faces, hp, rim, eps)
    _cap_contact(vp, faces, hp, hull.faces[artificial], eps)
    ht = hp[hull.faces]
    normals = np.cross(ht[:, 1]-ht[:, 0], ht[:, 2]-ht[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    planes = []
    for triangle, normal in zip(ht, normals):
        offset = -normal @ triangle[0]
        if np.max(hp @ normal+offset) > eps:
            raise ValueError('surface export requires an outward convex full hull')
        plane = np.r_[normal, offset]
        if not any(np.max(abs(plane-other)) <= eps for other in planes):
            planes.append(plane)
    planes = np.asarray(planes)
    triangles = vp[faces]
    lower, upper = triangles.min(axis=1), triangles.max(axis=1)
    nearby = (np.all(upper >= hp.min(axis=0)-eps, axis=1)
              & np.all(lower <= hp.max(axis=0)+eps, axis=1))
    uncut_faces = faces[~nearby]
    output = []
    registry = _Vertices(vp, eps)
    fragments = []
    for face in faces[nearby]:
        normal = np.cross(vp[face[1]]-vp[face[0]], vp[face[2]]-vp[face[0]])
        normal /= np.linalg.norm(normal)
        remaining = list(face)
        emitted = []
        coplanar = False
        for plane in planes:
            distances = vp[face] @ plane[:3]+plane[3]
            if np.max(abs(distances)) <= eps:
                coplanar = True
            remaining, outside = _split_polygon(remaining, plane[:3], plane[3], eps, registry)
            if len(outside) >= 3 and np.linalg.norm(_area_vector(outside, registry.points)) > eps**2:
                emitted.append(outside)
            if len(remaining) < 3 or np.linalg.norm(_area_vector(remaining, registry.points)) <= eps**2:
                remaining = []
                break
        if coplanar and len(remaining) >= 3 and np.linalg.norm(_area_vector(remaining, registry.points)) > eps**2:
            raise ValueError('coplanar hull/surface area contact has no unique exported waterline')
        fragments.extend((polygon, normal) for polygon in emitted)
    # Later clipping planes may subdivide an earlier fragment edge. A cut can
    # also lie on an original edge shared with a broadphase-rejected triangle:
    # early plane cuts extend outside the hull. Conform both kinds of neighbor
    # without sending those exterior triangles through the clipping operation.
    # Clipping is complete: no vertices are added during conformance or ears.
    # Reuse this exact array rather than rebuilding the full registry for every
    # polygon; keep the shared triangulator's indices and projection unchanged.
    registry_points = np.asarray(registry.points)
    registry_world = registry_points*scale+origin
    registry_world[:len(vertices)] = vertices
    if registry.cuts:
        cut_ids = np.array(sorted(registry.cuts), dtype=np.int64)
        cut_points = registry_points[cut_ids]
        cut_tree = cKDTree(cut_points)
    def conform(polygon):
        conforming = []
        for a, b in zip(polygon, np.roll(polygon, -1)):
            pa, pb = registry.points[a], registry.points[b]
            edge = pb-pa
            length = np.linalg.norm(edge)
            if length <= eps:
                raise ValueError('export fragment edge is unresolved at contact tolerance')
            conforming.append(a)
            if registry.cuts:
                candidates = cut_tree.query_ball_point((pa+pb)/2, length/2+eps)
                hits = []
                for k in candidates:
                    i = int(cut_ids[k])
                    t = (registry.points[i]-pa) @ edge/(length*length)
                    if eps/length < t < 1-eps/length and np.linalg.norm(registry.points[i]-pa-t*edge) <= eps:
                        hits.append((t, i))
                conforming.extend(i for _, i in sorted(hits))
        return _clean(conforming)

    if registry.cuts:
        uncut_triangles = vp[uncut_faces]
        affected = (np.all(uncut_triangles.max(axis=1) >= cut_points.min(axis=0)-eps, axis=1)
                    & np.all(uncut_triangles.min(axis=1) <= cut_points.max(axis=0)+eps, axis=1))
        output.extend(uncut_faces[~affected].tolist())
        for face in uncut_faces[affected]:
            polygon = conform(face)
            if len(polygon) == 3:
                output.append(list(face))
            else:
                triangle = vp[face]
                normal = np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0])
                normal /= np.linalg.norm(normal)
                output.extend(_export_ears(polygon, registry_points, normal, registry_world))
    else:
        output.extend(uncut_faces.tolist())
    for polygon, normal in fragments:
        output.extend(_export_ears(conform(polygon), registry_points, normal, registry_world))
    if not output:
        raise ValueError('surface export has no water triangles outside the hull')
    used, inverse = np.unique(np.asarray(output, dtype=np.int64), return_inverse=True)
    output_vertices = registry_world[used]
    output_faces = inverse.reshape(-1, 3)
    validate_surface_topology(output_vertices, output_faces)
    return output_vertices, output_faces
