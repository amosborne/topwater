"""Exact-arithmetic TWO-parent boat clipping reference, not a hot-path mesher.

This bounded code proof exercises a real chine and a transom coincident with
the shared vertical parent face. Plane/vertex identities are exact incidence
sets, never proximity welds. No volume mesh, fluid solve or evolved state.
"""
from fractions import Fraction as F
from itertools import combinations
from collections import Counter, defaultdict
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.spatial import ConvexHull

from topwater import HullMesh
from tests.reference_boat import reference_body
from .boat_chart_rebuild import affine, normal_work
from .p4_fluid_ordering import fresh


# Authored base-facet groups, retaining the original triangle mapping.
GROUPS = ((0, 1), (2,), (3, 4), (5,), (9, 10), (11,), (12, 13),
          (14,), (15, 16, 17), (6, 7, 8))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


class PlanePoints:
    """Tiny reference registry: all planes are frozen before intersections."""
    def __init__(self, *, max_points=None, max_intersections=None, capture_constraints=False,
                 cache_signs=False, integer_predicates=False):
        self.planes, self.labels, self.cache = [], [], {}
        self.points, self.incidence, self.vertex_by_owner = [], [], {}
        self.registered_intersections = {}
        self.relevant = [] if capture_constraints else None
        self.max_points, self.max_intersections = max_points, max_intersections
        self.sign_cache = {} if cache_signs else None
        self.sign_requests, self.sign_cache_hits = 0, 0
        self.integer_predicates = bool(integer_predicates)
        self.integer_geometry = None

    def _integer_owner(self):
        if self.integer_geometry is None:
            from .boat_integer_intersections import IntegerIntersections
            self.integer_geometry = IntegerIntersections(self.planes)
        self.integer_geometry.check_planes(self.planes)
        # One immutable owner after validation; subsequent checks are identity.
        self.planes = self.integer_geometry.frozen_planes
        return self.integer_geometry

    def plane(self, points, label):
        if self.cache or self.integer_geometry is not None:
            raise RuntimeError('freeze every plane before computing intersections')
        a, b, c = [tuple(x if isinstance(x, F) else F(float(x)) for x in p) for p in points]
        n = cross(sub(b, a), sub(c, a))
        first = next(x for x in n if x)
        coefficients = tuple(x/first for x in (*n, dot(n, a)))
        if coefficients not in self.planes:
            self.planes.append(coefficients)
            self.labels.append([])
        i = self.planes.index(coefficients)
        self.labels[i].append(label)
        return i, 1 if first > 0 else -1

    def value(self, i, point):
        plane = self.planes[i]
        return dot(plane[:3], point)-plane[3]

    def exact_signs(self, points):
        """Dense exact cold witness signs; no extra predicate cache."""
        points = list(points)
        signs = (self._integer_owner().signs(points) if self.integer_predicates else
                 [[int(np.sign(self.value(index, point))) for index in range(len(self.planes))]
                  for point in points])
        return np.asarray(signs, dtype=np.int8).reshape(len(points), len(self.planes))

    def _sign(self, triple, point, plane):
        """Exact cached predicate; masks retain signs, never large Fractions."""
        if self.integer_geometry is not None and self.planes is not self.integer_geometry.frozen_planes:
            self._integer_owner()
        self.sign_requests += 1
        bit = 1 << plane
        known, zero, positive = ((0, 0, 0) if self.sign_cache is None else
                                self.sign_cache.get(triple, (0, 0, 0)))
        if known & bit:
            self.sign_cache_hits += 1
            return 0 if zero & bit else (1 if positive & bit else -1)
        if self.integer_predicates:
            owner = self._integer_owner()
            sign = (owner.sign(triple, plane) if point is self.cache.get(triple)
                    else owner.sign_point(point, plane))
        else:
            value = self.value(plane, point)
            sign = (value > 0)-(value < 0)
        if self.sign_cache is not None:
            self.sign_cache[triple] = known | bit, zero | (bit if sign == 0 else 0), positive | (bit if sign > 0 else 0)
        return sign

    def intersection(self, triple):
        if self.integer_geometry is not None and self.planes is not self.integer_geometry.frozen_planes:
            self._integer_owner()
        triple = tuple(sorted(triple))
        if triple in self.cache:
            return self.cache[triple]
        if self.max_intersections is not None and len(self.cache) >= self.max_intersections:
            raise RuntimeError('reference intersection cap reached')
        if self.integer_predicates:
            point = self._integer_owner().intersection(triple)
            self.cache[triple] = point
            return point
        a, b, c = (self.planes[i] for i in triple)
        determinant = dot(a[:3], cross(b[:3], c[:3]))
        if not determinant:
            self.cache[triple] = None
            return None
        vectors = cross(b[:3], c[:3]), cross(c[:3], a[:3]), cross(a[:3], b[:3])
        point = tuple(sum(plane[3]*v[j] for plane, v in zip((a, b, c), vectors))/determinant
                      for j in range(3))
        self.cache[triple] = point
        return point

    def register(self, triple, point):
        triple = tuple(sorted(triple))
        if triple in self.registered_intersections:
            return self.registered_intersections[triple]
        # Only accepted polytope vertices need global incidence. Outside-plane
        # intersections cannot be required by the subsequent boundary splitting.
        owner = tuple(i for i in range(len(self.planes)) if not self._sign(triple, point, i))
        if owner not in self.vertex_by_owner:
            if self.max_points is not None and len(self.points) >= self.max_points:
                raise RuntimeError('reference point cap reached')
            self.vertex_by_owner[owner] = len(self.points)
            self.points.append(point)
            self.incidence.append(frozenset(owner))
            if self.relevant is not None:
                self.relevant.append(set())
        index = self.vertex_by_owner[owner]
        assert self.points[index] == point
        self.registered_intersections[triple] = index
        return index

    def polytope(self, constraints):
        vertices = set()
        for triple in combinations(sorted({i for i, _ in constraints}), 3):
            point = self.intersection(triple)
            if point is not None and all(sign*self._sign(triple, point, i) <= 0 for i, sign in constraints):
                vertices.add(self.register(triple, point))
        if self.relevant is not None:
            for vertex in vertices:
                self.relevant[vertex].update(i for i, _ in constraints)
        return vertices

    def polygon(self, vertices, plane, sign):
        ids = sorted(v for v in vertices if plane in self.incidence[v])
        if len(ids) < 3:
            return []
        # Drop a lower-dimensional intersection by EXACT incidence, not epsilon.
        a = self.points[ids[0]]
        if not any(any(cross(sub(self.points[b], a), sub(self.points[c], a)))
                   for b, c in combinations(ids[1:], 2)):
            return []
        xyz = np.array([[float(x) for x in self.points[i]] for i in ids])
        center = xyz.mean(axis=0)
        n = sign*np.array(self.planes[plane][:3], dtype=float)
        n /= np.linalg.norm(n)
        u = xyz[0]-center
        u /= np.linalg.norm(u)
        v = np.cross(n, u)
        order = np.argsort(np.arctan2((xyz-center)@v, (xyz-center)@u))
        return [ids[i] for i in order]

    def split_edges(self, polygon, candidates):
        out = []
        for a, b in zip(polygon, polygon[1:]+polygon[:1]):
            common = self.incidence[a]&self.incidence[b]
            direction = sub(self.points[b], self.points[a])
            length2 = dot(direction, direction)
            cuts = []
            for c in candidates:
                point = self.points[c]
                if (common.issubset(self.incidence[c])
                        and not any(cross(direction, sub(point, self.points[a])))):
                    t = dot(sub(point, self.points[a]), direction)/length2
                    if 0 <= t < 1:
                        cuts.append((t, c))
            out.extend(c for _, c in sorted(cuts))
        return out


def _solid_moments(registry, vertices):
    if len(vertices) < 4:
        return 0., np.zeros(3)
    points = [registry.points[i] for i in sorted(vertices)]
    a = points[0]
    if not any(dot(sub(b, a), cross(sub(c, a), sub(d, a)))
               for b, c, d in combinations(points[1:], 3)):
        return 0., np.zeros(3)
    xyz = np.array(points, dtype=float)
    hull = ConvexHull(xyz)
    triangles = xyz[hull.simplices]
    reverse = np.einsum('fi,fi->f', np.cross(triangles[:, 1]-triangles[:, 0],
        triangles[:, 2]-triangles[:, 0]), hull.equations[:, :3]) < 0
    triangles[reverse] = triangles[reverse, ::-1]
    volume = np.einsum('fi,fi->f', triangles[:, 0],
                       np.cross(triangles[:, 1], triangles[:, 2]))/6
    assert abs(volume.sum()-hull.volume) < 1e-12
    return hull.volume, (volume[:, None]*triangles.sum(axis=1)/4).sum(axis=0)


def hull_planes(registry, vertices, faces):
    """Validate authored planar groups AND convexity before using halfspaces."""
    hull = [registry.plane(vertices[faces[group[0]]], ('hull', group)) for group in GROUPS]
    for (plane, _), group in zip(hull, GROUPS):
        assert all(not registry.value(plane, tuple(p)) for p in vertices[faces[list(group)]].reshape(-1, 3)), 'noncoplanar facet group'
    assert all(sign*registry.value(i, tuple(p)) <= 0 for i, sign in hull for p in vertices), 'boat is not convex'
    return hull


def columns(registry, roof, parents, depth=2):
    cells, roofs, outer = [], [], []
    edge_counts = {}
    for parent in parents:
        for a, b in zip(parent, parent[1:]+parent[:1]):
            key = tuple(sorted((a, b)))
            edge_counts[key] = edge_counts.get(key, 0)+1
    for parent_id, parent in enumerate(parents):
        top = registry.plane(roof[list(parent)], ('roof', parent_id))
        bottom_points = roof[list(parent)].copy(); bottom_points[:, 2] = -depth
        bottom = registry.plane(bottom_points[::-1], ('bottom', parent_id))
        sides, exposed = [], [top, bottom]
        for a, b in zip(parent, parent[1:]+parent[:1]):
            p = roof[[a, b, b]].copy(); p[2, 2] = -depth
            side = registry.plane(p[[0, 2, 1]], ('background_edge', tuple(sorted((a, b)))))
            sides.append(side)
            if edge_counts[tuple(sorted((a, b)))] == 1:
                exposed.append(side)
        cells.append([top, bottom, *sides]); roofs.append(top); outer.append(exposed)
    return cells, roofs, outer


def cell_polygons(registry, hull, constraints, exposed, parent_id, body_seen, *, outside_bbox=False):
    """One column's exposed boundary; physical hull groups retain their IDs."""
    polygons, duplicates = [], 0
    if outside_bbox:
        outside = registry.polytope(constraints)
        for boundary, direction in exposed:
            polygon = registry.polygon(outside, boundary, direction)
            if polygon:
                polygons.append((polygon, boundary, direction, 'column', parent_id))
        return polygons, 0., np.zeros(3), 0
    for j, (plane, sign) in enumerate(hull):
        outside = registry.polytope([*constraints, *hull[:j], (plane, -sign)])
        for boundary, direction in exposed:
            polygon = registry.polygon(outside, boundary, direction)
            if polygon:
                polygons.append((polygon, boundary, direction, 'column', parent_id))
    solid = registry.polytope([*constraints, *hull])
    volume, first = _solid_moments(registry, solid)
    for group_id, (plane, sign) in enumerate(hull[:-1]):
        polygon = registry.polygon(solid, plane, -sign)
        if polygon:
            key = plane, tuple(sorted(polygon))
            if key in body_seen:
                duplicates += 1
                continue
            body_seen.add(key)
            polygons.append((polygon, plane, -sign, 'body', group_id))
    return polygons, volume, first, duplicates


def polygon_centroid(points, polygon):
    """Exact area centroid, invariant to extra collinear boundary cutpoints.

    Signed projected triangle areas share one fixed plane scale, so no square
    roots are needed. A vertex average would jump when a vanishing edge adds
    another vertex, changing finite fan triangles at an infinitesimal motion.
    """
    origin = points[polygon[0]]
    triangles = [(points[a], points[b]) for a, b in zip(polygon[1:-1], polygon[2:])]
    normals = [cross(sub(a, origin), sub(b, origin)) for a, b in triangles]
    normal = next(n for n in normals if any(n))
    axis = max(range(3), key=lambda j: abs(normal[j]))
    measure = sum(n[axis] for n in normals)
    return tuple(sum(n[axis]*(origin[j]+a[j]+b[j]) for n, (a, b) in zip(normals, triangles))
                 /(3*measure) for j in range(3))


def merge_body_polygons(registry, polygons, *, protected=(), allow_disconnected=False,
                        require_star=True):
    """Remove artificial submerged column cuts, retaining the exact wet hull.

    By default one star-shaped boundary loop per authored group is supported.
    The material-cell overlay can retain disconnected positive loops explicitly.
    Nonbody-required cutpoints and original hull corners remain shared;
    no proximity merge, convex hull replacement or change of physical coverage.
    """
    required = {i for p, _, _, kind, _ in polygons if kind != 'body' for i in p}|set(protected)
    hull_planes = {i for i, labels in enumerate(registry.labels)
                   if any(label[0] == 'hull' for label in labels)}
    points = registry.points
    pieces = defaultdict(list)
    other = []
    for item in polygons:
        (pieces[item[4]] if item[3] == 'body' else other).append(item)
    removed, merged = 0, []

    def moments(p):
        normal = tuple(sum(cross(points[a], points[b])[j]
                           for a, b in zip(p, p[1:]+p[:1])) for j in range(3))
        center = polygon_centroid(points, p)
        return normal, center

    for group, fragments in sorted(pieces.items()):
        plane, sign = fragments[0][1:3]
        edges = Counter()
        before = [F(0)]*4
        normal, _ = moments(fragments[0][0])
        axis = max(range(3), key=lambda j: abs(normal[j]))
        for polygon, p, s, _, _ in fragments:
            if (p, s) != (plane, sign):
                raise ValueError('original body group has inconsistent physical planes')
            n, c = moments(polygon)
            before[0] += n[axis]
            for j in range(3):
                before[j+1] += n[axis]*c[j]
            edges.update(zip(polygon, polygon[1:]+polygon[:1]))
        boundary = []
        for (a, b), count in edges.items():
            remaining = count-edges[b, a]
            if remaining > 1:
                raise ValueError('overlapping original-body polygon pieces')
            if remaining == 1:
                boundary.append((a, b))
        outgoing, incoming = {}, {}
        for a, b in boundary:
            if a in outgoing or b in incoming:
                raise ValueError('unsupported branching wet-facet boundary')
            outgoing[a], incoming[b] = b, a
        if set(outgoing) != set(incoming) or not outgoing:
            raise ValueError('open original-body facet boundary')
        loops, remaining = [], set(outgoing)
        while remaining:
            loop = [min(remaining)]
            remaining.remove(loop[0])
            while outgoing[loop[-1]] != loop[0]:
                successor = outgoing[loop[-1]]
                if successor not in remaining:
                    raise ValueError('unsupported repeated wet-facet boundary vertex')
                loop.append(successor)
                remaining.remove(successor)
            loops.append(loop)
        if len(loops) > 1 and not allow_disconnected:
            raise ValueError('unsupported multiple wet-facet boundary loops')
        # All incident physical groups simplify a straight shared chine the
        # same way. Actual physical corners carry at least three hull planes.
        after = [F(0)]*4
        for loop in loops:
            changed = True
            while changed:
                changed = False
                for k, vertex in enumerate(loop):
                    if vertex in required or len(registry.incidence[vertex]&hull_planes) >= 3:
                        continue
                    previous, following = points[loop[k-1]], points[loop[(k+1)%len(loop)]]
                    a, b = sub(previous, points[vertex]), sub(following, points[vertex])
                    if not any(cross(a, b)) and dot(a, b) < 0:
                        del loop[k]; removed += 1; changed = True
                        break
            n, center = moments(loop)
            if n[axis]*normal[axis] <= 0:
                raise ValueError('unsupported hole or reversed wet-facet boundary')
            after = [a+b for a, b in zip(after, [n[axis], *(n[axis]*c for c in center)])]
            if require_star:
                for a, b in zip(loop, loop[1:]+loop[:1]):
                    if dot(cross(sub(points[a], center), sub(points[b], center)), n) <= 0:
                        raise ValueError('wet-facet loop is not star-shaped about its area centroid')
            merged.append((loop, plane, sign, 'body', group))
        if before != after:
            raise ValueError('body union changed exact physical area or first moment')
    return other+merged, dict(body_polygons_before=sum(map(len, pieces.values())),
        body_polygons_after=len(merged), removed_collinear_body_points=removed,
        group_area_and_first_moment_preserved_exactly=True)


def polygon_ears(points, polygon, normal):
    """Exact simple-loop ears, retaining collinear shared boundary vertices."""
    loop, result = list(polygon), []
    def side(a, b, c):
        return dot(cross(sub(points[b], points[a]), sub(points[c], points[a])), normal)
    while len(loop) > 3:
        for i, b in enumerate(loop):
            a, c = loop[i-1], loop[(i+1) % len(loop)]
            if side(a, b, c) <= 0:
                continue
            if any(side(a, b, v) >= 0 and side(b, c, v) >= 0 and side(c, a, v) >= 0
                   for v in loop if v not in (a, b, c)):
                continue
            result.append([a, b, c])
            del loop[i]
            break
        else:
            raise ValueError('wet polygon has no nondegenerate exact ear')
    if side(*loop) <= 0:
        raise ValueError('degenerate final wet polygon ear')
    return [*result, loop]


def split_nonstar_polygons(registry, polygons):
    """Keep centroid fans where valid; exact ears partition concave exceptions.

    Retain every boundary vertex, including collinear shared cutpoints. Each
    returned triangle can use the existing centroid-fan/chart machinery.
    """
    points, result = registry.points, []
    for polygon, plane, sign, kind, parent in polygons:
        normal = tuple(sum(cross(points[a], points[b])[j]
                           for a, b in zip(polygon, polygon[1:]+polygon[:1])) for j in range(3))
        center = polygon_centroid(points, polygon)
        if all(dot(cross(sub(points[a], center), sub(points[b], center)), normal) > 0
               for a, b in zip(polygon, polygon[1:]+polygon[:1])):
            result.append((polygon, plane, sign, kind, parent))
            continue
        result.extend((ear, plane, sign, kind, parent)
                      for ear in polygon_ears(points, polygon, normal))
    return result


def triangulate(registry, polygons):
    # Every ACTUAL boundary cutpoint is propagated globally, including points
    # introduced by neighboring parents. Unused infinite-plane intersections
    # are not geometric boundary vertices and need not refine the output.
    candidates = sorted({v for polygon, *_ in polygons for v in polygon})
    polygons = [(registry.split_edges(p, candidates), plane, sign, kind, parent)
                for p, plane, sign, kind, parent in polygons]
    polygons, body_merge = merge_body_polygons(registry, polygons, require_star=False,
                                               allow_disconnected=True)
    merged_count = len(polygons)
    polygons = split_nonstar_polygons(registry, polygons)
    body_merge['additional_nonstar_polygon_pieces'] = len(polygons)-merged_count
    points = [tuple(p) for p in registry.points]
    faces, kinds, source, planes = [], [], [], []
    for polygon, plane, sign, kind, parent in polygons:
        center = polygon_centroid(points, polygon)
        ci = len(points); points.append(center)
        for a, b in zip(polygon, polygon[1:]+polygon[:1]):
            faces.append([ci, a, b]); kinds.append(kind); source.append(parent); planes.append(plane)
    active, inverse = np.unique(faces, return_inverse=True)
    mesh = HullMesh(np.array([points[i] for i in active], dtype=float), inverse.reshape(-1, 3))
    return mesh, dict(active=active, polygons=polygons, points=points,
                      kinds=np.array(kinds), source=np.array(source), planes=np.array(planes),
                      body_merge=body_merge)


def build(heights=(0., 0., 0., 0.)):
    """Closed fluid patch around the boat's actual transom/port-chine corner.

    The two background triangles share x=transom, so the vertical transom patch
    has finite area in BOTH prisms. It must be emitted once, not counted twice.
    Heights are authored P1 graph controls, not a generated or rescaled sea.
    """
    started = perf_counter()
    body = reference_body()
    # Rationalize SOURCE coordinates first, then apply the exact affine pose.
    # Rounding already-transformed vertices independently can destroy planes.
    translation = tuple(F(float(x)) for x in (-body.cg[0], -body.cg[1], -.39))
    vertices = np.array([tuple(F(float(x))+translation[j] for j, x in enumerate(p))
                         for p in body.mesh.vertices], dtype=object)
    xt, yc = vertices[0, 0], vertices[1, 1]
    half = F(1, 2)
    xy = np.array([[xt-half, yc-half], [xt, yc-half], [xt, yc+half], [xt+half, yc+half]], dtype=object)
    roof = np.column_stack((xy, [F(float(x)) for x in heights]))
    parents = ((0, 1, 2), (1, 3, 2))
    registry = PlanePoints()
    hull = hull_planes(registry, vertices, body.mesh.faces)
    cells, roofs, outer = columns(registry, roof, parents)
    polygons, solid_volume, solid_first, column_volume, column_first = [], 0., np.zeros(3), 0., np.zeros(3)
    duplicate_body_faces = 0
    body_seen = set()
    for parent_id, constraints in enumerate(cells):
        pieces, volume, first, duplicate = cell_polygons(registry, hull, constraints,
                                                         outer[parent_id], parent_id, body_seen)
        polygons.extend(pieces)
        duplicate_body_faces += duplicate
        solid_volume += volume; solid_first += first
        points = np.asarray(roof[list(parents[parent_id])], dtype=float)
        area = abs(np.linalg.det(points[1:, :2]-points[:1, :2]))/2
        bary = np.full((3, 3), 1/6); np.fill_diagonal(bary, 2/3)
        sample = bary@points
        column_volume += area*(sample[:, 2]+2).mean()
        column_first += area*np.r_[(sample[:, :2]*(sample[:, 2:]+2)).mean(axis=0),
                                   ((sample[:, 2]**2-4)/2).mean()]
    mesh, meta = triangulate(registry, polygons)
    active, polygons, points = meta['active'], meta['polygons'], meta['points']
    xyz, faces = mesh.vertices, mesh.faces
    triangles = xyz[faces]
    tetra = np.einsum('fi,fi->f', triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2]))/6
    volume = tetra.sum()
    first = (tetra[:, None]*triangles.sum(axis=1)/4).sum(axis=0)
    twist = np.array([.03, -.02, .01, .01, -.015, .02])
    center = np.array([0., 0., body.cg[2]-.39])
    work = normal_work(triangles, center, twist, affine)
    expected = np.array([.2, -.1, .3])@(twist[:3]*(column_volume-solid_volume)
        +np.cross(twist[3:], column_first-solid_first-center*(column_volume-solid_volume)))
    shared_edge = next(i for i, labels in enumerate(registry.labels)
                       if ('background_edge', (1, 2)) in labels)
    shared_roof = {shared_edge, *(i for i, _ in roofs)}
    shared_points = sum(shared_roof.issubset(registry.incidence[i])
                        for i in active if i < len(registry.points))
    return mesh, dict(qualified=False, seconds=perf_counter()-started, nodes=len(xyz), triangles=len(faces),
        body_polygon_merge=meta['body_merge'],
        exact_affine_pose_roundoff_difference=float(np.max(abs(np.asarray(vertices, dtype=float)
            -(body.mesh.vertices+np.array(translation, dtype=float))))),
        exact_vertex_owners=len(registry.vertex_by_owner), duplicate_coincident_body_faces_removed=duplicate_body_faces,
        volume_error=float(volume-(column_volume-solid_volume)),
        first_moment_max_error=float(np.max(abs(first-(column_first-solid_first)))),
        affine_rigid_work_error=float(work-expected),
        shared_roof_edge_owned_points=shared_points,
        emitted_physical_groups=sorted({parent for _, _, _, kind, parent in polygons if kind == 'body'}),
        roof_incidence_exact=all(not registry.value(plane, points[i])
            for polygon, plane, _, kind, _ in polygons if kind == 'column' and plane in {i for i, _ in roofs}
            for i in polygon),
        finite_transom_alias=any(any(label == ('hull', GROUPS[8]) for label in labels)
            and any(label[0] == 'background_edge' for label in labels) for labels in registry.labels),
        parent_count=2, source_physical_triangle_groups=[list(g) for g in GROUPS[:-1]])


def run(output):
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    report = dict(qualified=False, complete=False, cases=[],
        scope='Exact-rational adjacent-parent geometry reference; not full-boat meshing, dynamics or runtime qualification',
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        fixture_sha256=sha256(Path(__file__).resolve().parents[2].joinpath('reference_boat.py').read_bytes()).hexdigest())
    arrays = {}
    try:
        for name, heights in (('flat', (0.,)*4), ('nonflat', (-.01, .005, .025, 0.))):
            mesh, case = build(heights)
            report['cases'].append(dict(name=name, heights=heights, **case))
            arrays[name+'_points'], arrays[name+'_faces'] = mesh.vertices, mesh.faces
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        np.savez(output.with_suffix('.npz'), **arrays)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    print(json.dumps(run(parser.parse_args().output), indent=2))
