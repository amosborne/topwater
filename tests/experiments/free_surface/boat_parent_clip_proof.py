"""Exact-arithmetic TWO-parent boat clipping reference, not a hot-path mesher.

This bounded code proof exercises a real chine and a transom coincident with
the shared vertical parent face. Plane/vertex identities are exact incidence
sets, never proximity welds. No volume mesh, fluid solve or evolved state.
"""
from fractions import Fraction as F
from itertools import combinations
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
    def __init__(self):
        self.planes, self.labels, self.cache = [], [], {}
        self.points, self.incidence, self.vertex_by_owner = [], [], {}

    def plane(self, points, label):
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

    def intersection(self, triple):
        triple = tuple(sorted(triple))
        if triple in self.cache:
            return self.cache[triple]
        a, b, c = (self.planes[i] for i in triple)
        determinant = dot(a[:3], cross(b[:3], c[:3]))
        if not determinant:
            self.cache[triple] = None
            return None
        vectors = cross(b[:3], c[:3]), cross(c[:3], a[:3]), cross(a[:3], b[:3])
        point = tuple(sum(plane[3]*v[j] for plane, v in zip((a, b, c), vectors))/determinant
                      for j in range(3))
        # Complete incidence promotes roof-edge and hull-vertex degeneracies to
        # ONE owner regardless of the starting parent or intersection triple.
        owner = tuple(i for i in range(len(self.planes)) if not self.value(i, point))
        if owner not in self.vertex_by_owner:
            self.vertex_by_owner[owner] = len(self.points)
            self.points.append(point)
            self.incidence.append(frozenset(owner))
        index = self.vertex_by_owner[owner]
        assert self.points[index] == point
        self.cache[triple] = index
        return index

    def polytope(self, constraints):
        vertices = set()
        for triple in combinations(sorted({i for i, _ in constraints}), 3):
            index = self.intersection(triple)
            if index is not None and all(sign*self.value(i, self.points[index]) <= 0
                                         for i, sign in constraints):
                vertices.add(index)
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

    def split_edges(self, polygon):
        out = []
        for a, b in zip(polygon, polygon[1:]+polygon[:1]):
            common = self.incidence[a]&self.incidence[b]
            direction = sub(self.points[b], self.points[a])
            length2 = dot(direction, direction)
            cuts = []
            for c, point in enumerate(self.points):
                if common.issubset(self.incidence[c]):
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
    hull = [registry.plane(vertices[body.mesh.faces[group[0]]], ('hull', group)) for group in GROUPS]
    raw = list(map(tuple, vertices))
    assert all(sign*registry.value(i, p) <= 0 for i, sign in hull for p in raw), 'boat is not convex'
    cells, roofs, outer = [], [], []
    edge_counts = {}
    for parent in parents:
        for a, b in zip(parent, parent[1:]+parent[:1]):
            edge_counts[tuple(sorted((a, b)))] = edge_counts.get(tuple(sorted((a, b))), 0)+1
    for parent_id, parent in enumerate(parents):
        top = registry.plane(roof[list(parent)], ('roof', parent_id))
        bottom_points = roof[list(parent)].copy(); bottom_points[:, 2] = -2
        bottom = registry.plane(bottom_points[::-1], ('bottom', parent_id))
        sides, exposed = [], [top, bottom]
        for a, b in zip(parent, parent[1:]+parent[:1]):
            p = roof[[a, b, b]].copy(); p[2, 2] = -2
            side = registry.plane(p[[0, 2, 1]], ('background_edge', tuple(sorted((a, b)))))
            sides.append(side)
            if edge_counts[tuple(sorted((a, b)))] == 1:
                exposed.append(side)
        cells.append([top, bottom, *sides]); roofs.append(top); outer.append(exposed)
    polygons, solid_volume, solid_first, column_volume, column_first = [], 0., np.zeros(3), 0., np.zeros(3)
    duplicate_body_faces = 0
    body_seen = set()
    for parent_id, constraints in enumerate(cells):
        for j, (plane, sign) in enumerate(hull):
            # First violated hull plane partitions the cell exterior disjointly.
            outside = registry.polytope([*constraints, *hull[:j], (plane, -sign)])
            for boundary, direction in outer[parent_id]:
                polygon = registry.polygon(outside, boundary, direction)
                if polygon:
                    polygons.append((polygon, boundary, direction, 'column', parent_id))
        solid = registry.polytope([*constraints, *hull])
        volume, first = _solid_moments(registry, solid)
        solid_volume += volume; solid_first += first
        for group_id, (plane, sign) in enumerate(hull[:-1]):
            polygon = registry.polygon(solid, plane, -sign)
            if polygon:
                key = plane, tuple(sorted(polygon))
                if key in body_seen:
                    duplicate_body_faces += 1
                    continue  # Explicit one-parent ownership, same physical face.
                body_seen.add(key)
                polygons.append((polygon, plane, -sign, 'body', group_id))
        points = np.asarray(roof[list(parents[parent_id])], dtype=float)
        area = abs(np.linalg.det(points[1:, :2]-points[:1, :2]))/2
        bary = np.full((3, 3), 1/6); np.fill_diagonal(bary, 2/3)
        sample = bary@points
        column_volume += area*(sample[:, 2]+2).mean()
        column_first += area*np.r_[(sample[:, :2]*(sample[:, 2:]+2)).mean(axis=0),
                                   ((sample[:, 2]**2-4)/2).mean()]
    # All global intersections exist BEFORE any edge subdivision/triangulation.
    polygons = [(registry.split_edges(p), plane, sign, kind, parent)
                for p, plane, sign, kind, parent in polygons]
    points = [tuple(p) for p in registry.points]
    faces, kinds, source = [], [], []
    for polygon, plane, sign, kind, parent in polygons:
        center = tuple(sum(points[i][j] for i in polygon)/len(polygon) for j in range(3))
        ci = len(points); points.append(center)
        for a, b in zip(polygon, polygon[1:]+polygon[:1]):
            faces.append([ci, a, b]); kinds.append(kind); source.append(parent)
    active, inverse = np.unique(faces, return_inverse=True)
    xyz = np.array([points[i] for i in active], dtype=float)
    faces = inverse.reshape(-1, 3)
    mesh = HullMesh(xyz, faces)
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
