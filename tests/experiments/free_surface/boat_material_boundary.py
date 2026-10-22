"""Exact material-cell overlay for the canonical hull; cold reference only.

Fixed material triangles replace whole-facet moving centroid fans. World-parent
cuts are removed inside each material cell; actual waterline intersections are
shared with the free boundary. No coordinate welding or physical hull change.
The floating-point guarded template deliberately does not support this yet.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from time import perf_counter

import numpy as np

from topwater import HullMesh
from .boat_parent_clip_proof import (GROUPS, sub, cross, dot, polygon_centroid,
                                     merge_body_polygons, polygon_ears)


def material_cells(source, faces, refinement):
    """Conforming dyadic subdivision in ORIGINAL exact source coordinates."""
    cells = []
    for group, ids in enumerate(GROUPS[:-1]):
        for source_id in ids:
            triangles = [tuple(tuple(source[i]) for i in faces[source_id])]
            for _ in range(refinement):
                children = []
                for a, b, c in triangles:
                    ab, bc, ca = [tuple((u+v)/2 for u, v in zip(x, y))
                                  for x, y in ((a, b), (b, c), (c, a))]
                    children.extend(((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))
                triangles = children
            cells.extend((group, triangle) for triangle in triangles)
    return cells


def _clip(polygon, triangle):
    """Convex exact polygon intersected with one oriented exact triangle."""
    normal = cross(sub(triangle[1], triangle[0]), sub(triangle[2], triangle[0]))
    result = polygon
    for a, b in zip(triangle, (*triangle[1:], triangle[0])):
        if len(result) < 3:
            return []
        edge = sub(b, a)
        values = [dot(cross(edge, sub(p, a)), normal) for p in result]
        clipped = []
        for p, q, vp, vq in zip(result, result[1:]+result[:1], values, values[1:]+values[:1]):
            if vp >= 0:
                clipped.append(p)
            if (vp < 0 < vq) or (vq < 0 < vp):
                t = vp/(vp-vq)
                clipped.append(tuple(x+t*(y-x) for x, y in zip(p, q)))
        result = []
        for p in clipped:
            if not result or result[-1] != p:
                result.append(p)
        if len(result) > 1 and result[-1] == result[0]:
            result.pop()
    if len(result) < 3:
        return []
    area = tuple(sum(cross(p, q)[j] for p, q in zip(result, result[1:]+result[:1]))
                 for j in range(3))
    return result if any(area) else []


def triangulate_material(registry, polygons, source, source_faces, rr, cg, position, refinement):
    """Return the existing metadata contract, with retained material cell IDs."""
    if isinstance(refinement, bool) or not isinstance(refinement, int) or refinement not in (0, 1, 2, 3):
        raise ValueError('bounded material overlay supports levels 0..3')
    if registry.relevant is not None:
        raise ValueError('material overlay has no guarded template export yet')
    started = perf_counter()

    def budget():
        if perf_counter()-started > 30:
            raise RuntimeError('30 s material overlay budget exhausted between cells')

    def world(p):
        return tuple(sum(rr[j][k]*(p[k]-cg[k]) for k in range(3))+position[j] for j in range(3))

    cells = [(group, tuple(world(p) for p in triangle))
             for group, triangle in material_cells(source, source_faces, refinement)]
    by_group = defaultdict(list)
    for item in polygons:
        if item[3] == 'body':
            by_group[item[4]].append(item)
    # New nodes are registered by exact coordinates, never by rounded proximity.
    # No further plane intersections/polytope registrations follow this overlay.
    lookup = {tuple(p): i for i, p in enumerate(registry.points)}

    def register(point):
        point = tuple(point)
        if point not in lookup:
            if registry.max_points is not None and len(registry.points) >= registry.max_points:
                raise RuntimeError('reference material point cap reached')
            lookup[point] = len(registry.points)
            registry.points.append(point)
            registry.incidence.append(frozenset(i for i in range(len(registry.planes))
                                               if not registry.value(i, point)))
        return lookup[point]

    protected = set()
    material_ids = []
    fragments = [item for item in polygons if item[3] != 'body']
    for cell, (group, triangle) in enumerate(cells):
        budget()
        ids = tuple(register(p) for p in triangle)
        material_ids.append(ids)
        protected.update(ids)
        for polygon, plane, sign, kind, _ in by_group[group]:
            # Both pieces are convex; do this BEFORE merging across parents.
            clipped = _clip([registry.points[i] for i in polygon], triangle)
            if clipped:
                fragments.append(([register(p) for p in clipped], plane, sign, kind, cell))
    # Actual new seam nodes propagate into the free surface and neighboring
    # material cells. A physical plane alone no longer uniquely defines an edge,
    # so split_edges also checks exact collinearity.
    candidates = sorted({v for polygon, *_ in fragments for v in polygon})
    on_plane = defaultdict(list)
    for vertex in candidates:
        for plane in registry.incidence[vertex]:
            on_plane[plane].append(vertex)
    split = []
    for polygon, plane, sign, kind, owner in fragments:
        budget()
        split.append((registry.split_edges(polygon, on_plane[plane]), plane, sign, kind, owner))
    merged, merge_report = merge_body_polygons(registry, split, protected=protected,
                                                allow_disconnected=True, require_star=False)
    component_counts = Counter(owner for _, _, _, kind, owner in merged if kind == 'body')
    # Exact per-group area and first moments must match BEFORE triangulation.
    def moments(items, use_cells):
        result = defaultdict(lambda: [F(0)]*12)
        for polygon, _, _, kind, owner in items:
            if kind != 'body':
                continue
            group = cells[owner][0] if use_cells else owner
            p = [registry.points[i] for i in polygon]
            n = tuple(sum(cross(a, b)[j] for a, b in zip(p, p[1:]+p[:1])) for j in range(3))
            center = polygon_centroid(registry.points, polygon)
            values = [*n, *(n[i]*center[j] for i in range(3) for j in range(3))]
            result[group] = [a+b for a, b in zip(result[group], values)]
        return dict(result)
    if moments(polygons, False) != moments(merged, True):
        raise ValueError('material overlay changed exact physical body area/first moments')
    points = [tuple(p) for p in registry.points]
    faces, kinds, sources, planes, cell_ids = [], [], [], [], []
    fully_wet, direct, nonstar = [], [], []
    for polygon, plane, sign, kind, owner in merged:
        budget()
        if kind == 'body' and set(polygon) == set(material_ids[owner]):
            # Exact coverage of the entire material triangle; three immutable
            # corners, no moving centroid or inherited world-parent cuts.
            triangles = [polygon]
            fully_wet.append(owner)
        elif kind == 'body' and len(polygon) == 3:
            triangles = [polygon]
        else:
            center = polygon_centroid(points, polygon)
            normal = tuple(sign*x for x in registry.planes[plane][:3])
            star = all(dot(cross(sub(points[a], center), sub(points[b], center)), normal) > 0
                       for a, b in zip(polygon, polygon[1:]+polygon[:1]))
            if not star and kind == 'body':
                triangles = polygon_ears(points, polygon, normal)
                nonstar.append(owner)
            else:
                if not star:
                    raise ValueError('nonbody reference polygon lost its interior fan')
                ci = len(points)
                points.append(center)
                triangles = [[ci, a, b] for a, b in zip(polygon, polygon[1:]+polygon[:1])]
        if kind == 'body' and len(triangles) == 1:
            direct.append(owner)
        for triangle in triangles:
            faces.append(triangle)
            kinds.append(kind)
            sources.append(cells[owner][0] if kind == 'body' else owner)
            planes.append(plane)
            if kind == 'body':
                cell_ids.append(owner)
    active, inverse = np.unique(faces, return_inverse=True)
    mesh = HullMesh(np.array([points[i] for i in active], dtype=float), inverse.reshape(-1, 3))
    return mesh, dict(active=active, polygons=merged, points=points, kinds=np.array(kinds),
        source=np.array(sources), planes=np.array(planes), body_merge=merge_report,
        body_material_cell=np.array(cell_ids), material_cells=cells,
        fully_wet_material_cells=np.array(fully_wet, dtype=int),
        material_report=dict(refinement=refinement, seconds=perf_counter()-started,
            authored_material_cells=len(cells), fully_wet_cells=len(fully_wet),
            direct_body_triangles=len(direct), nonstar_fringe_loops=len(nonstar),
            disconnected_material_cells=sum(count > 1 for count in component_counts.values()),
            per_group_area_and_first_moments_exact=True))
