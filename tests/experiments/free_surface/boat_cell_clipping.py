"""Exact known-column clipping falsifier; the reference remains the default.

Only full-dimensional convex columns are incrementally clipped. Nonempty
lower-dimensional results use the original rational enumerator, including its
point/line/planar ownership. No provisional vertex enters the global registry.
"""
from dataclasses import dataclass
from itertools import combinations

import numpy as np

from .boat_parent_clip_proof import cross, dot, sub, _solid_moments


def _value(registry, plane, point):
    row = registry.planes[plane]
    return dot(row[:3], point)-row[3]


def _full_dimension(points):
    if len(points) < 4:
        return False
    origin = points[0]
    directions = [sub(point, origin) for point in points[1:]]
    first = next((entry for entry in directions if any(entry)), None)
    if first is None:
        return False
    normal = next((cross(first, entry) for entry in directions if any(cross(first, entry))), None)
    return normal is not None and any(dot(normal, entry) for entry in directions)


def _intersection(registry, triple):
    """Local rational construction, without cache/registry side effects."""
    a, b, c = [registry.planes[index] for index in triple]
    determinant = dot(a[:3], cross(b[:3], c[:3]))
    if not determinant:
        return None
    vectors = cross(b[:3], c[:3]), cross(c[:3], a[:3]), cross(a[:3], b[:3])
    return tuple(sum(row[3]*vector[k] for row, vector in zip((a, b, c), vectors))/determinant
                 for k in range(3))


@dataclass(frozen=True)
class Column:
    points: tuple
    constraints: tuple


def known_column(registry, constraints):
    """Six vertices from top/bottom crossed with three vertical side pairs.

    Return None for unsupported/degenerate initial cells, causing reference
    fallback. The existing columns() owner supplies this exact ordering.
    """
    if len(constraints) != 5 or any(sign not in (-1, 1) for _, sign in constraints):
        return None
    top, bottom = [registry.planes[index] for index, _ in constraints[:2]]
    if (constraints[0][1]*top[2] <= 0 or constraints[1][1]*bottom[2] >= 0
            or any(bottom[:2]) or any(registry.planes[index][2] for index, _ in constraints[2:])):
        return None
    sides = [index for index, _ in constraints[2:]]
    points = []
    for cap, _ in constraints[:2]:
        for pair in combinations(sides, 2):
            point = _intersection(registry, (cap, *pair))
            if point is None or any(sign*_value(registry, index, point) > 0 for index, sign in constraints):
                return None
            points.append(point)
    points = tuple(sorted(set(points)))
    if len(points) != 6 or not _full_dimension(points):
        return None
    return Column(points, tuple(constraints))


def clip(registry, column, constraint):
    """Intersect a full-dimensional bounded column with one closed halfspace.

    None denotes a nonempty lower-dimensional/unsupported result. Empty tuples
    remain empty; they are not silently substituted for tangent contacts.
    """
    if column is None:
        return None
    plane, sign = constraint
    if sign not in (-1, 1):
        raise ValueError('halfspace sign must be plus or minus one')
    constraints = (*column.constraints, constraint)
    if not column.points:
        return Column((), constraints)
    values = {point: sign*_value(registry, plane, point) for point in column.points}
    kept = {point for point, value in values.items() if value <= 0}
    if all(value <= 0 for value in values.values()):
        return Column(column.points, constraints)
    active = {point: {index for index, _ in column.constraints if not _value(registry, index, point)}
              for point in column.points}
    for a, b in combinations(column.points, 2):
        if values[a]*values[b] >= 0:
            continue
        common = sorted(active[a]&active[b])
        if not any(any(cross(registry.planes[i][:3], registry.planes[j][:3]))
                   for i, j in combinations(common, 2)):
            continue
        fraction = values[a]/(values[a]-values[b])
        kept.add(tuple(x+fraction*(y-x) for x, y in zip(a, b)))
    points = tuple(sorted(kept))
    if points and not _full_dimension(points):
        return None
    return Column(points, constraints)


def register_column(registry, column, constraints):
    """Preserve the legacy accepted-triple registration and relevance order."""
    if column is None:
        return registry.polytope(constraints)
    candidates = []
    for point in column.points:
        if any(sign*_value(registry, index, point) > 0 for index, sign in constraints):
            raise ValueError('clipped point violates its exact halfspace constraints')
        active = sorted({index for index, _ in constraints if not _value(registry, index, point)})
        for triple in combinations(active, 3):
            a, b, c = [registry.planes[index][:3] for index in triple]
            if dot(a, cross(b, c)):
                candidates.append((triple, point))
    vertices = set()
    for triple, point in sorted(candidates):
        # Reuse the rational source owner for accepted vertices only. This also
        # preserves cache freeze semantics and independently checks interpolation.
        exact = registry.intersection(triple)
        if exact != point:
            raise ValueError('incremental point differs from rational triple owner')
        vertices.add(registry.register(triple, exact))
    if len(vertices) != len(column.points):
        raise ValueError('clipped vertex lacks a nonsingular active owner triple')
    if registry.relevant is not None:
        for vertex in vertices:
            registry.relevant[vertex].update(index for index, _ in constraints)
    return vertices


def cell_polygons(registry, hull, constraints, exposed, parent_id, body_seen, *, outside_bbox=False):
    """Explicit adapter for ONE existing cell_polygons owner, not a mesher."""
    column = known_column(registry, constraints)
    polygons, duplicates = [], 0
    if outside_bbox:
        outside = register_column(registry, column, constraints)
        for boundary, direction in exposed:
            polygon = registry.polygon(outside, boundary, direction)
            if polygon:
                polygons.append((polygon, boundary, direction, 'column', parent_id))
        return polygons, 0., np.zeros(3), 0
    inside = column
    for j, (plane, sign) in enumerate(hull):
        current = [*constraints, *hull[:j], (plane, -sign)]
        outside = register_column(registry, clip(registry, inside, (plane, -sign)), current)
        for boundary, direction in exposed:
            polygon = registry.polygon(outside, boundary, direction)
            if polygon:
                polygons.append((polygon, boundary, direction, 'column', parent_id))
        inside = clip(registry, inside, (plane, sign))
    solid = register_column(registry, inside, [*constraints, *hull])
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
