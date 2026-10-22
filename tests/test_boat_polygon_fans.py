"""Exact concave wet coverage with shared collinear boundary cutpoints."""
from fractions import Fraction as F
from types import SimpleNamespace

from .experiments.free_surface.boat_parent_clip_proof import (
    split_nonstar_polygons, polygon_centroid, cross)


def test_concave_partition_preserves_area_moment_and_boundary_cutpoints():
    xy = [(0, 0), (F(3, 2), 0), (3, 0), (3, 3), (2, 3),
          (2, 1), (1, 1), (1, 3), (0, 3)]
    points = [tuple(map(F, (x, y, 0))) for x, y in xy]
    polygon = list(range(len(points)))
    pieces = split_nonstar_polygons(SimpleNamespace(points=points), [(polygon, 0, 1, 'body', 4)])
    assert len(pieces) == len(points)-2
    edges, area, first = {}, F(0), [F(0)]*3
    for loop, plane, sign, kind, parent in pieces:
        assert (plane, sign, kind, parent) == (0, 1, 'body', 4)
        measure = sum(cross(points[a], points[b])[2] for a, b in zip(loop, loop[1:]+loop[:1]))/2
        assert measure > 0
        center = polygon_centroid(points, loop)
        area += measure
        first = [a+measure*b for a, b in zip(first, center)]
        for a, b in zip(loop, loop[1:]+loop[:1]):
            edges[a, b] = edges.get((a, b), 0)+1
    assert area == 7
    assert first == [F(21, 2), F(19, 2), 0]
    boundary = {(a, b) for (a, b), count in edges.items() if count-edges.get((b, a), 0) == 1}
    assert boundary == set(zip(polygon, polygon[1:]+polygon[:1]))
