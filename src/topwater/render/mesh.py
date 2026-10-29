"""Bounded geometric intersections for projected mesh annotations, not physics."""
import numpy as np
from scipy.spatial import ConvexHull


def projected_footprint(points):
    """Convex XY footprint for an opaque schematic hull interior."""
    xy = np.asarray(points)[:, :2]
    return xy[ConvexHull(xy).vertices]


def mesh_edges(faces):
    """Return unique triangle edges and open-boundary edges, independent of shape."""
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    return unique, unique[counts == 1]


# Compatibility import: producers and renderers share the same contact kernel.
from ..geometry.waterline import waterline_segments
