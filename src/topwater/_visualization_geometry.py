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


def _cuts(triangles, distances, tolerance):
    """Vertices/edge crossings where triangles meet a normalized plane."""
    other = np.roll(triangles, -1, axis=1)
    d_other = np.roll(distances, -1, axis=1)
    crossing = ((distances < -tolerance) & (d_other > tolerance)) | ((distances > tolerance) & (d_other < -tolerance))
    fraction = np.divide(distances, distances - d_other, out=np.zeros_like(distances), where=crossing)
    points = np.concatenate([triangles, triangles + fraction[:, :, None] * (other - triangles)], axis=1)
    valid = np.concatenate([abs(distances) <= tolerance, crossing], axis=1)
    return points, valid


def waterline_segments(hull_points, hull_faces, water_points, water_faces):
    """Intersect actual physical hull triangles with a piecewise-planar water mesh.

    Each candidate pair is cut by the other's plane. Overlap of the resulting
    collinear intervals gives its intersection segment. Candidate arrays are
    bounded to one physical face at a time, never faces times water faces.
    Point contacts have zero line length. Coplanar area contact is ambiguous and
    rejected explicitly instead of inventing a waterline around a fake cap.
    """
    hull = np.asarray(hull_points)[hull_faces]
    water = np.asarray(water_points)[water_faces]
    lower, upper = water.min(axis=1), water.max(axis=1)
    scale = max(float(np.linalg.norm(np.ptp(hull_points, axis=0))), 1e-12)
    tolerance = scale * 1e-9
    output = []
    for triangle in hull:
        selected = np.all((upper >= triangle.min(axis=0) - tolerance) & (lower <= triangle.max(axis=0) + tolerance), axis=1)
        other = water[selected]
        if not len(other):
            continue
        normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
        normal /= np.linalg.norm(normal)
        other_normal = np.cross(other[:, 1] - other[:, 0], other[:, 2] - other[:, 0])
        magnitude = np.linalg.norm(other_normal, axis=1)
        if np.any(magnitude <= 1e-15 * scale**2):
            raise ValueError('Degenerate water triangle in hull contact region')
        other_normal /= magnitude[:, None]
        water_distances = (other - triangle[0]) @ normal
        direction = np.cross(normal, other_normal)
        length = np.linalg.norm(direction, axis=1)
        parallel = length <= 1e-12
        if np.any(parallel & (abs(water_distances).max(axis=1) <= tolerance)):
            raise ValueError('Coplanar hull/water area contact has no unique waterline')
        valid_pair = ~parallel
        other, other_normal = other[valid_pair], other_normal[valid_pair]
        water_distances = water_distances[valid_pair]
        direction = direction[valid_pair] / length[valid_pair, None]
        if not len(other):
            continue
        hull_copies = np.broadcast_to(triangle, other.shape)
        hull_distances = np.einsum('kvi,ki->kv', hull_copies - other[:, :1], other_normal)
        hp, hv = _cuts(hull_copies, hull_distances, tolerance)
        wp, wv = _cuts(other, water_distances, tolerance)
        hs = np.einsum('kvi,ki->kv', hp, direction)
        ws = np.einsum('kvi,ki->kv', wp, direction)
        lo = np.maximum(np.where(hv, hs, np.inf).min(axis=1), np.where(wv, ws, np.inf).min(axis=1))
        hi = np.minimum(np.where(hv, hs, -np.inf).max(axis=1), np.where(wv, ws, -np.inf).max(axis=1))
        keep = (hi - lo > tolerance) & np.isfinite(lo) & np.isfinite(hi)
        if np.any(keep):
            reference = wp[np.arange(len(wp)), np.argmax(wv, axis=1)][keep]
            along = np.einsum('ki,ki->k', reference, direction[keep])
            output.append(reference[:, None] + (np.stack([lo[keep], hi[keep]], axis=1) - along[:, None])[:, :, None] * direction[keep, None])
    return np.concatenate(output) if output else np.empty((0, 2, 3))
