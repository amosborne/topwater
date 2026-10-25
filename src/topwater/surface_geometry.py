"""Shared physical-hull contact geometry; never trim or change saved water.

Water vertices inside a posed hull remain in the rectangular producer mesh but
are outside its physical water-air domain. The independent intersection contour
uses physical hull faces only. Artificial opening caps never define waterlines.
"""

import numpy as np

from ._validation import array, indices
from .inputs import HullMesh


def _mesh(vertices, faces, name):
    points = array(vertices, (None, 3), name+' points')
    faces = np.asarray(faces)
    if (faces.ndim != 2 or faces.shape[1:] != (3,) or faces.dtype.kind not in 'iu'
            or not len(faces) or np.any(faces < 0) or np.any(faces >= len(points))):
        raise ValueError(name+' faces must be nonempty integer triangles with valid indices')
    return points, faces


def _cuts(triangles, distances, tolerance):
    """Vertices/edge crossings where triangles meet a normalized plane."""
    other = np.roll(triangles, -1, axis=1)
    d_other = np.roll(distances, -1, axis=1)
    crossing = ((distances < -tolerance) & (d_other > tolerance)) | ((distances > tolerance) & (d_other < -tolerance))
    fraction = np.divide(distances, distances - d_other, out=np.zeros_like(distances), where=crossing)
    points = np.concatenate([triangles, triangles + fraction[:, :, None] * (other - triangles)], axis=1)
    valid = np.concatenate([abs(distances) <= tolerance, crossing], axis=1)
    return points, valid


def waterline_segments(hull_points, hull_faces, water_points, water_faces, *, artificial_faces=()):
    """Intersect actual physical hull triangles with a piecewise-planar water mesh.

    Each candidate pair is cut by the other's plane. Overlap of the resulting
    collinear intervals gives its intersection segment. Candidate arrays are
    bounded to one physical face at a time, never faces times water faces.
    Point contacts have zero line length. Coplanar area contact is ambiguous and
    rejected explicitly instead of inventing a waterline around a fake cap.
    hull_faces must be physical, or artificial_faces must identify the internal
    cap triangles to exclude. No water vertices, values or face indices change.
    """
    hull_points, hull_faces = _mesh(hull_points, hull_faces, 'hull')
    water_points, water_faces = _mesh(water_points, water_faces, 'water')
    artificial = indices(artificial_faces, len(hull_faces), 'artificial_faces')
    mask = np.ones(len(hull_faces), dtype=bool)
    mask[artificial] = False
    hull = np.asarray(hull_points)[hull_faces[mask]]
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



def validate_opening_clearance(water_points, water_faces, hull_points, hull_faces, *,
                               opening_rim=(), artificial_faces=(), tolerance_m=None):
    """Retain existing fatal open-boat compatibility checks, without remeshing.

    Coordinates are world metres. Supply the full posed outward hull including
    its INTERNAL cap and paired rim/cap indices. Exact rim/cap contact is fatal,
    including point/edge tangency. Complete rim-edge coverage must be locally
    single-valued and dry. A sealed hull supplies neither rim nor cap indices.
    This checks local opening compatibility, not a new fluid boundary condition.
    Input vertex values and face indices are never changed or returned altered.
    """
    water, water_faces = _mesh(water_points, water_faces, 'water')
    hull = HullMesh(hull_points, hull_faces)
    origin = hull.vertices.mean(axis=0)
    scale = float(np.max(np.ptp(hull.vertices, axis=0)))
    eps = 1e-10 if tolerance_m is None else float(tolerance_m)/scale
    if not np.isfinite(eps) or eps <= 0:
        raise ValueError('opening clearance tolerance_m must be positive and finite')
    rim = indices(opening_rim, len(hull.vertices), 'opening_rim')
    artificial = indices(artificial_faces, len(hull.faces), 'artificial_faces')
    if bool(len(rim)) != bool(len(artificial)) or (len(rim) and len(rim) < 3):
        raise ValueError('open hull compatibility needs opening_rim and artificial_faces together')
    if len(rim):
        HullMesh(hull.vertices, hull.faces, opening_rim=rim, artificial_faces=artificial)
    normalized_hull = (hull.vertices-origin)/scale
    normalized_water = (water-origin)/scale
    _opening_clearance(normalized_water, water_faces, normalized_hull, rim, eps)
    _cap_contact(normalized_water, water_faces, normalized_hull, hull.faces[artificial], eps)
