"""Float integration-reference cuts, never production contact geometry.

The six original P1 parents retain their affine trace. No tolerance welding,
source-mesh replacement, or physical hidden-water state is introduced.
"""
import numpy as np


def _clip(polygon, plane, inside):
    # Same strict crossing convention as boat_spatial_work.clip_polygon;
    # kept dependency-free so this geometry reference imports no solver/GPU.
    result = []
    if not len(polygon):
        return np.empty((0, 2))
    values = polygon @ plane[:2] + plane[2]
    if not inside:
        values = -values
    for a, b, va, vb in zip(polygon, np.roll(polygon, -1, axis=0),
                            values, np.roll(values, -1)):
        if va <= 0:
            result.append(a)
        if (va < 0 < vb) or (vb < 0 < va):
            result.append(a+(b-a)*(va/(va-vb)))
    return np.asarray(result).reshape(-1, 2)


def _planes(polygon):
    edge = np.roll(polygon, -1, axis=0)-polygon
    normal = np.column_stack((edge[:, 1], -edge[:, 0]))
    length = np.linalg.norm(normal, axis=1)
    if np.any(length == 0):
        raise ValueError('duplicate polygon vertex')
    normal /= length[:, None]
    return np.column_stack((normal, -np.sum(normal*polygon, axis=1)))


def _area(polygon):
    if len(polygon) < 3:
        return 0.
    relative = polygon-polygon[0]
    return float(np.sum(relative[:, 0]*np.roll(relative[:, 1], -1)
                        -relative[:, 1]*np.roll(relative[:, 0], -1))/2)


def cut_hat(center, hx, hy, footprint_ccw):
    """Cut one complete uniform BL-to-TR P1 hat by a convex footprint.

    Return interior/exterior positive triangles, parent IDs, affine nodal
    values/gradients, and exact affine triangle moments (in float arithmetic).
    Waterline normals point OUT OF THE HOLE; exterior-domain normals are their
    negatives. Lower-dimensional clipping remnants are explicitly retained.
    """
    center, footprint = np.asarray(center, float), np.asarray(footprint_ccw, float)
    if (center.shape != (2,) or footprint.ndim != 2 or footprint.shape[1] != 2
            or len(footprint) < 3 or not np.isfinite(center).all()
            or not np.isfinite(footprint).all() or not np.isfinite([hx, hy]).all()
            or min(hx, hy) <= 0):
        raise ValueError('finite center, positive spacings, convex CCW footprint required')
    edge = np.roll(footprint, -1, axis=0)-footprint
    following = np.roll(edge, -1, axis=0)
    if np.any(edge[:, 0]*following[:, 1]-edge[:, 1]*following[:, 0] <= 0):
        raise ValueError('footprint must be strictly convex and CCW')
    planes = _planes(footprint)
    # Local turns alone do not exclude winding/self-intersection.
    if np.any(footprint @ planes[:, :2].T + planes[:, 2] >
              64*np.finfo(float).eps*max(1., np.max(np.abs(footprint)))):
        raise ValueError('footprint violates its convex supporting planes')
    nodes = np.array([[x*hx, y*hy] for y in (-1, 0, 1) for x in (-1, 0, 1)])+center
    faces = []
    for row in range(2):
        for col in range(2):
            a = row*3+col
            for ids in ((a, a+1, a+4), (a, a+4, a+3)):
                if 4 in ids:
                    faces.append(ids)
    faces = np.asarray(faces)
    parents = nodes[faces]
    parent_values = (faces == 4).astype(float)
    gradients = np.linalg.solve(parents[:, 1:]-parents[:, :1],
        (parent_values[:, 1:]-parent_values[:, :1])[..., None])[..., 0]
    parts = {'interior': [], 'exterior': []}
    tangencies = []

    def store(kind, polygon, parent):
        if not len(polygon):
            return
        area = _area(polygon)
        if area < 0:
            raise ValueError('clipping reversed a polygon')
        if area == 0:
            tangencies.append(dict(kind=kind, parent_id=parent, points=polygon.copy()))
            return
        middle = polygon.mean(axis=0)
        for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
            triangle = np.array([middle, a, b])
            if _area(triangle) <= 0:
                raise ValueError('degenerate positive-polygon fan; no vertex deletion allowed')
            parts[kind].append((parent, triangle))

    for parent, triangle in enumerate(parents):
        remaining = triangle
        for plane in planes:
            store('exterior', _clip(remaining, plane, False), parent)
            remaining = _clip(remaining, plane, True)
        store('interior', remaining, parent)

    def pack(records):
        ids = np.array([p for p, _ in records], dtype=int)
        triangles = np.array([t for _, t in records]).reshape(-1, 3, 2)
        values = parent_values[ids, :1]+np.sum(
            (triangles-parents[ids, :1])*gradients[ids, None, :], axis=2)
        areas = np.array([_area(t) for t in triangles])
        return dict(triangles=triangles, parent_id=ids, values=values,
                    gradients=gradients[ids], areas=areas,
                    hat_moments=areas*values.mean(axis=1))

    interior, exterior = pack(parts['interior']), pack(parts['exterior'])
    exterior_area = np.bincount(exterior['parent_id'], weights=exterior['areas'], minlength=6)
    segments, segment_ids, normals = [], [], []
    segment_owners = {}
    for parent, triangle in enumerate(parents):
        if exterior_area[parent] == 0:
            continue
        for edge_id, (a, b) in enumerate(zip(footprint, np.roll(footprint, -1, axis=0))):
            low, high = 0., 1.
            for plane in _planes(triangle):
                value, slope = a @ plane[:2]+plane[2], (b-a) @ plane[:2]
                if slope == 0:
                    if value > 0:
                        high = -1.
                        break
                elif slope > 0:
                    high = min(high, -value/slope)
                else:
                    low = max(low, -value/slope)
            if low < high:
                segment = np.array([a+low*(b-a), a+high*(b-a)])
                if np.linalg.norm(segment[1]-segment[0]) == 0:
                    raise ValueError('nonempty clipped segment collapsed')
                key = tuple(sorted(map(tuple, segment)))
                if key in segment_owners:
                    raise ValueError('ambiguous duplicate waterline segment ownership')
                segment_owners[key] = parent
                segments.append(segment)
                segment_ids.append(parent)
                normals.append(planes[edge_id, :2])
            elif low == high:
                tangencies.append(dict(kind='waterline', parent_id=parent,
                    points=(a+low*(b-a))[None, :]))
    segment_ids = np.asarray(segment_ids, int)
    segments = np.asarray(segments).reshape(-1, 2, 2)
    segment_values = parent_values[segment_ids, :1]+np.sum(
        (segments-parents[segment_ids, :1])*gradients[segment_ids, None, :], axis=2)
    return dict(parents=parents, parent_values=parent_values, parent_gradients=gradients,
        interior=interior, exterior=exterior, tangencies=tangencies,
        waterline=dict(segments=segments, parent_id=segment_ids,
            hole_outward_normals=np.asarray(normals).reshape(-1, 2),
            values=segment_values, gradients=gradients[segment_ids]),
        physical_hat_mass=float(exterior['hat_moments'].sum()),
        scope='Float integration reference; no certified moving-contact topology guard.')
