"""Strict physical hull ingestion for the convex, single-opening workflow.

The file mesh contains physical panels only. An inferred opening is
closed internally for the existing HullMesh API; that cap is never physical.
This is intentionally narrower than the generic mathematical HullMesh class.
"""

import numpy as np

from .._validation import array
from ..inputs import HullMesh


class UnsupportedClosureError(ValueError):
    """Valid physical hull needs an internal cap operation not yet supported."""


def _envelope_cap(points, faces, normals, eps):
    """Cap absent whole convex support facets; do not invent polygon booleans."""
    from scipy.spatial import ConvexHull
    from ..geometry.triangulation import _ears

    envelope = ConvexHull(points)
    planes = []
    for equation in envelope.equations:
        if not any(np.max(abs(equation-plane)) <= 8*eps for plane in planes):
            planes.append(equation)
    triangles = points[faces]
    areas = np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0],
                                   triangles[:, 2]-triangles[:, 0]), axis=1)/2
    cap = []
    for plane in planes:
        normal, offset = plane[:3], plane[3]
        on = np.flatnonzero(abs(points @ normal+offset) <= eps)
        axis = np.eye(3)[np.argmin(abs(normal))]
        u = np.cross(axis, normal); u /= np.linalg.norm(u)
        v = np.cross(normal, u)
        projected = points @ np.stack((u, v), axis=1)
        polygon = ConvexHull(projected[on])
        full_area = polygon.volume  # SciPy's 2-D hull volume is polygon area.
        physical = (np.max(abs(triangles @ normal+offset), axis=1) <= eps)
        physical &= (normals @ normal >= 1-8*eps)
        occupied = float(areas[physical].sum())
        if np.isclose(occupied, full_area, rtol=1e-9, atol=eps**2):
            continue
        if occupied > eps**2:
            raise UnsupportedClosureError(
                'physical hull has a nonplanar opening cutting a partially occupied convex facet; '
                'the internal closure backend needs constrained polygon difference for this valid hull')
        corners = on[polygon.vertices]
        loop = []
        for a, b in zip(corners, np.roll(corners, -1)):
            edge = points[b]-points[a]
            t = (points[on]-points[a]) @ edge/(edge @ edge)
            distance = np.linalg.norm(points[on]-points[a]-t[:, None]*edge, axis=1)
            selected = np.flatnonzero((distance <= eps) & (t >= -eps) & (t < 1-eps))
            loop.extend(on[selected[np.argsort(t[selected])]].tolist())
        if set(loop) != set(on):
            raise UnsupportedClosureError('internal convex-facet closure cannot preserve an interior input vertex')
        cap.extend(_ears(loop, points, normal))
    if not cap:
        raise UnsupportedClosureError('internal convex-envelope closure produced no missing facets')
    return np.asarray(cap, dtype=np.int64)


def _cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _inside(point, triangle, eps):
    signs = np.array([_cross2(b-a, point-a)
                      for a, b in zip(triangle, np.roll(triangle, -1, axis=0))])
    return np.all(signs >= -eps) or np.all(signs <= eps)


def _coplanar_contacts(a, b, normal, eps):
    """Intersection witnesses, including containment and collinear edges."""
    drop = int(np.argmax(abs(normal)))
    aa, bb = np.delete(a, drop, axis=1), np.delete(b, drop, axis=1)
    contacts = [p for p, q in zip(a, aa) if _inside(q, bb, eps)]
    contacts.extend(p for p, q in zip(b, bb) if _inside(q, aa, eps))
    for i in range(3):
        p, q = aa[i], aa[(i+1) % 3]
        u = q-p
        for j in range(3):
            r, s = bb[j], bb[(j+1) % 3]
            v = s-r
            det = _cross2(u, v)
            if abs(det) <= eps:
                continue  # Overlap endpoints are included by the containment tests.
            t, w = _cross2(r-p, v)/det, _cross2(r-p, u)/det
            if -eps <= t <= 1+eps and -eps <= w <= 1+eps:
                contacts.append(a[i] + np.clip(t, 0., 1.)*(a[(i+1) % 3]-a[i]))
    return contacts


def _triangle_contacts(a, b, na, nb, eps):
    if np.linalg.norm(np.cross(na, nb)) <= eps:
        if np.max(abs((a-b[0]) @ nb)) > eps:
            return []
        return _coplanar_contacts(a, b, na, eps)
    contacts = []
    for source, target, normal in ((a, b, nb), (b, a, na)):
        drop = int(np.argmax(abs(normal)))
        projected = np.delete(target, drop, axis=1)
        distances = (source-target[0]) @ normal
        for i in range(3):
            j = (i+1) % 3
            if abs(distances[i]) <= eps and _inside(np.delete(source[i], drop), projected, eps):
                contacts.append(source[i])
            if distances[i] * distances[j] < 0:
                fraction = distances[i]/(distances[i]-distances[j])
                point = source[i] + fraction*(source[j]-source[i])
                if _inside(np.delete(point, drop), projected, eps):
                    contacts.append(point)
    return contacts


def _check_intersections(points, faces, normals, eps):
    triangles = points[faces]
    lower, upper = triangles.min(axis=1), triangles.max(axis=1)
    for i, a in enumerate(triangles):
        candidates = np.flatnonzero(np.all(upper[i] >= lower[i+1:]-eps, axis=1)
                                    & np.all(upper[i+1:] >= lower[i]-eps, axis=1)) + i+1
        for j in candidates:
            common = np.intersect1d(faces[i], faces[j])
            for point in _triangle_contacts(a, triangles[j], normals[i], normals[j], eps):
                if len(common) == 1:
                    allowed = np.linalg.norm(point-points[common[0]]) <= 8*eps
                elif len(common) == 2:
                    start, end = points[common]
                    edge = end-start
                    t = np.clip((point-start) @ edge/(edge @ edge), 0., 1.)
                    allowed = np.linalg.norm(point-start-t*edge) <= 8*eps
                else:
                    allowed = False
                if not allowed:
                    raise ValueError('physical hull has intersecting faces or a geometric T-junction')


def physical_hull(vertices_m, faces):
    """Return an internal HullMesh from one uncapped convex physical hull.

    Rejects closed/capped input, multiple openings, self intersections, vertex
    pinches, inconsistent winding and nonconvex geometry. Planar openings and
    nonplanar openings comprising whole missing convex facets can be closed.
    Other valid openings raise UnsupportedClosureError, a backend limitation.
    Collinear boundary vertices are preserved. Validation
    uses coordinates normalized by the hull span, with 1e-10 geometric tolerance.
    No welding, repair, normal reversal or physical panel replacement occurs.
    """
    vertices = array(vertices_m, (None, 3), 'vertices_m')
    raw = np.asarray(faces)
    if raw.ndim != 2 or raw.shape[1:] != (3,) or raw.dtype.kind not in 'iu':
        raise ValueError('faces must be an integer (triangle, 3) array')
    if len(vertices) < 4 or len(raw) < 3 or np.any(raw < 0) or np.any(raw >= len(vertices)):
        raise ValueError('physical hull needs at least four vertices and valid triangle indices')
    faces = np.array(raw, dtype=np.int64, order='C', copy=True)
    if len(np.unique(vertices, axis=0)) != len(vertices):
        raise ValueError('vertices must be unique; no implicit welding is performed')
    if len(np.unique(faces)) != len(vertices):
        raise ValueError('physical hull contains unused vertices')
    if len(np.unique(np.sort(faces, axis=1), axis=0)) != len(faces):
        raise ValueError('physical hull contains duplicate faces')
    span = float(np.max(np.ptp(vertices, axis=0)))
    if span <= 0:
        raise ValueError('physical hull must enclose a three-dimensional volume')
    points = (vertices-vertices[0])/span
    triangles = points[faces]
    normal = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    lengths = np.linalg.norm(normal, axis=1)
    eps = 1e-10
    if np.any(lengths <= 64*np.finfo(float).eps):
        raise ValueError('physical hull contains degenerate triangles')
    normal /= lengths[:, None]
    edges = {}
    incident = [[] for _ in vertices]
    for i, face in enumerate(faces):
        for v in face:
            incident[v].append(i)
        for u, v in zip(face, np.roll(face, -1)):
            edges.setdefault(tuple(sorted((int(u), int(v)))), []).append((i, int(u), int(v)))
    boundary, adjacency = [], [set() for _ in faces]
    for owners in edges.values():
        if len(owners) == 1:
            boundary.append(owners[0][1:])
        elif len(owners) == 2:
            a, b = owners
            if a[1:] != b[1:][::-1]:
                raise ValueError('physical hull faces have inconsistent winding')
            adjacency[a[0]].add(b[0]); adjacency[b[0]].add(a[0])
        else:
            raise ValueError('physical hull has non-manifold edges')
    if not boundary:
        raise ValueError('boat must have one uncapped opening; closed or capped hull supplied')
    reached, stack = set(), [0]
    while stack:
        i = stack.pop()
        if i not in reached:
            reached.add(i); stack.extend(adjacency[i]-reached)
    if len(reached) != len(faces):
        raise ValueError('physical hull must be one connected face shell')
    next_vertex = {}
    incoming = {}
    for u, v in boundary:
        if u in next_vertex or v in incoming:
            raise ValueError('opening boundary has a branch or pinched vertex')
        next_vertex[u] = v; incoming[v] = u
    if set(next_vertex) != set(incoming):
        raise ValueError('opening boundary must form a closed loop')
    loop, start = [], min(next_vertex)
    current = start
    while current not in loop:
        loop.append(current); current = next_vertex[current]
    if current != start or len(loop) != len(boundary):
        raise ValueError('boat must have exactly one opening loop')
    if len(vertices)-len(edges)+len(faces) != 1:
        raise ValueError('physical hull must be a topological disk')
    # At each vertex its face link must be one path (rim) or one cycle (interior).
    for v, owners in enumerate(incident):
        link = {}
        for i in owners:
            a, b = [int(x) for x in faces[i] if x != v]
            link.setdefault(a, set()).add(b); link.setdefault(b, set()).add(a)
        degrees = [len(neighbors) for neighbors in link.values()]
        expected_ends = 2 if v in next_vertex else 0
        if any(d not in (1, 2) for d in degrees) or degrees.count(1) != expected_ends:
            raise ValueError('physical hull has a non-manifold vertex')
        seen, todo = set(), [next(iter(link))]
        while todo:
            node = todo.pop()
            if node not in seen:
                seen.add(node); todo.extend(link[node]-seen)
        if len(seen) != len(link):
            raise ValueError('physical hull has a pinched vertex')
    # Convexity is tested against every supplied physical support plane.
    for triangle, n in zip(triangles, normal):
        if np.max((points-triangle[0]) @ n) > eps:
            raise ValueError('outward convex physical hull required')
    _check_intersections(points, faces, normal, eps)
    # The cap follows the opposite direction to the physical boundary.
    rim = np.array(loop[::-1], dtype=np.int64)
    center = points[rim].mean(axis=0)
    _, _, basis = np.linalg.svd(points[rim]-center, full_matrices=False)
    if np.max(abs((points[rim]-center) @ basis[-1])) > eps:
        cap = _envelope_cap(points, faces, normal, eps)
        closed_faces = np.vstack((faces, cap))
        return HullMesh(vertices, closed_faces, opening_rim=rim,
                        artificial_faces=np.arange(len(faces), len(closed_faces)))
    area_vector = np.cross(points[rim]-center, points[np.roll(rim, -1)]-center).sum(axis=0)
    if np.linalg.norm(area_vector) <= eps:
        raise ValueError('opening must enclose a planar convex area')
    cap_normal = area_vector/np.linalg.norm(area_vector)
    if np.max((points-center) @ cap_normal) > eps:
        raise ValueError('opening must bound the convex hull')
    for a, b in zip(points[rim], points[np.roll(rim, -1)]):
        if np.min(np.cross(b-a, points[rim]-a) @ cap_normal) < -eps:
            raise ValueError('opening boundary must be a convex simple polygon')
    cap_vertex = len(vertices)
    cap = np.column_stack((np.full(len(rim), cap_vertex), rim, np.roll(rim, -1)))
    closed_vertices = np.vstack((vertices, vertices[rim].mean(axis=0)))
    closed_faces = np.vstack((faces, cap))
    return HullMesh(closed_vertices, closed_faces, opening_rim=rim,
                    artificial_faces=np.arange(len(faces), len(closed_faces)))
