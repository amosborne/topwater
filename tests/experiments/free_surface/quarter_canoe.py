"""Canonical quarter PLC and exact reflected tetrahedral discretization.

Static rest only. Neither mesh symmetry nor a parity solve validates a moving
nonlinear boundary or outgoing-wave accuracy. Meshing is explicit, never done
by importing this module. No coordinate rounding or node welding is used.
"""
from dataclasses import dataclass
from itertools import permutations, product

import numpy as np
from scipy.spatial import ConvexHull, Delaunay

from topwater._clipping import clip_flat_water
from topwater.frames import rotation_matrix
from tests.experiments.capytaine_mesh_cost import subdivide
from .moving_canoe import make_domain
from .pocket import boundary
from .p4_canoe import tetrahedralize, triangle_keys
from . import p4_tetra


POCKET = np.array([2*np.pi/3, np.pi/3])
PARITY = np.array([[-1, 1], [1, -1], [1, 1], [1, -1], [-1, 1], [-1, -1]])


@dataclass(frozen=True)
class QuarterPLC:
    points: np.ndarray
    faces: np.ndarray
    markers: np.ndarray
    body_parent: np.ndarray
    exterior_points: np.ndarray
    exterior_tets: np.ndarray
    center: np.ndarray
    parent_reflections: np.ndarray


def _positive(points, tets):
    tets = np.array(tets, dtype=int, copy=True)
    determinant = np.linalg.det(points[tets[:, 1:]]-points[tets[:, :1]])
    if not np.all(np.isfinite(determinant)) or np.any(determinant == 0):
        raise ValueError('degenerate tetrahedral geometry')
    negative = determinant < 0
    tets[negative, :2] = tets[negative, 1::-1]
    return tets


def _notched_plane(grid, curve, axis, level, normal):
    """Triangulate this canoe's convex corner notch; reject lost constraints.

    The body cross-section/footprint occupies the corner next to the origin.
    Convexity and every prescribed segment are checked; this is deliberately
    not a general constrained-triangulation implementation.
    """
    axes = [a for a in range(3) if a != axis]
    curve = np.unique(curve, axis=0)
    origin = np.zeros(3); origin[axis] = level
    polygon = np.vstack((curve[:, axes], origin[axes]))
    hull = ConvexHull(polygon)
    signed = grid[:, axes]@hull.equations[:, :2].T+hull.equations[:, 2]
    outside = np.any(signed > 64*np.finfo(float).eps, axis=1)
    points = np.unique(np.vstack((grid[outside], curve)), axis=0)
    planar = points[:, axes]
    triangles = Delaunay(planar).simplices
    center = planar[triangles].mean(axis=1)
    keep = np.any(center@hull.equations[:, :2].T+hull.equations[:, 2]
                  > 64*np.finfo(float).eps, axis=1)
    triangles = triangles[keep]
    p = points[triangles]
    reverse = np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0])@normal < 0
    triangles[reverse] = triangles[reverse, ::-1]
    return points[triangles]


def make_plc(x, y, z, *, subdivisions=1):
    """Build quarter near-boundary and untouched structured exterior.

    x/y start at zero and include the exact pocket planes; z includes -4,-1,0.
    Markers: body1, free surface2, buried interface3, symmetry x4/y5.
    """
    x, y, z = [np.asarray(a, float) for a in (x, y, z)]
    for values in (x, y, z):
        if values.ndim != 1 or len(values) < 2 or not np.all(np.isfinite(values)) or np.any(np.diff(values) <= 0):
            raise ValueError('strictly increasing finite coordinate vectors required')
    if x[0] != 0 or y[0] != 0 or z[0] != -4 or z[-1] != 0 or -1 not in z or POCKET[0] not in x or POCKET[1] not in y:
        raise ValueError('retain exact pocket planes, quarter axes, depth4 and flat surface')
    if isinstance(subdivisions, bool) or subdivisions not in (0, 1, 2):
        raise ValueError('use numerical hull subdivision 0, 1 or 2')
    domain = make_domain(12, False, 0)
    body, rest = domain.body, domain.rest
    world = (body.mesh.vertices-body.cg)@rotation_matrix(rest.orientation).T+rest.position
    physical_ids = np.asarray(body.mesh.physical_faces, dtype=int)
    physical_lookup = {tuple(sorted(map(tuple, world[body.mesh.faces[i]]))): int(i)
                       for i in physical_ids}
    parent_reflections = np.full((4, len(body.mesh.faces)), -1, int)
    for row, (sx, sy) in enumerate(product((1, -1), repeat=2)):
        for i in physical_ids:
            reflected = world[body.mesh.faces[i]]*np.array([sx, sy, 1])
            identity = tuple(sorted(map(tuple, reflected)))
            if identity not in physical_lookup:
                raise ValueError('original physical facets lack exact reflection symmetry')
            parent_reflections[row, i] = physical_lookup[identity]
    wet, parents, *_ = clip_flat_water(world, body.mesh.faces, body.mesh.physical_faces, 0., True)
    for _ in range(subdivisions):
        wet = subdivide(wet); parents = np.tile(parents, 4)
    for axis in (0, 1):
        if np.any((wet[:, :, axis].min(axis=1) < 0)&(wet[:, :, axis].max(axis=1) > 0)):
            raise ValueError('physical wet facet crosses a proposed symmetry plane')
    selected = np.all(wet[:, :, :2] >= 0, axis=(1, 2))
    wet, parents = wet[selected], parents[selected]
    if not len(wet):
        raise ValueError('empty canonical quarter hull')
    # Body outward -> fluid outward. Preserve parent face IDs in this order.
    wet = wet[:, ::-1]
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    numbering = np.arange(len(points)).reshape(xx.shape)
    steps = np.zeros((6, 4, 3), int)
    steps[:, 1:] = np.cumsum(np.eye(3, dtype=int)[list(permutations(range(3)))], axis=1)
    offsets = steps@np.array([len(y)*len(z), len(z), 1])
    tets = (numbering[:-1, :-1, :-1].ravel()[:, None, None]+offsets).reshape(-1, 4)
    tets = _positive(points, tets)
    near = np.all(points[tets, 0] <= POCKET[0], axis=1)&np.all(points[tets, 1] <= POCKET[1], axis=1)&np.all(points[tets, 2] >= -1, axis=1)
    near_faces = boundary(points, tets[near])
    triangle = points[near_faces]
    top = np.all(triangle[:, :, 2] == 0, axis=1)
    sx = np.all(triangle[:, :, 0] == 0, axis=1)
    sy = np.all(triangle[:, :, 1] == 0, axis=1)
    pieces = [wet, triangle[~(top|sx|sy)]]
    labels = [np.ones(len(wet), int), np.full(np.sum(~(top|sx|sy)), 3)]
    for axis, level, mask, marker, normal in (
            (2, 0., top, 2, [0, 0, 1]), (0, 0., sx, 4, [-1, 0, 0]),
            (1, 0., sy, 5, [0, -1, 0])):
        curve = wet.reshape(-1, 3); curve = curve[curve[:, axis] == level]
        planar = _notched_plane(points[np.unique(near_faces[mask])], curve, axis, level, np.array(normal))
        pieces.append(planar); labels.append(np.full(len(planar), marker))
    triangles = np.concatenate(pieces)
    vertices, inverse = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
    faces = inverse.reshape(-1, 3); markers = np.concatenate(labels)
    edges = np.sort(np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]])), axis=1)
    if np.any(np.unique(edges, axis=0, return_counts=True)[1] != 2):
        raise ValueError('quarter PLC lost a curve/interface edge; no automatic repair')
    # Winding as well as incidence: every oriented edge must cancel exactly.
    directed = np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    if set(map(tuple, directed)) != set(map(tuple, directed[:, ::-1])):
        raise ValueError('quarter PLC boundary winding is inconsistent')
    return QuarterPLC(vertices, faces, markers, parents, points, tets[~near], rest.position,
                      parent_reflections)


def make_quarter(plc):
    """Explicit optional TetGen preprocessing, with exact interface checks."""
    local, tets, faces, markers, report = tetrahedralize(plc.points, plc.faces, plc.markers)
    for marker in (4, 5):
        if triangle_keys(local[faces[markers == marker]]) != triangle_keys(plc.points[plc.faces[plc.markers == marker]]):
            raise ValueError('mesher changed a symmetry-plane facet')
    if not np.array_equal(local[:len(plc.points)], plc.points):
        raise ValueError('mesher changed the exact input vertex prefix')
    points = list(plc.exterior_points); lookup = {tuple(p): i for i, p in enumerate(points)}
    mapping = []
    for p in local:
        identity = tuple(p)
        if identity not in lookup:
            lookup[identity] = len(points); points.append(p)
        mapping.append(lookup[identity])
    points = np.asarray(points); mapping = np.asarray(mapping)
    interface = mapping[faces[markers == 3]]
    key = lambda f: {tuple(sorted(map(int, row))) for row in f}
    if not key(interface).issubset(key(boundary(points, plc.exterior_tets))):
        raise ValueError('quarter buried triangular interface changed')
    tets = np.vstack((plc.exterior_tets, mapping[tets]))
    # Prune all inactive original grid vertices before constructing P4 nodes.
    active, inverse = np.unique(tets, return_inverse=True)
    remap = np.full(len(points), -1, int); remap[active] = np.arange(len(active))
    points = points[active]; tets = _positive(points, inverse.reshape(-1, 4))
    body = remap[mapping[faces[markers == 1]]]
    parent_by_triangle = {tuple(sorted(map(tuple, tri))): int(parent) for tri, parent in
                          zip(plc.points[plc.faces[plc.markers == 1]], plc.body_parent)}
    parents = np.array([parent_by_triangle[tuple(sorted(map(tuple, tri)))] for tri in points[body]])
    mesh = p4_tetra.make_mesh(points, tets)
    physical = boundary(points, tets)
    return mesh, physical, body, parents, report


def mirror(quarter):
    """Return exact four-copy mesh and P4 maps, including quadrature swaps.

    An odd reflection swaps corners0/1 for positive orientation. Its integration
    rule must swap those same barycentric entries to mirror the discrete rule.
    """
    points = []; lookup = {}; vertex_maps = []; tets = []; swapped = []
    signs = list(product((1, -1), repeat=2))
    for sx, sy in signs:
        mapping = []
        for vertex, p in enumerate(quarter.points[:quarter.vertex_count]):
            identity = (vertex, sx if p[0] != 0 else 0, sy if p[1] != 0 else 0)
            if identity not in lookup:
                lookup[identity] = len(points); points.append(p*np.array([sx, sy, 1]))
            mapping.append(lookup[identity])
        mapping = np.asarray(mapping); vertex_maps.append(mapping)
        reflected = mapping[quarter.tets].copy(); odd = sx*sy < 0
        if odd:
            reflected[:, :2] = reflected[:, 1::-1]
        tets.append(reflected); swapped.extend([odd]*len(reflected))
    full = p4_tetra.make_mesh(np.asarray(points), np.vstack(tets))
    maps = []
    for mapping in vertex_maps:
        nodes = np.empty(len(quarter.points), int)
        for identity, node in quarter.node_by_key.items():
            transformed = tuple(sorted((int(mapping[v]), count) for v, count in identity))
            nodes[node] = full.node_by_key[transformed]
        maps.append(nodes)
    if np.any(np.linalg.det(full.points[full.tets[:, 1:]]-full.points[full.tets[:, :1]]) <= 0):
        raise ValueError('reflected tetrahedron lost positive orientation')
    return full, np.asarray(maps), np.asarray(swapped), np.asarray(signs)


def mirrored_body(full, maps, quarter_body, quarter_parents, parent_reflections):
    """Fresh fluid-outward boundary faces and reflected original parent IDs."""
    parent_by_key = {}
    for row, mapping in enumerate(maps):
        for face, parent in zip(mapping[quarter_body], quarter_parents):
            identity = tuple(sorted(map(int, face)))
            if identity in parent_by_key:
                raise ValueError('body face lies in a symmetry plane or is duplicated')
            transformed = parent_reflections[row, parent]
            if transformed < 0:
                raise ValueError('missing reflected original physical parent')
            parent_by_key[identity] = transformed
    body = []; parents = []
    for face in boundary(full.points, full.tets):
        identity = tuple(sorted(map(int, face)))
        if identity in parent_by_key:
            body.append(face); parents.append(parent_by_key.pop(identity))
    if parent_by_key:
        raise ValueError('mirrored body does not belong to the physical boundary')
    return np.asarray(body), np.asarray(parents)
