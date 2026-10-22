"""Contained affine tetra-P4 accuracy probe, not a production PDE backend.

Geometry stays affine. Potential has the complete degree-four polynomial
space on each tetrahedron; neighboring tetrahedra share all trace nodes.
"""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from numpy.polynomial import Polynomial, legendre
from scipy.sparse import coo_matrix


@lru_cache(None)
def indices(dimension):
    return np.array([a for a in np.ndindex(*([5]*dimension)) if sum(a) == 4])


def basis(barycentric):
    """P4 nodal values and barycentric partial derivatives, triangles or tets."""
    barycentric = np.asarray(barycentric)
    alpha = indices(barycentric.shape[1])
    polynomials = []
    for count in range(5):
        polynomial = Polynomial([1.])
        for j in range(count):
            polynomial *= Polynomial([-j, 4])/(count-j)
        polynomials.append(polynomial)
    values = np.ones((len(barycentric), len(alpha)))
    derivatives = np.zeros((*values.shape, barycentric.shape[1]))
    for a, counts in enumerate(alpha):
        factors = np.array([polynomials[count](barycentric[:, i])
                            for i, count in enumerate(counts)]).T
        values[:, a] = np.prod(factors, axis=1)
        for i, count in enumerate(counts):
            derivatives[:, a, i] = (
                polynomials[count].deriv()(barycentric[:, i])
                * np.prod(np.delete(factors, i, axis=1), axis=1))
    return values, derivatives


@lru_cache(None)
def reference_gradient():
    # Degree-six gradient products times the Duffy Jacobian have maximum
    # coordinate degree eight, exactly integrated by five-point Gauss.
    x, w = legendre.leggauss(5)
    x, w = (x+1)/2, w/2
    u, v, t = np.meshgrid(x, x, x, indexing='ij')
    wu, wv, wt = np.meshgrid(w, w, w, indexing='ij')
    xyz = np.column_stack((u.ravel(), ((1-u)*v).ravel(),
                           ((1-u)*(1-v)*t).ravel()))
    barycentric = np.column_stack((1-xyz.sum(axis=1), xyz))
    _, derivative = basis(barycentric)
    return (derivative[:, :, 1:]-derivative[:, :, :1],
            (wu*wv*wt*(1-u)**2*(1-v)).ravel())


def key(vertices, alpha):
    return tuple(sorted((int(vertex), int(count))
                        for vertex, count in zip(vertices, alpha) if count))


@dataclass(frozen=True)
class Mesh:
    points: np.ndarray
    tets: np.ndarray
    elements: np.ndarray
    node_by_key: dict
    vertex_count: int


def make_mesh(points, tets):
    """Share nodes by topological barycentric weights, without coordinate hashing."""
    points, tets = np.asarray(points), np.asarray(tets)
    coordinates = list(points)
    lookup = {((i, 4),): i for i in range(len(points))}
    elements = np.empty((len(tets), 35), dtype=np.int32)
    for t, vertices in enumerate(tets):
        for a, alpha in enumerate(indices(4)):
            identity = key(vertices, alpha)
            if identity not in lookup:
                lookup[identity] = len(coordinates)
                coordinates.append(alpha@points[vertices]/4)
            elements[t, a] = lookup[identity]
    return Mesh(np.asarray(coordinates), tets, elements, lookup, len(points))


def face_nodes(mesh, faces):
    return np.array([[mesh.node_by_key[key(vertices, alpha)]
                      for alpha in indices(3)] for vertices in faces], dtype=np.int32)


def stiffness(mesh):
    """Exact physical integral of grad(P4) dot grad(P4), affine geometry."""
    gradient, weights = reference_gradient()
    # The brick experiment has only six element shapes. Cache exact geometry,
    # not coordinates rounded into a tolerance-dependent class.
    local_by_edges = {}
    matrices = []
    for vertices in mesh.tets:
        coordinates = mesh.points[vertices]
        edges = coordinates[1:]-coordinates[:1]
        identity = edges.tobytes()
        if identity not in local_by_edges:
            physical = gradient@np.linalg.inv(edges).T
            local_by_edges[identity] = np.einsum(
                'qai,qbi,q->ab', physical, physical, weights)*abs(np.linalg.det(edges))
        matrices.append(local_by_edges[identity])
    nodes = mesh.elements
    return coo_matrix((np.asarray(matrices).ravel(),
                       (np.repeat(nodes, 35, axis=1).ravel(),
                        np.tile(nodes, (1, 35)).ravel())),
                      shape=(len(mesh.points), len(mesh.points))).tocsr()


def restrict_trace(mesh, interface_faces):
    """R maps retained coefficients to a trace exactly P1 on every interface face.

    Every higher-order interface node, including shared edges, is eliminated
    by its barycentric interpolation of the original triangle vertices.
    Retained columns are indexed by the returned full-node IDs.
    """
    constraints = {}
    for face in interface_faces:
        for alpha in indices(3):
            identity = key(face, alpha)
            node = mesh.node_by_key[identity]
            if node < mesh.vertex_count:
                continue
            constraint = tuple((int(vertex), float(count)/4)
                               for vertex, count in identity)
            if node in constraints and constraints[node] != constraint:
                raise ValueError('inconsistent shared-edge P1 trace restriction')
            constraints[node] = constraint
    retained = np.setdiff1d(np.arange(len(mesh.points)), list(constraints))
    column = np.full(len(mesh.points), -1, dtype=int)
    column[retained] = np.arange(len(retained))
    rows, cols, values = list(retained), list(column[retained]), [1.]*len(retained)
    for node, constraint in constraints.items():
        for vertex, weight in constraint:
            rows.append(node); cols.append(column[vertex]); values.append(weight)
    restriction = coo_matrix((values, (rows, cols)),
                             shape=(len(mesh.points), len(retained))).tocsr()
    return restriction, retained
