"""Fixed reference extensions of boundary displacement into the P4 fluid mesh.

This is an ALE coordinate choice, not a fluid equation or an operator-update
approximation. The fixed reference map supplies exact tangent/transpose actions.
It does not guarantee valid deformed elements; physical geometry guards remain
the caller's responsibility.
"""
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from .p4_tetra import indices, face_nodes
from .pocket import boundary


@dataclass(frozen=True)
class VolumeMap:
    reference: np.ndarray
    boundary: np.ndarray
    interior: np.ndarray
    coupling: object
    factor: object

    def lift(self, boundary_displacement):
        """Map boundary displacement/velocity to all allocated geometry nodes."""
        value = np.asarray(boundary_displacement, dtype=float)
        if value.shape != (len(self.boundary), 3) or not np.all(np.isfinite(value)):
            raise ValueError('finite boundary vectors (boundary_nodes,3) required')
        result = np.zeros_like(self.reference)
        result[self.boundary] = value
        if len(self.interior):
            result[self.interior] = self.factor.solve(-self.coupling@value)
        return result

    def pullback(self, nodal_force):
        """Exact transpose lift: return the force on boundary coordinates."""
        value = np.asarray(nodal_force, dtype=float)
        if value.shape != self.reference.shape or not np.all(np.isfinite(value)):
            raise ValueError('finite force vectors matching allocated geometry nodes required')
        result = value[self.boundary].copy()
        if len(self.interior):
            result -= self.coupling.T@self.factor.solve(value[self.interior], trans='T')
        return result


def make_map(mesh):
    """Positive reference springs between adjacent P4 barycentric lattice nodes.

    Physical boundary nodes are prescribed. Unused input vertices have no
    geometry displacement or force path. No dense volume-to-surface Jacobian
    is formed. The factor belongs only to this fixed reference geometry map.
    """
    alpha = indices(4)
    adjacency = np.sum(abs(alpha[:, None]-alpha[None, :]), axis=2) == 2
    local_edges = np.column_stack(np.where(np.triu(adjacency, 1)))
    edges = np.unique(np.sort(mesh.elements[:, local_edges].reshape(-1, 2), axis=1), axis=0)
    squared = np.sum(np.diff(mesh.points[edges], axis=1)[:, 0]**2, axis=1)
    if not np.all(np.isfinite(squared)) or np.any(squared <= 0):
        raise ValueError('reference P4 graph requires distinct finite connected points')
    weight = 1/squared
    a, b = edges.T
    laplacian = coo_matrix((np.r_[weight, weight, -weight, -weight],
        (np.r_[a, b, a, b], np.r_[a, b, b, a])),
        shape=(len(mesh.points), len(mesh.points))).tocsr()
    prescribed = np.unique(face_nodes(mesh, boundary(mesh.points, mesh.tets)))
    active = np.unique(mesh.elements)
    interior = np.setdiff1d(active, prescribed)
    coupling = laplacian[interior][:, prescribed].tocsr()
    factor = splu(laplacian[interior][:, interior].tocsc()) if len(interior) else None
    return VolumeMap(mesh.points.copy(), prescribed, interior, coupling, factor)


@dataclass(frozen=True)
class HierarchicalVolumeMap:
    """P1 FEM corner extension plus the genuine P4 boundary residual.

    If the prescribed boundary displacement is a P1 trace, the result is
    piecewise affine in every tetrahedron. Independent higher-order boundary
    modes remain independent and are extended by the original fixed graph.
    Neither extension guarantees positive current Jacobians.
    """
    reference: np.ndarray
    boundary: np.ndarray
    corners: VolumeMap
    graph: VolumeMap
    corner_boundary_rows: np.ndarray
    barycentric_lift: object

    def lift(self, boundary_displacement):
        value = np.asarray(boundary_displacement, dtype=float)
        if value.shape != (len(self.boundary), 3) or not np.all(np.isfinite(value)):
            raise ValueError('finite boundary vectors (boundary_nodes,3) required')
        base = self.barycentric_lift@self.corners.lift(value[self.corner_boundary_rows])
        return base+self.graph.lift(value-base[self.boundary])

    def pullback(self, nodal_force):
        value = np.asarray(nodal_force, dtype=float)
        if value.shape != self.reference.shape or not np.all(np.isfinite(value)):
            raise ValueError('finite force vectors matching allocated geometry nodes required')
        residual_force = self.graph.pullback(value)
        base_force = value.copy()
        base_force[self.boundary] -= residual_force
        corner_force = self.corners.pullback(self.barycentric_lift.T@base_force)
        result = residual_force.copy()
        result[self.corner_boundary_rows] += corner_force
        return result


def make_hierarchical_map(mesh):
    """Two fixed reference factors; no dense volume-to-boundary derivative.

    The reference FEM metric accounts for thin-cell altitudes, unlike edge
    springs alone. Its harmonic corner values need not be convex weights.
    The P4 graph acts only on the remaining high-order boundary residual.
    """
    from .pocket import stiffness

    graph = make_map(mesh)
    vertex_points = mesh.points[:mesh.vertex_count]
    prescribed = np.unique(boundary(vertex_points, mesh.tets))
    active = np.unique(mesh.tets)
    interior = np.setdiff1d(active, prescribed)
    operator = stiffness(vertex_points, mesh.tets)
    coupling = operator[interior][:, prescribed].tocsr()
    factor = splu(operator[interior][:, interior].tocsc()) if len(interior) else None
    corners = VolumeMap(vertex_points.copy(), prescribed, interior, coupling, factor)
    rows = np.searchsorted(graph.boundary, prescribed)
    if not np.array_equal(graph.boundary[rows], prescribed):
        raise RuntimeError('physical corner boundary is absent from the P4 boundary')
    row, column, value = [], [], []
    for identity, node in mesh.node_by_key.items():
        for vertex, numerator in identity:
            row.append(node); column.append(vertex); value.append(numerator/4)
    lift = coo_matrix((value, (row, column)),
                      shape=(len(mesh.points), mesh.vertex_count)).tocsr()
    return HierarchicalVolumeMap(mesh.points.copy(), graph.boundary, corners,
                                 graph, rows, lift)
