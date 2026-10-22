"""Reference graph extension of P4 boundary displacement into the fluid mesh.

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
