"""Physical-exterior P1 mass and algebraic hidden fill, CPU reference only."""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from .compact_p1_trace import CompactP1Trace


def _clip(polygon, plane, inside):
    """Strict float integration cut carrying original barycentric values."""
    if not len(polygon):
        return np.empty((0, 5))
    values = polygon[:, :2] @ plane[:2]+plane[2]
    if not inside:
        values = -values
    result = []
    for a, b, va, vb in zip(polygon, np.roll(polygon, -1, axis=0), values, np.roll(values, -1)):
        if va <= 0:
            result.append(a)
        if (va < 0 < vb) or (vb < 0 < va):
            result.append(a+(b-a)*(va/(va-vb)))
    return np.asarray(result).reshape(-1, 5)


def _signed_area(points):
    if len(points) < 3:
        return 0.
    r = points[:, :2]-points[0, :2]
    return float(np.sum(r[:, 0]*np.roll(r[:, 1], -1)-r[:, 1]*np.roll(r[:, 0], -1))/2)


class PhysicalP1Patch:
    """Structured BLTR patch; only positive physical-mass interior nodes evolve.

    ``E`` expands active coefficients to ALL row-major nodes. Perimeter nodes
    are fixed zero; exact-zero-support interior nodes are algebraic, using a
    five-point harmonic fill when requested. Fill never enters physical mass.
    ``condense_maps`` accepts maps with all nodal IDs in their LAST dimension.
    Float cut geometry is a verification reference, not a topology certificate.
    """
    def __init__(self, x, y, footprint_ccw, *, harmonic_fill=True):
        source = CompactP1Trace(x, y, np.zeros((len(y), len(x))))
        self.x, self.y, self.panels, self.faces = source.x, source.y, source.panels, source.faces
        nx, ny = len(self.x), len(self.y)
        xx, yy = np.meshgrid(self.x, self.y)
        self.nodes = np.column_stack((xx.ravel(), yy.ravel()))
        footprint = np.asarray(footprint_ccw, dtype=float)
        if (footprint.ndim != 2 or footprint.shape[1] != 2 or len(footprint) < 3
                or not np.isfinite(footprint).all()):
            raise ValueError('finite convex CCW footprint required')
        edge = np.roll(footprint, -1, axis=0)-footprint
        following = np.roll(edge, -1, axis=0)
        if np.any(edge[:, 0]*following[:, 1]-edge[:, 1]*following[:, 0] <= 0):
            raise ValueError('strictly convex CCW footprint required')
        normals = np.column_stack((edge[:, 1], -edge[:, 0]))
        planes = np.column_stack((normals, -np.sum(normals*footprint, axis=1)))
        # Convex support, not merely positive local turns (star polygons fail).
        support = footprint @ normals.T+planes[:, 2]
        validation_roundoff = 64*np.finfo(float).eps*max(1., np.max(np.abs(footprint)))*np.linalg.norm(normals, axis=1)
        if np.any(support > validation_roundoff):
            raise ValueError('footprint is not convex')
        cuts, parents, hidden_cuts, hidden_parents, degeneracies = [], [], [], [], []

        def retain(polygon, parent, hidden=False):
            if not len(polygon):
                return
            area = _signed_area(polygon)
            if area < 0:
                raise ValueError('cut polygon orientation reversed')
            if area == 0:
                degeneracies.append((parent, hidden, polygon.copy()))
                return
            center = polygon.mean(axis=0)
            for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
                triangle = np.array([center, a, b])
                if _signed_area(triangle) <= 0:
                    raise ValueError('degenerate positive cut fan; no vertex deletion')
                (hidden_cuts if hidden else cuts).append(triangle)
                (hidden_parents if hidden else parents).append(parent)

        for parent, triangle in enumerate(self.panels.triangles[:, :, :2]):
            polygon = np.column_stack((triangle, np.eye(3)))
            side = triangle @ planes[:, :2].T+planes[:, 2]
            if np.any(np.all(side >= 0, axis=0)):
                # Entire triangle exterior to at least one convex half-plane.
                cuts.append(polygon)
                parents.append(parent)
            elif np.all(side <= 0):
                hidden_cuts.append(polygon)
                hidden_parents.append(parent)
                continue
            else:
                remaining = polygon
                for plane in planes:
                    retain(_clip(remaining, plane, False), parent)
                    remaining = _clip(remaining, plane, True)
                retain(remaining, parent, hidden=True)
        cut = np.asarray(cuts).reshape(-1, 3, 5)
        self.cut_triangles = cut[:, :, :2]
        self.basis_values = cut[:, :, 2:]
        self.parent_id = np.asarray(parents, dtype=int)
        self.cut_areas = np.array([_signed_area(t) for t in cut])
        hidden_cut = np.asarray(hidden_cuts).reshape(-1, 3, 5)
        # Hidden pieces are integration-subtraction geometry ONLY, never mass
        # or evolved water. They retain the same original-parent basis.
        self.hidden_triangles = hidden_cut[:, :, :2]
        self.hidden_basis_values = hidden_cut[:, :, 2:]
        self.hidden_parent_id = np.asarray(hidden_parents, dtype=int)
        self.hidden_areas = np.array([_signed_area(t) for t in hidden_cut])
        self.degenerate_pieces = degeneracies
        if np.any(self.basis_values < 0) or not np.isfinite(cut).all():
            raise ValueError('invalid interpolated cut basis')
        ids = self.faces[self.parent_id]
        b = self.basis_values
        sums = b.sum(axis=1)
        local_mass = self.cut_areas[:, None, None]/12*(
            np.einsum('tvi,tvj->tij', b, b)+sums[:, :, None]*sums[:, None, :])
        if np.any((local_mass == 0) & (np.einsum('tvi,tvj->tij', b, b) > 0)):
            raise ValueError('positive cut mass underflow')
        rows = np.broadcast_to(ids[:, :, None], local_mass.shape).ravel()
        cols = np.broadcast_to(ids[:, None, :], local_mass.shape).ravel()
        full_mass = sparse.coo_matrix((local_mass.ravel(), (rows, cols)), shape=(nx*ny, nx*ny)).tocsr()
        self.interior_nodes = np.arange(nx*ny).reshape(ny, nx)[1:-1, 1:-1].ravel()
        diagonal = full_mass.diagonal()[self.interior_nodes]
        self.active_nodes = self.interior_nodes[diagonal > 0]
        self.covered_nodes = self.interior_nodes[diagonal == 0]
        self.mass = full_mass[self.active_nodes][:, self.active_nodes].tocsr()
        if not np.isfinite(self.mass.data).all() or np.any(diagonal < 0):
            raise ValueError('invalid physical mass')
        self.E = self._extension(harmonic_fill)
        self.physical_node_moments = np.bincount(ids.ravel(),
            weights=(self.cut_areas[:, None]*sums/3).ravel(), minlength=nx*ny)
        for array in (self.nodes, self.cut_triangles, self.basis_values, self.parent_id,
                      self.cut_areas, self.interior_nodes, self.active_nodes, self.covered_nodes,
                      self.physical_node_moments, self.hidden_triangles,
                      self.hidden_basis_values, self.hidden_parent_id, self.hidden_areas):
            array.setflags(write=False)

    def _extension(self, harmonic_fill):
        nx, ny = len(self.x), len(self.y)
        n = nx*ny
        active = {int(node): i for i, node in enumerate(self.active_nodes)}
        hidden = {int(node): i for i, node in enumerate(self.covered_nodes)}
        E = sparse.coo_matrix((np.ones(len(active)), (self.active_nodes, np.arange(len(active)))),
                              shape=(n, len(active))).tocsr()
        if not harmonic_fill or not hidden or not active:
            return E
        hr, hc, hv, br, bc, bv = [], [], [], [], [], []
        for node, row in hidden.items():
            iy, ix = divmod(node, nx)
            dx0, dx1 = self.x[ix]-self.x[ix-1], self.x[ix+1]-self.x[ix]
            dy0, dy1 = self.y[iy]-self.y[iy-1], self.y[iy+1]-self.y[iy]
            links = [(node-1, 2/(dx0*(dx0+dx1))), (node+1, 2/(dx1*(dx0+dx1))),
                     (node-nx, 2/(dy0*(dy0+dy1))), (node+nx, 2/(dy1*(dy0+dy1)))]
            hr.append(row); hc.append(row); hv.append(sum(weight for _, weight in links))
            for neighbor, weight in links:
                if neighbor in hidden:
                    hr.append(row); hc.append(hidden[neighbor]); hv.append(-weight)
                elif neighbor in active:
                    br.append(row); bc.append(active[neighbor]); bv.append(weight)
        A = sparse.coo_matrix((hv, (hr, hc)), shape=(len(hidden), len(hidden))).tocsc()
        boundary = np.unique(bc)
        if not len(boundary):
            return E
        B = sparse.coo_matrix((bv, (br, bc)), shape=(len(hidden), len(active))).tocsc()
        fill = splu(A).solve(B[:, boundary].toarray())
        if not np.isfinite(fill).all():
            raise ValueError('nonfinite harmonic extension')
        rr, cc = np.indices(fill.shape)
        return E+sparse.coo_matrix((fill.ravel(), (self.covered_nodes[rr.ravel()], boundary[cc.ravel()])),
                                    shape=E.shape).tocsr()

    def expand(self, coefficients):
        """Return full nodal values, shape (full_node_count, ...)."""
        values = np.asarray(coefficients)
        if values.ndim not in (1, 2) or values.shape[0] != len(self.active_nodes) or not np.isfinite(values).all():
            raise ValueError('finite active coefficient vector/matrix required')
        return self.E @ values

    def condense_maps(self, full_map):
        """Apply algebraic fill to a full-node potential/normal map."""
        values = np.asarray(full_map)
        if values.ndim != 2 or values.shape[1] != len(self.nodes) or not np.isfinite(values).all():
            raise ValueError('finite matrix with full nodal columns required')
        return (self.E.T @ values.T).T
