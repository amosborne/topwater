"""Fixed flat physical-exterior weak DtN, bounded CPU reference operator.

Evaluate contact rows by exterior-area work plus the waterline term. No hidden
physical mass, source remeshing, diagonal/row repair or symmetrization.
"""
from time import perf_counter
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from tests.laplace_reference import triangle_samples
from .cut_hat_geometry import cut_hat
from .physical_p1_patch import PhysicalP1Patch
from .plane_single_layer import plane_single_layer


class PhysicalWeakFlux:
    """Full-plane FFT action with directly integrated physical contact rows.

    ``kernel_owner`` is the caller's precomputed full-plane WeakPatch.
    Its interior dimensions/spacings must match patch.
    The physical footprint is checked against patch's saved cut geometry by
    rebuilding ONLY the integration geometry, never changing the source field.
    """
    def __init__(self, patch, footprint, kernel_owner, *, order=4,
                 maximum_pairs=500_000_000, maximum_bytes=512*1024**2):
        started = perf_counter()
        if not isinstance(patch, PhysicalP1Patch) or order not in (4, 8, 16):
            raise ValueError('physical patch and quadrature order 4,8,16 required')
        if any(isinstance(v, bool) or not isinstance(v, int) or v <= 0
               for v in (maximum_pairs, maximum_bytes)):
            raise ValueError('positive integer work/storage budgets required')
        self.patch = patch
        self.kernel_owner = kernel_owner
        self.order = order
        self.maximum_bytes = maximum_bytes
        x, y = patch.x, patch.y
        ny, nx = len(y)-2, len(x)-2
        values = kernel_owner.kernel
        hx, hy = kernel_owner.spacing
        if tuple(kernel_owner.shape) != (ny, nx):
            raise ValueError('kernel owner has a different interior shape')
        hx, hy = float(hx), float(hy)
        if (not np.isfinite([hx, hy]).all() or min(hx, hy) <= 0
                or not np.allclose(np.diff(x), hx, rtol=0, atol=32*np.finfo(float).eps*max(1., abs(x).max()))
                or not np.allclose(np.diff(y), hy, rtol=0, atol=32*np.finfo(float).eps*max(1., abs(y).max()))):
            raise ValueError('kernel must own the uniform patch spacings')
        values = np.asarray(values)
        if values.shape != (2*ny-1, 2*nx-1) or np.iscomplexobj(values) or not np.isfinite(values).all():
            raise ValueError('finite real signed kernel with matching patch interior required')
        # The supplied hole must own exactly the same cut-domain geometry.
        check = PhysicalP1Patch(x, y, footprint, harmonic_fill=False)
        for name in ('cut_triangles', 'parent_id', 'basis_values', 'hidden_triangles', 'hidden_parent_id'):
            if not np.array_equal(getattr(check, name), getattr(patch, name)):
                raise ValueError('footprint differs from physical patch cut owner')
        del check
        count, active_count = len(patch.interior_nodes), len(patch.active_nodes)
        full_to_active = np.full(len(patch.nodes), -1, dtype=int)
        full_to_active[patch.active_nodes] = np.arange(active_count)
        # Only active hats that overlap the hidden integration region require
        # correction. Covered-only target supports are algebraic, not rows.
        hidden_ids = patch.faces[patch.hidden_parent_id]
        touched = full_to_active[hidden_ids]
        affected = np.unique(touched[touched >= 0])
        geometry_bytes = sum(getattr(patch, name).nbytes for name in
            ('cut_triangles', 'basis_values', 'parent_id', 'cut_areas',
             'hidden_triangles', 'hidden_basis_values', 'hidden_parent_id', 'hidden_areas'))
        estimated_peak = 2*geometry_bytes+8*(2*len(affected)*len(patch.panels)
            +64*len(patch.panels)+3*len(affected)*len(patch.nodes)
            +2*len(affected)*active_count+32*len(patch.panels)+4*len(patch.nodes))
        if estimated_peak > maximum_bytes:
            raise MemoryError('weak physical constructor correction/pair scratch bound exceeds byte budget')
        interior_position = np.full(len(patch.nodes), -1, dtype=int)
        interior_position[patch.interior_nodes] = np.arange(count)
        self._active_interior = interior_position[patch.active_nodes]
        contact_x = np.zeros((len(affected), len(patch.panels)))
        contact_y = np.zeros_like(contact_x)
        active_to_affected = np.full(active_count, -1, dtype=int)
        active_to_affected[affected] = np.arange(len(affected))
        triangle = patch.panels.triangles[:, :, :2]
        gradient = np.linalg.solve(triangle[:, 1:]-triangle[:, :1],
            np.broadcast_to((np.eye(3)[1:]-np.eye(3)[:1]), (len(triangle), 2, 3)))
        pairs = 0
        scalar_counts = dict(fallback_pairs=0, exact_zero_distance_products=0)

        def single(points):
            nonlocal pairs
            work = len(points)*len(patch.panels)
            if pairs+work > maximum_pairs:
                raise RuntimeError('physical weak source-pair budget exhausted')
            pairs += work
            # A is a continuous scalar single layer on the COMMON plane.
            # Unlike a gradient trace, its finite edge/vertex values are valid.
            values,diagnostics=plane_single_layer(points,patch.panels)
            for name in scalar_counts:
                scalar_counts[name]+=diagnostics[name]
            return values

        # Direct physical area: no full-hat minus hidden-area cancellation.
        # Retain every positive piece supporting any affected test gradient.
        for cut_id, parent in enumerate(patch.parent_id):
            active = full_to_active[patch.faces[parent]]
            local = np.flatnonzero(active >= 0)
            local = local[active_to_affected[active[local]] >= 0]
            if not len(local):
                continue
            tri = np.column_stack((patch.cut_triangles[cut_id], np.zeros(3)))
            points, weights = triangle_samples(tri[None], order)
            integrated = np.zeros(len(patch.panels))
            for first in range(0, len(weights), 64):
                p = points[0, first:first+64, :2]
                integrated += patch.cut_areas[cut_id]*weights[first:first+64] @ single(p)
            for vertex in local:
                row = active_to_affected[active[vertex]]
                contact_x[row] += 2*gradient[parent, 0, vertex]*integrated
                contact_y[row] += 2*gradient[parent, 1, vertex]*integrated
        area_pairs = pairs
        # Keep original parent provenance for the test function. Scalar S
        # itself is defined at every common-plane query, including boundaries.
        nodes, weights = leggauss(order)
        t, weights = (nodes+1)/2, weights/2
        for affected_id, active in enumerate(affected):
            center = patch.nodes[patch.active_nodes[active]]
            geometry = cut_hat(center, hx, hy, footprint)
            center_row, center_col = divmod(patch.active_nodes[active], len(x))
            cells = geometry['parent_cells']+np.array([center_col, center_row])
            original_parents = (cells[:, 1]*(len(x)-1)+cells[:, 0]
                +geometry['parent_orientation']*(len(x)-1)*(len(y)-1))
            line = geometry['waterline']
            for segment, parent_local, normal, nodal in zip(line['segments'], line['parent_id'],
                    line['hole_outward_normals'], line['values'], strict=True):
                point = segment[:1]+t[:, None]*(segment[1:]-segment[:1])
                # Exact structured source-parent provenance from the hat's
                # affine triangle coordinates; no nearest-node assignment.
                target_triangle = geometry['parents'][parent_local]
                parent = original_parents[parent_local]
                if np.max(abs(triangle[parent]-target_triangle)) > 64*np.finfo(float).eps*max(1., abs(target_triangle).max()):
                    raise ValueError('waterline original-parent provenance differs')
                value = (1-t)*nodal[0]+t*nodal[1]
                integrated = (np.linalg.norm(segment[1]-segment[0])*weights*value) @ single(point)
                contact_x[affected_id] += 2*normal[0]*integrated
                contact_y[affected_id] += 2*normal[1]*integrated
        source_rows = np.repeat(np.arange(len(triangle)), 3)
        source_cols = patch.faces.ravel()
        gx = sparse.coo_matrix((gradient[:, 0].ravel(), (source_rows, source_cols)),
                               shape=(len(triangle), len(patch.nodes))).tocsr()
        gy = sparse.coo_matrix((gradient[:, 1].ravel(), (source_rows, source_cols)), shape=gx.shape).tocsr()
        contact = (gx.T @ contact_x.T+gy.T @ contact_y.T).T
        del contact_x, contact_y, gx, gy
        self.contact_rows = np.ascontiguousarray((patch.E.T @ contact.T).T)
        del contact
        if not np.isfinite(self.contact_rows).all():
            raise ValueError('nonfinite physical weak contact rows')
        self.affected_rows = affected
        self.shape = (active_count, active_count)
        for array in (self.contact_rows, self.affected_rows, self._active_interior):
            array.setflags(write=False)
        retained = self.contact_rows.nbytes+self.affected_rows.nbytes+self._active_interior.nbytes
        self.diagnostics = dict(assembly_seconds=perf_counter()-started,
            active_nodes=active_count, covered_nodes=len(patch.covered_nodes),
            contact_rows=len(affected), source_pairs=pairs, physical_area_pairs=area_pairs,
            waterline_pairs=pairs-area_pairs, retained_operator_bytes=retained,
            constructor_array_bound_bytes=estimated_peak, order=order,
            scalar_source_diagnostics=scalar_counts,
            storage_scope='Owned contact rows/maps; excludes caller-owned patch/kernel, FFT workspace and apply scratch',
            scope='Fixed flat physical weak loads; direct exterior contact rows replace FFT outputs, no left E transpose')

    def apply(self, active_coefficients):
        values = np.asarray(active_coefficients)
        if (values.ndim not in (1, 2) or values.shape[0] != self.shape[1]
                or np.iscomplexobj(values) or not np.isfinite(values).all()):
            raise ValueError('finite active coefficient vector/matrix required')
        vector = values.ndim == 1
        if vector:
            values = values[:, None]
        cases = values.shape[1]
        if cases < 1:
            raise ValueError('at least one coefficient column required')
        scratch = 8*cases*(len(self.patch.nodes)+len(self.patch.interior_nodes)
            +2*self.shape[0]+len(self.affected_rows))
        if scratch+self.diagnostics['retained_operator_bytes'] > self.maximum_bytes:
            raise MemoryError('physical weak host apply scratch exceeds byte budget; excludes FFT workspace')
        full = self.patch.expand(values)
        interior = full[self.patch.interior_nodes].reshape(*self.kernel_owner.shape, cases)
        weak = np.asarray(self.kernel_owner.apply_host(interior))
        if weak.shape != interior.shape or not np.isfinite(weak).all():
            raise FloatingPointError('full-plane weak action returned invalid loads')
        weak = weak.reshape(-1, cases)
        result = weak[self._active_interior].copy()
        result[self.affected_rows] = self.contact_rows @ values
        return result[:, 0] if vector else result
