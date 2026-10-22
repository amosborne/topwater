"""Exact-current FP64 quadrature-metric action candidate, no fluid solver.

Caches six weighted reference-metric coefficients per cell/q point. No local
35-by-35 matrices, inverse, or geometry update occurs inside an action. This
cache describes ONE current geometry; lagging it changes the fluid operator.
"""
from time import perf_counter
import numpy as np

from .p4_geometry import quadrature, InvalidGeometryError


class MetricAction:
    def __init__(self, coordinates, nodes, node_count, *, order=8):
        import cupy as cp
        import cupyx

        started = perf_counter()
        coordinates = np.asarray(coordinates, dtype=np.float64)
        nodes = np.asarray(nodes)
        if (coordinates.ndim != 3 or coordinates.shape[1:] != (35, 3)
                or nodes.shape != coordinates.shape[:2] or nodes.dtype.kind not in 'iu'
                or not np.all(np.isfinite(coordinates)) or node_count < 1
                or node_count >= 2**31 or np.any(nodes < 0) or np.any(nodes >= node_count)):
            raise ValueError('finite current cells[N,35,3] with matching valid integer node IDs required')
        self.node_count = node_count;self.cp = cp
        _, derivative, weights = quadrature(order)
        stream = cp.cuda.get_current_stream();stream.synchronize()
        upload = perf_counter()
        x = cp.asarray(coordinates);d = cp.asarray(derivative);w = cp.asarray(weights)
        self.nodes = cp.asarray(nodes, dtype=cp.int32)
        self.derivative = d.transpose(1, 0, 2).reshape(35, -1)
        self.metric = cp.empty((len(x), len(weights), 6), dtype=cp.float64)
        minimum = cp.empty(len(x), dtype=cp.float64)
        stream.synchronize();upload = perf_counter()-upload
        compute = perf_counter()
        with cupyx.errstate(linalg='raise'):
            for start in range(0, len(x), 256):
                stop = min(start+256, len(x))
                jac = d[None].transpose(0, 1, 3, 2)@x[start:stop, None]
                determinant = cp.linalg.det(jac)
                if bool(cp.any((determinant <= 0)|~cp.isfinite(determinant))):
                    raise InvalidGeometryError('nonpositive or invalid sampled P4 geometry Jacobian')
                inverse = cp.linalg.inv(jac)
                tensor = (inverse.swapaxes(-1, -2)@inverse)*(determinant*w)[..., None, None]
                for column, (i, j) in enumerate(((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))):
                    self.metric[start:stop, :, column] = tensor[..., i, j]
                minimum[start:stop] = cp.min(determinant, axis=1)
        stream.synchronize();compute = perf_counter()-compute
        if not bool(cp.all(cp.isfinite(self.metric))):
            raise FloatingPointError('nonfinite current reference metric')
        self.minimum_sampled_jacobian = cp.asnumpy(minimum)
        self.setup_timing = dict(upload_seconds=upload, metric_seconds=compute,
                                 total_seconds=perf_counter()-started)
        self.persistent_array_bytes = int(self.metric.nbytes+self.nodes.nbytes+self.derivative.nbytes)

    def action(self, values, *, chunk_size=None):
        """Full nodal host vector in/out; None means ONE whole-cell batch."""
        cp = self.cp;started = perf_counter()
        values = np.asarray(values, dtype=np.float64)
        if values.shape != (self.node_count,) or not np.all(np.isfinite(values)):
            raise ValueError('finite full nodal vector required')
        count = len(self.nodes)
        chunk = count if chunk_size is None else chunk_size
        if not isinstance(chunk, int) or chunk < 1:
            raise ValueError('positive integer action chunk required')
        stream = cp.cuda.get_current_stream();stream.synchronize()
        upload = perf_counter();x = cp.asarray(values)
        result = cp.zeros(self.node_count, dtype=cp.float64)
        stream.synchronize();upload = perf_counter()-upload
        compute = perf_counter()
        for start in range(0, count, chunk):
            stop = min(start+chunk, count)
            gradient = (x[self.nodes[start:stop]]@self.derivative).reshape(stop-start, -1, 3)
            m = self.metric[start:stop]
            gx, gy, gz = gradient[..., 0], gradient[..., 1], gradient[..., 2]
            flux = cp.stack((m[..., 0]*gx+m[..., 3]*gy+m[..., 4]*gz,
                             m[..., 3]*gx+m[..., 1]*gy+m[..., 5]*gz,
                             m[..., 4]*gx+m[..., 5]*gy+m[..., 2]*gz), axis=-1)
            local = flux.reshape(stop-start, -1)@self.derivative.T
            cp.add.at(result, self.nodes[start:stop].ravel(), local.ravel())
        stream.synchronize();compute = perf_counter()-compute
        download = perf_counter();host = cp.asnumpy(result)
        stream.synchronize();download = perf_counter()-download
        if not np.all(np.isfinite(host)):
            raise FloatingPointError('nonfinite current metric action')
        return host, dict(upload_seconds=upload, action_seconds=compute,
                          download_seconds=download, total_seconds=perf_counter()-started)
