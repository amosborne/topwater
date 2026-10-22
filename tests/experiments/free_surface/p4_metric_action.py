"""Exact-current FP64 quadrature-metric action candidate, no fluid solver.

Caches six weighted reference-metric coefficients per cell/q point. No local
35-by-35 matrices, inverse, or geometry update occurs inside an action. This
cache describes ONE current geometry; lagging it changes the fluid operator.
"""
from time import perf_counter
import numpy as np

from .p4_geometry import quadrature, InvalidGeometryError


class TangentMetricError(ValueError):
    """The approximate tangent tensor lost positivity; refresh, not geometry failure."""


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
        self.coordinates = x
        self.inverse = cp.empty((len(x), len(weights), 3, 3), dtype=cp.float64)
        self.measure = cp.empty((len(x), len(weights)), dtype=cp.float64)
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
                self.inverse[start:stop] = inverse
                self.measure[start:stop] = determinant*w
                tensor = (inverse.swapaxes(-1, -2)@inverse)*self.measure[start:stop, :, None, None]
                for column, (i, j) in enumerate(((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))):
                    self.metric[start:stop, :, column] = tensor[..., i, j]
                minimum[start:stop] = cp.min(determinant, axis=1)
        stream.synchronize();compute = perf_counter()-compute
        if not bool(cp.all(cp.isfinite(self.metric))):
            raise FloatingPointError('nonfinite current reference metric')
        self.minimum_sampled_jacobian = cp.asnumpy(minimum)
        self.setup_timing = dict(upload_seconds=upload, metric_seconds=compute,
                                 total_seconds=perf_counter()-started)
        self.persistent_array_bytes = int(sum(array.nbytes for array in (self.metric, self.nodes, self.derivative, self.coordinates, self.inverse, self.measure)))

    def action(self, values, *, chunk_size=None):
        """Exact BASE-geometry action; full nodal host vector in/out."""
        return self._action(values, self.metric, chunk_size=chunk_size)

    def _flux(self, gradient, metric):
        cp = self.cp
        gx, gy, gz = gradient[..., 0], gradient[..., 1], gradient[..., 2]
        return cp.stack((metric[..., 0]*gx+metric[..., 3]*gy+metric[..., 4]*gz,
                         metric[..., 3]*gx+metric[..., 1]*gy+metric[..., 5]*gz,
                         metric[..., 4]*gx+metric[..., 5]*gy+metric[..., 2]*gz), axis=-1)

    def _action(self, values, metric, *, chunk_size=None):
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
            flux = self._flux(gradient, metric[start:stop])
            local = flux.reshape(stop-start, -1)@self.derivative.T
            cp.add.at(result, self.nodes[start:stop].ravel(), local.ravel())
        stream.synchronize();compute = perf_counter()-compute
        download = perf_counter();host = cp.asnumpy(result)
        stream.synchronize();download = perf_counter()-download
        if not np.all(np.isfinite(host)):
            raise FloatingPointError('nonfinite current metric action')
        return host, dict(upload_seconds=upload, action_seconds=compute,
                          download_seconds=download, total_seconds=perf_counter()-started)


    def kinetic(self, values):
        """Exact base energy/full fixed-phi physical shape force, same q8 cache."""
        return self._kinetic(values)

    def _kinetic(self, values, tangent_metric=None):
        cp = self.cp;started = perf_counter()
        values = np.asarray(values, dtype=np.float64)
        if values.shape != (self.node_count,) or not np.all(np.isfinite(values)):
            raise ValueError('finite full nodal potential required')
        x = cp.asarray(values);force = cp.zeros((self.node_count, 3), dtype=cp.float64)
        energies = cp.empty(len(self.nodes), dtype=cp.float64)
        for start in range(0, len(self.nodes), 256):
            stop = min(start+256, len(self.nodes));n = stop-start
            reference = (x[self.nodes[start:stop]]@self.derivative).reshape(n, -1, 3)
            inverse, measure = self.inverse[start:stop], self.measure[start:stop]
            velocity = (inverse@reference[..., None])[..., 0]
            squared = cp.sum(velocity*velocity, axis=-1)
            dual = (inverse.swapaxes(-1, -2)@velocity[..., None])[..., 0]
            tensor = measure[..., None, None]*(.5*squared[..., None, None]*inverse
                                               -velocity[..., :, None]*dual[..., None, :])
            local = 1000*(self.derivative[None]@tensor.swapaxes(-1, -2).reshape(n, -1, 3))
            for axis in range(3):
                cp.add.at(force[:, axis], self.nodes[start:stop].ravel(), local[..., axis].ravel())
            if tangent_metric is None:
                energies[start:stop] = 500*cp.sum(measure*squared, axis=1)
            else:
                flux = self._flux(reference, tangent_metric[start:stop])
                energies[start:stop] = 500*cp.sum(reference*flux, axis=(1, 2))
        host_force = cp.asnumpy(force);energy = float(cp.sum(energies))
        if not np.isfinite(energy) or not np.all(np.isfinite(host_force)):
            raise FloatingPointError('nonfinite cached kinetic energy or shape force')
        return dict(kinetic_energy=energy, kinetic_gradient=host_force,
                    total_seconds=perf_counter()-started)

    def tangent(self, coordinates):
        """Full first-order operator at new X, with its matching BASE shape force.

        Current geometry positivity is independent of tangent tensor positivity.
        No automatic refresh, geometry rejection backoff, or energy correction.
        """
        cp = self.cp;started = perf_counter()
        coordinates = np.asarray(coordinates, dtype=np.float64)
        if coordinates.shape != self.coordinates.shape or not np.all(np.isfinite(coordinates)):
            raise ValueError('finite current coordinates matching cached cells required')
        current = cp.asarray(coordinates);metric = cp.empty_like(self.metric)
        derivative = self.derivative.reshape(35, -1, 3).transpose(1, 0, 2)
        minimum = cp.empty(len(self.nodes), dtype=cp.float64)
        pivots = cp.empty((len(self.nodes), 3), dtype=cp.float64)
        strain = cp.empty(len(self.nodes), dtype=cp.float64)
        for start in range(0, len(self.nodes), 256):
            stop = min(start+256, len(self.nodes))
            dj = derivative[None].transpose(0, 1, 3, 2)@(current[start:stop]-self.coordinates[start:stop])[:, None]
            current_j = derivative[None].transpose(0, 1, 3, 2)@current[start:stop, None]
            determinant = cp.linalg.det(current_j)
            if bool(cp.any((determinant <= 0)|~cp.isfinite(determinant))):
                raise InvalidGeometryError('nonpositive or invalid CURRENT sampled P4 Jacobian')
            inverse = self.inverse[start:stop]
            c = inverse@dj  # Transpose of grad_X0(delta X), under row-J convention.
            trace = cp.trace(c, axis1=-2, axis2=-1)
            tensor = (1+trace)[..., None, None]*cp.eye(3)-c-c.swapaxes(-1, -2)
            first = tensor[..., 0, 0]
            second = first*tensor[..., 1, 1]-tensor[..., 0, 1]**2
            third = cp.linalg.det(tensor)
            if bool(cp.any((first <= 0)|(second <= 0)|(third <= 0)|~cp.isfinite(third))):
                raise TangentMetricError('tangent tensor is not positive definite; refresh approximation at this state')
            reference = (inverse.swapaxes(-1, -2)@tensor@inverse)*self.measure[start:stop, :, None, None]
            for column, (i, j) in enumerate(((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))):
                metric[start:stop, :, column] = reference[..., i, j]
            minimum[start:stop] = cp.min(determinant, axis=1)
            pivots[start:stop] = cp.stack((cp.min(first, axis=1), cp.min(second, axis=1), cp.min(third, axis=1)), axis=1)
            strain[start:stop] = cp.max(cp.sqrt(cp.sum(c*c, axis=(-1, -2))), axis=1)
        if not bool(cp.all(cp.isfinite(metric))):
            raise FloatingPointError('nonfinite tangent reference metric')
        return TangentAction(self, metric, dict(
            minimum_current_jacobian=cp.asnumpy(minimum),
            minimum_tangent_leading_minors=cp.asnumpy(pivots),
            maximum_displacement_gradient_frobenius=float(cp.max(strain)),
            total_seconds=perf_counter()-started))


class TangentAction:
    """K0+DK0[deltaX] and matching derivative; force is NOT exact-current force."""
    def __init__(self, base, metric, diagnostics):
        self.base = base;self.metric = metric;self.diagnostics = diagnostics

    def action(self, values, *, chunk_size=None):
        return self.base._action(values, self.metric, chunk_size=chunk_size)

    def kinetic(self, values):
        """Tangent quadratic energy; base-geometry/current-phi shape gradient."""
        return self.base._kinetic(values, self.metric)
