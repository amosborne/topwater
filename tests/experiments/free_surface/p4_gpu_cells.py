"""Optional FP64 chunked cell candidate; no global operator or fluid solver.

The SAME current-cell q8 metrics supply local stiffness and kinetic shape work.
CuPy is imported only when this explicitly experimental function is invoked.
"""
from time import perf_counter
import numpy as np

from .p4_geometry import quadrature, InvalidGeometryError


def cell_terms(coordinates, potential, *, order=8, chunk_size=64):
    """Host arrays in/out; return unchanged physical cell quantities plus costs."""
    import cupy as cp
    import cupyx

    started = perf_counter()
    coordinates = np.asarray(coordinates, dtype=np.float64)
    potential = np.asarray(potential, dtype=np.float64)
    if (coordinates.ndim != 3 or coordinates.shape[1:] != (35, 3)
            or potential.shape != coordinates.shape[:2]
            or not np.all(np.isfinite(coordinates)) or not np.all(np.isfinite(potential))
            or not isinstance(chunk_size, int) or chunk_size < 1 or chunk_size > 1024):
        raise ValueError('finite cells[N,35,3], potential[N,35], chunk1..1024 required')
    _, derivative, weights = quadrature(order)
    stream = cp.cuda.get_current_stream();stream.synchronize()
    upload = perf_counter()
    x = cp.asarray(coordinates);phi = cp.asarray(potential)
    d = cp.asarray(derivative);w = cp.asarray(weights)
    count = len(x);nq = len(weights)
    stiffness = cp.empty((count, 35, 35), dtype=cp.float64)
    force = cp.empty((count, 35, 3), dtype=cp.float64)
    energy = cp.empty(count, dtype=cp.float64)
    minimum = cp.empty(count, dtype=cp.float64)
    stream.synchronize();upload = perf_counter()-upload
    compute = perf_counter()
    with cupyx.errstate(linalg='raise'):
        for start in range(0, count, chunk_size):
            stop = min(start+chunk_size, count);n = stop-start
            jac = d[None].transpose(0, 1, 3, 2)@x[start:stop, None]
            determinant = cp.linalg.det(jac)
            if bool(cp.any((determinant <= 0)|~cp.isfinite(determinant))):
                raise InvalidGeometryError('nonpositive or invalid sampled P4 geometry Jacobian')
            gradient = d[None]@cp.linalg.inv(jac).swapaxes(-1, -2)
            measure = determinant*w[None]
            velocity = cp.einsum('ba,bqai->bqi', phi[start:stop], gradient)
            squared = cp.sum(velocity*velocity, axis=2)
            directional = cp.einsum('bqai,bqi->bqa', gradient, velocity)
            first = gradient.reshape(n, nq, 105).transpose(0, 2, 1)@(.5*measure*squared)[..., None]
            second = (directional*measure[..., None]).transpose(0, 2, 1)@velocity
            force[start:stop] = 1000*(first.reshape(n, 35, 3)-second)
            energy[start:stop] = 500*cp.sum(measure*squared, axis=1)
            minimum[start:stop] = cp.min(determinant, axis=1)
            flat = gradient.transpose(0, 2, 1, 3).reshape(n, 35, nq*3)
            stiffness[start:stop] = (flat*cp.repeat(measure, 3, axis=1)[:, None])@flat.transpose(0, 2, 1)
    stream.synchronize();compute = perf_counter()-compute
    download = perf_counter()
    result = dict(element_stiffness=cp.asnumpy(stiffness), kinetic_energy=cp.asnumpy(energy),
                  kinetic_gradient=cp.asnumpy(force), minimum_sampled_jacobian=cp.asnumpy(minimum))
    stream.synchronize();download = perf_counter()-download
    if not all(np.all(np.isfinite(value)) for value in result.values()):
        raise FloatingPointError('nonfinite FP64 cell outputs')
    result['timing'] = dict(upload_seconds=upload, compute_seconds=compute,
                            download_seconds=download, total_seconds=perf_counter()-started)
    return result
