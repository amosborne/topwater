"""One current-geometry directional energy ACTION, reusing current V/H.

Differentiate the same quadrature and symmetrization as prepare_energy. No
finite-difference fluid operators, derivative factorization or full dE/dW.
"""
from pathlib import Path
import numpy as np
from cupyx import errstate
from cupyx.scipy.linalg import solve_triangular

from .cuda_operator import cp
from .cuda_panel_geometry import directional_geometry

_folder = Path(__file__).parent
_source = '\n'.join((_folder/name).read_text() for name in
                     ('cuda_dual.cuh', 'cuda_pair.cuh', 'cuda_directional.cu'))
_pairs = cp.RawKernel(_source, 'directional_pairs', options=('--std=c++11',))


def energy_direction_action(points, faces, moving, potential, *, assembly):
    """Return (dE/dX[moving])@potential, with nodal potential held fixed.

    V is the production SYMMETRIZED single layer. q=inv(V) H phi is current.
    dE phi = dW phi + dH.T q + H.T inv(V)(dH phi-dV q).
    W=C.T V C; its product rule is applied via three face-curl components.
    """
    points, moving, potential = [np.asarray(value, dtype=float) for value in
                                (points, moving, potential)]
    if (not assembly.uniform_quadrature
            or moving.shape != points.shape or potential.shape != (len(points),)
            or not np.isfinite(moving).all() or not np.isfinite(potential).all()
            or not np.array_equal(faces, assembly.faces)
            or not np.array_equal(points[faces], assembly.triangles)
            or np.any(assembly.normal_sign != 1)):
        raise ValueError('directional action requires uniform quadrature and identical current fluid geometry')
    values, derivatives, curl, dcurl = directional_geometry(points, faces, moving,
                                                order=assembly.quadrature_order)
    geometry = []
    for value, derivative in zip(values[:-1], derivatives[:-1], strict=True):
        geometry.extend(cp.asarray(np.ascontiguousarray(a)) for a in (value, derivative))
    weights = cp.asarray(values[-1])
    nf, nv, nq = len(faces), len(points), len(values[-1])
    indices = cp.asarray(np.ascontiguousarray(faces, dtype=np.int64))
    dv = cp.empty((nf, nf), dtype=cp.float64)
    dh = cp.zeros((nf, nv), dtype=cp.float64)
    _pairs(((nf*nf+127)//128,), (128,), (*geometry, weights, indices,
        np.int32(nf), np.int32(nv), np.int32(nq), dv, dh))
    dv = (dv+dv.T)*.5
    phi = cp.asarray(potential)
    curl, dcurl = cp.asarray(curl), cp.asarray(dcurl)
    cphi = cp.sum(curl*phi[indices][:, :, None], axis=1)
    dcphi = cp.sum(dcurl*phi[indices][:, :, None], axis=1)
    factor, trace = assembly.factor, assembly.mixed_trace

    def inverse(rhs):
        with errstate(linalg='raise'):
            return solve_triangular(factor.T, solve_triangular(factor, rhs, lower=True), lower=False)

    auxiliary = inverse(trace@phi)
    # Raw V was intentionally discarded by the primal owner. Its current
    # Cholesky factor applies V without reconstructing another dense matrix.
    vcphi = factor@(factor.T@cphi)
    variation = dv@cphi+factor@(factor.T@dcphi)
    local = cp.sum(dcurl*vcphi[:, None, :]+curl*variation[:, None, :], axis=2)
    action = cp.zeros(nv, dtype=cp.float64)
    cp.add.at(action, indices.ravel(), local.ravel())
    action += dh.T@auxiliary+trace.T@inverse(dh@phi-dv@auxiliary)
    result = action.get()
    if not np.isfinite(result).all():
        raise FloatingPointError('nonfinite current directional energy action')
    return result, dict(quadrature_order=assembly.quadrature_order,
        derivative_trace_constant_norm=float(cp.max(cp.abs(dh@cp.ones(nv))).get()),
        derivative_array_bytes=dict(single_layer=dv.nbytes, mixed_trace=dh.nbytes),
        full_energy_derivative_assembled=False)
