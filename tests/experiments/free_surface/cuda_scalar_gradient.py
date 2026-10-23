"""Contained forward reference for the exactly contracted stationary scalar.

The nonlinear sample formula remains in cuda_pair.cuh for both raw and scalar
adapters. This is a correctness gate for a later sample-local pullback, not a
new driver backend or an assumed speed improvement. Current phi/auxiliary and
all curl/gradient geometry dependencies retain the existing owners.
"""
from pathlib import Path
import numpy as np
from .cuda_operator import cp
from .cuda_panel_geometry import geometry_jacobians
from .cuda_shape_gradient import _source as _shape_source

_source=_shape_source+'\n'+Path(__file__).with_name('cuda_scalar_shape.cu').read_text()
_kernels={}


def vertex_gradient(points,faces,potential,*,assembly,density=1000.,order=3,batch_width=1):
    """Same current energy gradient, contracting before forming local H columns."""
    points=np.asarray(points,dtype=float);faces=np.asarray(faces);phi=np.asarray(potential,dtype=float)
    if (points.ndim!=2 or points.shape[1]!=3 or not len(points)
            or not np.isfinite(points).all() or faces.ndim!=2 or faces.shape[1]!=3
            or not len(faces) or faces.dtype.kind not in 'iu' or np.any(faces<0)
            or np.any(faces>=len(points)) or phi.shape!=(len(points),) or not np.isfinite(phi).all()):
        raise ValueError('finite current points/potential and valid triangular faces required')
    density=float(density)
    if not np.isfinite(density) or density<=0 or batch_width not in (1,3):
        raise ValueError('positive finite density and derivative width one or three required')
    auxiliary=assembly.auxiliary(points,faces,phi,order=order)
    geometry,weights,curl,curl_jac=geometry_jacobians(points,faces,order=order)
    nf,nv,nq=len(faces),len(points),len(weights)
    c=cp.asarray(np.einsum('fci,fc->fi',curl,phi[faces]))
    dc=cp.asarray(np.ascontiguousarray(np.einsum('fsci,fc->fsi',curl_jac,phi[faces])))
    result=cp.zeros((nv,3),dtype=cp.float64)
    if batch_width not in _kernels:
        _kernels[batch_width]=cp.RawKernel(_source,'contracted_vertex_gradient',
            options=('--std=c++11',f'-DSHAPE_WIDTH={batch_width}'))
    _kernels[batch_width](((nf*nf*(18//batch_width)+127)//128,),(128,),
        (*geometry,weights,cp.asarray(np.ascontiguousarray(faces,dtype=np.int64)),c,dc,
         cp.asarray(phi),auxiliary,np.int32(nf),np.int32(nv),np.int32(nq),np.float64(density),result))
    result=result.get()
    if not np.isfinite(result).all():
        raise FloatingPointError('nonfinite contracted-scalar vertex derivative')
    return result
