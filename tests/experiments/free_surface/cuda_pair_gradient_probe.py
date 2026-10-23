"""Tiny selected-pair diagnostic of the SAME production gradient kernels.

Explicit compile macro changes only output/restricted work selection; physical
sample and pullback owners are unchanged. Never allocate nf*nf*18 outputs.
"""
import numpy as np
from .cuda_operator import cp
from .cuda_shape_gradient import _geometry_jacobians
from .cuda_scalar_gradient import _source as _forward_source
from .cuda_reverse_gradient import _source as _reverse_source

_forward=cp.RawKernel(_forward_source,'contracted_vertex_gradient',
    options=('--std=c++11','-DSHAPE_WIDTH=1','-DPAIR_DIAGNOSTIC'))
_reverse=cp.RawKernel(_reverse_source,'reverse_vertex_gradient',
    options=('--std=c++11','-DPAIR_DIAGNOSTIC'))


def selected_pair_derivatives(points,faces,potential,pairs,*,assembly,density=1000.,order=3):
    """Return (forward,reverse) arrays, rows=(target9,source9,value)."""
    points=np.asarray(points,dtype=float);faces=np.asarray(faces);phi=np.asarray(potential,dtype=float)
    pairs=np.asarray(pairs)
    if (pairs.ndim!=2 or pairs.shape[1]!=2 or pairs.dtype.kind not in 'iu'
            or not 1<=len(pairs)<=64 or np.any(pairs<0) or np.any(pairs>=len(faces))):
        raise ValueError('one to sixty-four valid integer face pairs required')
    if (phi.shape!=(len(points),) or not np.isfinite(phi).all()
            or not np.isfinite(density) or density<=0):
        raise ValueError('finite current nodal potential and positive density required')
    z=assembly.auxiliary(points,faces,phi,order=order)
    geometry,weights,curl,curl_jac=_geometry_jacobians(points,faces,order=order)
    nf,nv,nq=len(faces),len(points),len(weights)
    device_faces=cp.asarray(np.ascontiguousarray(faces,dtype=np.int64));device_phi=cp.asarray(phi)
    c=cp.asarray(np.einsum('fci,fc->fi',curl,phi[faces]))
    dc=cp.asarray(np.ascontiguousarray(np.einsum('fsci,fc->fsi',curl_jac,phi[faces])))
    g=cp.ascontiguousarray(cp.einsum('fci,fc->fi',geometry[16],device_phi[device_faces]))
    dg=cp.ascontiguousarray(cp.einsum('fsci,fc->fsi',geometry[17],device_phi[device_faces]))
    selected=cp.asarray(np.ascontiguousarray(pairs[:,0]*nf+pairs[:,1],dtype=np.int64))
    forward=cp.empty((len(pairs),19),dtype=cp.float64);reverse=cp.empty_like(forward)
    prefix=(*geometry,weights,device_faces,c,dc)
    suffix=(device_phi,z,np.int32(nf),np.int32(nv),np.int32(nq),np.float64(density))
    _forward(((len(pairs)*18+127)//128,),(128,),(*prefix,*suffix,forward,selected,np.int32(len(pairs))))
    _reverse(((len(pairs)+127)//128,),(128,),(*prefix,g,dg,*suffix,reverse,selected,np.int32(len(pairs))))
    forward,reverse=forward.get(),reverse.get()
    if not np.isfinite(forward).all() or not np.isfinite(reverse).all():
        raise FloatingPointError('nonfinite selected-pair diagnostic')
    return forward,reverse
