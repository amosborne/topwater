"""Experimental fixed sample-local analytic stationary scalar pullback.

No runtime tape, approximate operator, or new triangle geometry owner. The
shared sample primitive records primal values once; the specific pullback
chains feature adjoints through the existing exact geometry Jacobians.
"""
from pathlib import Path
import numpy as np
from .cuda_operator import cp
from .cuda_shape_gradient import _geometry_jacobians,_source as _shape_source

_source=_shape_source+'\n'+Path(__file__).with_name('cuda_scalar_reverse.cu').read_text()
_kernel=cp.RawKernel(_source,'reverse_vertex_gradient',options=('--std=c++11',))


def vertex_gradient(points,faces,potential,*,assembly,density=1000.,order=3):
    points=np.asarray(points,dtype=float);faces=np.asarray(faces);phi=np.asarray(potential,dtype=float)
    if (points.ndim!=2 or points.shape[1]!=3 or not len(points) or not np.isfinite(points).all()
            or faces.ndim!=2 or faces.shape[1]!=3 or not len(faces) or faces.dtype.kind not in 'iu'
            or np.any(faces<0) or np.any(faces>=len(points)) or phi.shape!=(len(points),)
            or not np.isfinite(phi).all()):
        raise ValueError('finite current points/potential and valid triangular faces required')
    density=float(density)
    if not np.isfinite(density) or density<=0:
        raise ValueError('positive finite density required')
    auxiliary=assembly.auxiliary(points,faces,phi,order=order)
    geometry,weights,curl,curl_jac=_geometry_jacobians(points,faces,order=order)
    nf,nv,nq=len(faces),len(points),len(weights)
    c=cp.asarray(np.einsum('fci,fc->fi',curl,phi[faces]))
    dc=cp.asarray(np.ascontiguousarray(np.einsum('fsci,fc->fsi',curl_jac,phi[faces])))
    device_faces=cp.asarray(np.ascontiguousarray(faces,dtype=np.int64));device_phi=cp.asarray(phi)
    local_phi=device_phi[device_faces]
    # Linear contractions of the EXISTING source P1 gradient/Jacobian owner.
    g=cp.ascontiguousarray(cp.einsum('fci,fc->fi',geometry[16],local_phi))
    dg=cp.ascontiguousarray(cp.einsum('fsci,fc->fsi',geometry[17],local_phi))
    result=cp.zeros((nv,3),dtype=cp.float64)
    _kernel(((nf*nf+127)//128,),(128,),
        (*geometry,weights,device_faces,c,dc,g,dg,device_phi,auxiliary,
         np.int32(nf),np.int32(nv),np.int32(nq),np.float64(density),result))
    result=result.get()
    if not np.isfinite(result).all():
        raise FloatingPointError('nonfinite sample-reverse vertex gradient')
    return result
