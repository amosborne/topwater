"""Experimental exact vertex derivative of the current CUDA P1/P0 energy.

The numerical pair formula is shared with primal assembly. Eighteen forward
coordinate lanes, grouped in shared-primal batches, contract each pair's stationary scalar directly into vertex
atomics; there is no pair-by-coordinate output tensor. This differentiates the
existing quadrature approximation, not an independently qualified continuum
energy. CuPy remains optional, exactly as for cuda_operator.
"""
from pathlib import Path
import numpy as np
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from .cuda_operator import cp, _pairs
from cupyx import errstate
from cupyx.scipy.linalg import solve_triangular

# Inline the explicit shared include: NVRTC on Windows need not find source
# headers on a compiler include path. The numeric body has one source owner.
_folder=Path(__file__).parent
_source=(_folder/'cuda_shape.cu').read_text().replace(
    '#include "cuda_pair.cuh"',(_folder/'cuda_pair.cuh').read_text()).replace(
    '#include "cuda_dual.cuh"',(_folder/'cuda_dual.cuh').read_text()).replace(
    '#include "cuda_wide_dual.cuh"',(_folder/'cuda_wide_dual.cuh').read_text())
_kernel=cp.RawKernel(_source,'vertex_gradient',options=('--std=c++11',))
_kernels={1:_kernel}
DEFAULT_BATCH_WIDTH=3


def _geometric(points,faces,moving,*,order=3):
    tri=points[faces];dt=moving[faces]
    panel=TrianglePanels(tri,1)
    e=tri[:,1:]-tri[:,:1];de=dt[:,1:]-dt[:,:1]
    n=panel.geometric_normals
    da_vector=np.cross(de[:,0],e[:,1])+np.cross(e[:,0],de[:,1])
    d_twice=np.sum(n*da_vector,axis=1)
    dn=(da_vector-n*d_twice[:,None])/(2*panel.areas[:,None])
    edge=np.roll(tri,-1,axis=1)-tri;d_edge=np.roll(dt,-1,axis=1)-dt
    length=panel.edge_lengths;tangent=panel.tangents
    dl=np.sum(tangent*d_edge,axis=2)
    d_tangent=(d_edge-tangent*dl[:,:,None])/length[:,:,None]
    d_outward=np.cross(d_tangent,n[:,None])+np.cross(tangent,dn[:,None])
    opposite=np.roll(tri,-2,axis=1)-np.roll(tri,-1,axis=1)
    d_opposite=np.roll(dt,-2,axis=1)-np.roll(dt,-1,axis=1)
    gradient=np.cross(n[:,None],opposite)/(2*panel.areas[:,None,None])
    d_gradient=(np.cross(dn[:,None],opposite)+np.cross(n[:,None],d_opposite)-
        gradient*d_twice[:,None,None])/(2*panel.areas[:,None,None])
    curl=np.cross(n[:,None],gradient)
    d_curl=np.cross(dn[:,None],gradient)+np.cross(n[:,None],d_gradient)
    samples,weights=triangle_samples(tri,order)
    bary,_=triangle_samples(np.eye(3)[None],order)
    d_samples=np.einsum('qc,fci->fqi',bary[0],dt)
    values=[tri,n,panel.edge_outward,tangent,length,panel.areas,panel.centroids,
            panel.normal_sign,gradient,samples,weights]
    derivatives=[dt,dn,d_outward,d_tangent,dl,d_twice/2,dt.mean(axis=1),
                 np.zeros(len(faces)),d_gradient,d_samples,np.zeros_like(weights)]
    return values,derivatives,curl,d_curl


def _geometry_jacobians(points,faces,*,order=3):
    values,_,curl,_=_geometric(points,faces,np.zeros_like(points),order=order)
    isolated=points[faces].reshape(-1,3)
    numbering=np.arange(len(isolated)).reshape(-1,3)
    jacobians=[];curl_jacobians=[]
    for seed in range(9):
        direction=np.zeros((len(faces),3,3));direction[:,seed//3,seed%3]=1.
        _,derivative,_,dc=_geometric(isolated,numbering,direction.reshape(-1,3),order=order)
        jacobians.append(derivative[:-1]);curl_jacobians.append(dc)
    # Layout: face, nine independent local coordinates, remaining entries.
    packed=[]
    for index,value in enumerate(values[:-1]):
        packed.extend((cp.asarray(np.ascontiguousarray(value)),
            cp.asarray(np.ascontiguousarray(np.stack([entry[index] for entry in jacobians],axis=1)))))
    return packed,cp.asarray(values[-1]),curl,np.stack(curl_jacobians,axis=1)


def vertex_gradient(points,faces,potential,*,assembly=None,density=1000.,order=3,batch_width=DEFAULT_BATCH_WIDTH):
    """Return d(.5*rho*phi.T*E*phi)/d(points), phi fixed.

    Faces have the same complete fluid-outward winding as TrianglePanels(..., 1).
    A matching current ``prepare_energy`` assembly shares its factor/H; otherwise
    they are assembled here. The explicit order must match the retained assembly.
    Auxiliary q=V^-1 H phi is held stationary in the
    derivative. Density defaults to the legacy freshwater control. Source and target vertex
    contributions, including shared vertices and the self jump, accumulate
    once each. No normal-flux, row-sum or gauge correction is applied.
    """
    points=np.asarray(points,dtype=float)
    faces=np.asarray(faces)
    potential=np.asarray(potential,dtype=float)
    if (points.ndim!=2 or points.shape[1]!=3 or len(points)==0 or
            not np.all(np.isfinite(points)) or faces.ndim!=2 or faces.shape[1]!=3 or
            len(faces)==0 or faces.dtype.kind not in 'iu' or np.any(faces<0) or
            np.any(faces>=len(points)) or potential.shape!=(len(points),) or
            not np.all(np.isfinite(potential))):
        raise ValueError('finite points/potential and valid integer triangular faces required')
    density=float(density)
    if not np.isfinite(density) or density <= 0:
        raise ValueError('finite positive density required')
    if batch_width not in (1, 3, 6):
        raise ValueError('derivative batch width must be one, three or six')
    # Validate reuse BEFORE allocating or running the coordinate-derivative work.
    auxiliary = None if assembly is None else assembly.auxiliary(points,faces,potential,order=order)
    geometry,weights,curl,curl_jac=_geometry_jacobians(points,faces,order=order)
    nf,nv,nq=len(faces),len(points),len(weights)
    face_array=cp.asarray(np.ascontiguousarray(faces,dtype=np.int64))
    phi=cp.asarray(potential)
    if auxiliary is None:
        arguments=geometry[::2].copy()
        arguments.insert(8,face_array)
        arguments.append(weights)
        v=cp.empty((nf,nf),dtype=cp.float64)
        h=cp.zeros((nf,nv),dtype=cp.float64)
        _pairs(((nf*nf+127)//128,),(128,),
            (*arguments,np.int32(nf),np.int32(nv),np.int32(nq),v,h))
        v=(v+v.T)*.5
        with errstate(linalg='raise'):
            factor=cp.linalg.cholesky(v)
            auxiliary=solve_triangular(factor,h@phi,lower=True)
            auxiliary=solve_triangular(factor.T,auxiliary,lower=False)
    c=cp.asarray(np.einsum('fci,fc->fi',curl,potential[faces]))
    dc=cp.asarray(np.ascontiguousarray(np.einsum('fsci,fc->fsi',curl_jac,potential[faces])))
    gradient=cp.zeros((nv,3),dtype=cp.float64)
    if batch_width not in _kernels:
        _kernels[batch_width]=cp.RawKernel(_source,'vertex_gradient',
            options=('--std=c++11',f'-DSHAPE_WIDTH={batch_width}'))
    _kernels[batch_width](((nf*nf*(18//batch_width)+127)//128,),(128,),
        (*geometry,weights,face_array,c,dc,phi,auxiliary,
         np.int32(nf),np.int32(nv),np.int32(nq),np.float64(density),gradient))
    result=gradient.get()
    if not np.all(np.isfinite(result)):
        raise FloatingPointError('nonfinite CUDA energy vertex gradient')
    return result
