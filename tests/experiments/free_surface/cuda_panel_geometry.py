"""Shared current triangle geometry and coordinate Jacobians for CUDA energy kernels."""
import numpy as np
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from .cuda_operator import cp


def directional_geometry(points,faces,moving,*,order=3):
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


def geometry_jacobians(points,faces,*,order=3):
    values,_,curl,_=directional_geometry(points,faces,np.zeros_like(points),order=order)
    isolated=points[faces].reshape(-1,3)
    numbering=np.arange(len(isolated)).reshape(-1,3)
    jacobians=[];curl_jacobians=[]
    for seed in range(9):
        direction=np.zeros((len(faces),3,3));direction[:,seed//3,seed%3]=1.
        _,derivative,_,dc=directional_geometry(isolated,numbering,direction.reshape(-1,3),order=order)
        jacobians.append(derivative[:-1]);curl_jacobians.append(dc)
    # Layout: face, nine independent local coordinates, remaining entries.
    packed=[]
    for index,value in enumerate(values[:-1]):
        packed.extend((cp.asarray(np.ascontiguousarray(value)),
            cp.asarray(np.ascontiguousarray(np.stack([entry[index] for entry in jacobians],axis=1)))))
    return packed,cp.asarray(values[-1]),curl,np.stack(curl_jacobians,axis=1)


