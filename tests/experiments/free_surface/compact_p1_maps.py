"""Blocked body evaluation maps of the existing compact P1 Poisson lift.

Full nodal columns retain the original source-parent IDs. Outer perimeter
coefficients are fixed zero by the caller. No C-source trace is interpolated
into this space and no physical hidden-water mass is introduced.
"""
import numpy as np
from .panel_gradient_reference import checked_influence


def body_maps(trace,points,normals):
    points=np.asarray(points);normals=np.asarray(normals)
    if (points.ndim!=2 or points.shape[1]!=3 or normals.shape!=points.shape
            or np.iscomplexobj(points) or np.iscomplexobj(normals)
            or not np.isfinite(points).all() or not np.isfinite(normals).all()
            or np.any(points[:,2]>=0)):
        raise ValueError('finite real body points strictly below the source plane and matching normals required')
    panels,faces=trace.panels,trace.faces
    inverse=np.linalg.inv(panels.triangles[:,1:,:2]-panels.triangles[:,:1,:2])
    gradients=np.concatenate((-inverse.sum(axis=2)[:,None,:],inverse.transpose(0,2,1)),axis=1)
    count=len(trace.x)*len(trace.y)
    potential=np.zeros((len(points),count));normal=np.zeros_like(potential)
    for first in range(0,len(points),32):
        last=min(first+32,len(points));p=points[first:last];n=normals[first:last]
        influence=checked_influence(p,panels)
        affine=np.einsum('mti,tji->mtj',p[:,None,:2]-panels.triangles[None,:,0,:2],gradients)
        affine[:,:,0]+=1.
        tangent=np.einsum('mti,tji->mtj',influence.gradient_single_layer[:,:,:2],gradients)
        value=-2*(affine*influence.double_layer[:,:,None]+p[:,2,None,None]*tangent)
        derivative=-2*(influence.double_layer[:,:,None]*np.einsum('mi,tji->mtj',n[:,:2],gradients)
            +n[:,2,None,None]*tangent)
        for row in range(last-first):
            np.add.at(potential[first+row],faces.ravel(),value[row].ravel())
            np.add.at(normal[first+row],faces.ravel(),derivative[row].ravel())
    return potential,normal
