"""Planar polygon triangulation shared by physical hull ingestion."""
import numpy as np


def _cross(a,b):
    return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def _ears(indices,vertices,normal):
    """Planar polygon triangulation retaining every boundary vertex."""
    axis=np.eye(3)[np.argmin(abs(normal))]
    u=np.cross(axis,normal);u/=np.linalg.norm(u);v=np.cross(normal,u)
    p=vertices@np.stack((u,v),axis=1)
    loop=list(indices)
    if sum(_cross(p[loop[i]],p[loop[(i+1)%len(loop)]]) for i in range(len(loop)))<0:loop.reverse()
    first=loop.index(min(loop));loop=loop[first:]+loop[:first]
    scale=max(np.ptp(p[loop],axis=0));eps=64*np.finfo(float).eps*scale**2
    triangles=[]
    while len(loop)>3:
        for i in range(len(loop)):
            a,b,c=loop[i-1],loop[i],loop[(i+1)%len(loop)]
            if _cross(p[b]-p[a],p[c]-p[a])<=eps:continue
            other=[j for j in loop if j not in (a,b,c)]
            q=p[other]
            inside=(_cross(p[b]-p[a],q-p[a])>=-eps)&(_cross(p[c]-p[b],q-p[b])>=-eps)&(_cross(p[a]-p[c],q-p[c])>=-eps)
            if np.any(inside):continue
            triangles.append([a,b,c]);loop.pop(i);break
        else:raise ValueError('wet facet polygon cannot be triangulated without losing a boundary vertex')
    if _cross(p[loop[1]]-p[loop[0]],p[loop[2]]-p[loop[0]])<=eps:
        raise ValueError('degenerate final wet-facet triangle')
    return triangles+[loop]
