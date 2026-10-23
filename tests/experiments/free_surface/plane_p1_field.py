"""Bounded CPU integration reference for one compact COMMON-PLANE P1 field.

A=sum grad(f_t)*S_t and Df=-2 sum grad(f_t) dot tangent_grad(S_t).
No curved-graph claim and no arbitrary boundary-gradient limit is introduced.
"""
from time import perf_counter
import numpy as np
from topwater.rankine import single_layer_interior,single_layer_tangent
from .compact_p1_trace import CompactP1Trace


class PlaneP1Field:
    def __init__(self,x,y,values,maximum_pairs=500_000_000):
        values=np.asarray(values)
        if values.ndim!=3 or np.iscomplexobj(values) or not np.isfinite(values).all():
            raise ValueError('finite real (y,x,columns) nodal values required')
        if np.any(values[[0,-1],:]) or np.any(values[:,[0,-1]]):
            raise ValueError('compact field must be exactly zero on the patch perimeter')
        owner=CompactP1Trace(x,y,values[:,:,0])
        self.panels=owner.panels;self.faces=owner.faces
        self.x=owner.x;self.y=owner.y
        nodal=np.array(values.reshape(-1,values.shape[-1])[self.faces],copy=True)
        edges=self.panels.triangles[:,1:,:2]-self.panels.triangles[:,:1,:2]
        self.gradient=np.linalg.solve(edges,nodal[:,1:]-nodal[:,:1])
        self.gradient.setflags(write=False)
        self.maximum_pairs=maximum_pairs;self.pairs=0;self.seconds=0.

    def source_parents(self,row,column,geometry):
        """Original integer cell/orientation provenance, never nearest queries."""
        cells=geometry['parent_cells']+np.array([column,row])
        nx=len(self.x)-1;ny=len(self.y)-1
        if np.any(cells<0) or np.any(cells[:,0]>=nx) or np.any(cells[:,1]>=ny):
            raise ValueError('complete hat parents must lie in the source patch')
        indices=cells[:,1]*nx+cells[:,0]+geometry['parent_orientation']*nx*ny
        allowance=256*np.finfo(float).eps*max(1.,float(abs(geometry['parents']).max()))
        if np.max(abs(self.panels.triangles[indices,:,:2]-geometry['parents']))>allowance:
            raise ValueError('cut parent provenance differs from the source mesh')
        return indices

    def _evaluate(self,points,owners,tangent):
        points=np.asarray(points,dtype=float);owners=np.asarray(owners)
        if points.ndim!=2 or points.shape[1:]!=(2,) or owners.shape!=(len(points),):
            raise ValueError('one source owner per planar integration point required')
        pairs=len(points)*len(self.panels)
        if self.pairs+pairs>self.maximum_pairs:
            raise RuntimeError('plane P1 integration-reference source-pair budget exhausted')
        self.pairs+=pairs
        started=perf_counter()
        columns=self.gradient.shape[-1]
        result=np.empty((len(points),columns) if tangent else (len(points),2,columns))
        for first in range(0,len(points),64):
            last=min(first+64,len(points))
            targets=np.column_stack((points[first:last],np.zeros(last-first)))
            if tangent:
                gradient=single_layer_tangent(targets,self.panels,owners[first:last])
                result[first:last]=-2*np.einsum('mni,nic->mc',gradient[:,:,:2],self.gradient)
            else:
                single=single_layer_interior(targets,self.panels,owners[first:last])
                result[first:last]=(single@self.gradient.reshape(len(self.panels),-1)).reshape(-1,2,columns)
        self.seconds+=perf_counter()-started
        return result

    def vector_potential(self,points,owners):
        return self._evaluate(points,owners,False)

    def normal_flux(self,points,owners):
        return self._evaluate(points,owners,True)
