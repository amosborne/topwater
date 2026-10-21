"""Experimental fused Galerkin integration; no certified error estimates.

This specializes the constant-triangle edge/solid-angle formulas and their
first source moments. It omits production cancellation guards/fallbacks and
must not replace topwater._rankine_kernel without further qualification.
"""
from math import sqrt,atan2,log1p,pi
import numpy as np
from numba import njit


@njit(cache=True)
def symmetrize_single_layer(v):
    """Measure raw quadrature disagreement and average in the owned buffer."""
    largest, disagreement = 0., 0.
    for i in range(len(v)):
        largest = max(largest, abs(v[i, i]))
        for j in range(i):
            left, right = v[i, j], v[j, i]
            largest = max(largest, abs(left), abs(right))
            disagreement = max(disagreement, abs(left-right))
            average = (left+right)/2
            v[i, j], v[j, i] = average, average
    return disagreement/largest


@njit(cache=True)
def assemble_pairs(triangles,normals,outward,tangents,lengths,areas,centroids,signs,
                   faces,gradient,points,weights,nv):
    nf=len(triangles);v=np.zeros((nf,nf));hmat=np.zeros((nf,nv))
    r=np.empty((3,3));distance=np.empty(3);edge_gradient=np.empty(3)
    for target_face in range(nf):
        for source_face in range(nf):
            for point in range(len(weights)):
                x=points[target_face,point]
                for edge in range(3):
                    squared=0.
                    for axis in range(3):
                        r[edge,axis]=x[axis]-triangles[source_face,edge,axis]
                        squared+=r[edge,axis]**2
                    distance[edge]=sqrt(squared)
                height=0.
                for axis in range(3):
                    height+=r[0,axis]*normals[source_face,axis]
                    edge_gradient[axis]=0.
                self_pair=target_face==source_face
                if self_pair:
                    height=0.
                angle=0.
                if not self_pair:
                    denominator=distance[0]*distance[1]*distance[2]
                    for edge in range(3):
                        dot=0.
                        for axis in range(3):
                            dot+=r[edge,axis]*r[(edge+1)%3,axis]
                        denominator+=dot*distance[(edge+2)%3]
                    angle=2*atan2(2*areas[source_face]*height,denominator)
                terms=0.
                for edge in range(3):
                    d=0.;a=0.
                    for axis in range(3):
                        d-=r[edge,axis]*outward[source_face,edge,axis]
                        a+=r[edge,axis]*tangents[source_face,edge,axis]
                    length=lengths[source_face,edge];end=length-a;perp=d*d+height*height
                    denominator=perp/(distance[edge]+abs(a))+perp/(distance[(edge+1)%3]+abs(end))+2*(max(-a,0.)+max(a-length,0.))
                    logarithm=log1p(2*length/denominator)
                    terms+=d*logarithm
                    for axis in range(3):
                        edge_gradient[axis]+=outward[source_face,edge,axis]*logarithm
                single=(terms-height*angle)/(4*pi)
                double=angle*signs[source_face]/(4*pi)
                weight=areas[target_face]*weights[point]
                v[target_face,source_face]+=weight*single
                for corner in range(3):
                    linear=double/3
                    if not self_pair:
                        for axis in range(3):
                            offset=x[axis]-centroids[source_face,axis]-height*normals[source_face,axis]
                            moment=offset*double-height*edge_gradient[axis]*signs[source_face]/(4*pi)
                            linear+=moment*gradient[source_face,corner,axis]
                    hmat[target_face,faces[source_face,corner]]+=weight*linear
        for corner in range(3):
            hmat[target_face,faces[target_face,corner]]+=areas[target_face]/6
    return v,hmat
