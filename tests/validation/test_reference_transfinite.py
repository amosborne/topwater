"""Geometry gates for fixed affine triangular patches, without native solves."""
import numpy as np
import pytest
from .reference_transfinite import lattice_triangles,verify_patch


SPEC=dict(barycentric_lattice_absolute=1e-7,child_area_relative=1e-8,child_area_absolute_m2=1e-14)


@pytest.mark.parametrize('n',[2,3,4,5])
def test_affine_grid_preserves_patch_and_connectivity(n):
    parent=np.array([[.5,0.,-1.],[2.,1.,-.5],[-1.,2.,-.5]])
    uv=lattice_triangles(n)/n
    children=parent[0]+uv[...,0,None]*(parent[1]-parent[0])+uv[...,1,None]*(parent[2]-parent[0])
    assert verify_patch(parent,children[::-1],n,SPEC)<1e-12


@pytest.mark.parametrize('defect',['displaced','off_plane','reversed','missing','duplicate'])
def test_affine_gate_rejects_fallback_or_defective_mesh(defect):
    parent=np.array([[0.,0.,-1.],[2.,0.,-1.],[0.,1.,-1.]])
    uv=lattice_triangles(3)/3
    children=parent[0]+uv[...,0,None]*(parent[1]-parent[0])+uv[...,1,None]*(parent[2]-parent[0])
    if defect=='displaced':children[0,1,0]+=.01
    if defect=='off_plane':children[0,1,2]+=.01
    if defect=='reversed':children[0]=children[0,[0,2,1]]
    if defect=='missing':children=children[:-1]
    if defect=='duplicate':children[0]=children[1]
    with pytest.raises((AssertionError,ValueError)):
        verify_patch(parent,children,3,SPEC)
