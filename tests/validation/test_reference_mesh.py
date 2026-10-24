"""Geometry invariants for the reference-only planar mesher."""
import numpy as np
import pytest
from .reference_mesh import facets,half_surface,moments,verify_surface,farthest_first,facet_width


def wet_box():
    vertices=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                       [-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0]],float)
    quads=((0,3,2,1),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
    return vertices[np.array([t for a,b,c,d in quads for t in ((a,b,c),(a,c,d))])]


def test_coplanar_reconstruction_and_reflection_keep_surface():
    original=wet_box();half=half_surface(original)
    reflected=half[:,[0,2,1]].copy();reflected[:,:,1]*=-1
    assert len(facets(original))==5
    quality=verify_surface(original,np.concatenate((half,reflected)))
    assert quality['facets']==5
    np.testing.assert_allclose(moments(original),moments(np.concatenate((half,reflected))),atol=1e-14)


def test_canonical_vertex_order_preserves_winding_and_maximizes_radius():
    original=np.array([[[0.,0.,-1.],[3.,0.,-1.],[0.,1.,-1.]]])
    canonical=farthest_first(original)
    for shift in range(3):
        np.testing.assert_array_equal(farthest_first(np.roll(original,shift,axis=1)),canonical)
    np.testing.assert_allclose(moments(canonical),moments(original),atol=1e-14)
    radii=np.linalg.norm(canonical-canonical.mean(axis=1)[:,None],axis=2)
    np.testing.assert_array_equal(radii[:,0],radii.max(axis=1))


def test_facet_width_uses_actual_narrow_geometry():
    polygon=np.array([[0.,0.,-1.],[1.,0.,-1.],[0.,.003,-1.]])
    expected=.003/np.sqrt(1+.003**2)
    assert facet_width(polygon,np.array([0.,0.,1.]))==pytest.approx(expected)
    moved=polygon+np.array([100.,-100.,50.])
    assert facet_width(moved,np.array([0.,0.,1.]))==pytest.approx(expected)


@pytest.mark.parametrize('defect',['missing','reversed','moved','overlap'])
def test_surface_proof_rejects_geometry_defects(defect):
    original=wet_box();changed=original.copy()
    if defect=='missing':changed=changed[:-1]
    if defect=='reversed':changed[0]=changed[0,[0,2,1]]
    if defect=='moved':changed[0,0,2]-=.001
    if defect=='overlap':changed=np.concatenate((changed,changed[:1]))
    with pytest.raises((AssertionError,ValueError)):
        verify_surface(original,changed)
