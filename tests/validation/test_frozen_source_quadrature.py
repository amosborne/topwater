import numpy as np
import pytest
from .frozen_source_quadrature import integration_cells,owned_points,pressure_matrix,telescoping,duffy_rule,mirrored_gauss_targets


def test_descendant_centroids_are_owned_and_preserve_area_weights():
    t=np.array([[[0.,0.,-1.],[2.,0.,-1.],[0.,3.,-1.]]])
    for level in range(4):
        cells,owner=integration_cells(t,level)
        assert len(cells)==4**level and not owner.any()
        proof=owned_points(cells.mean(axis=1),owner,t,1e-10)
        assert proof['minimum_barycentric_weight']>0
        area=np.linalg.norm(np.cross(cells[:,1]-cells[:,0],cells[:,2]-cells[:,0]),axis=1)/2
        np.testing.assert_allclose(area,3/4**level,atol=1e-14)


@pytest.mark.parametrize('point',[[2,3,-1],[0,0,-1],[.3,.3,-.9],[.3,.3,.1]])
def test_invalid_or_edge_owned_targets_fail(point):
    t=np.array([[[0.,0.,-1.],[2.,0.,-1.],[0.,3.,-1.]]])
    with pytest.raises(ValueError):owned_points(np.array([point]),np.array([0]),t,1e-10)


def test_pressure_moment_uses_target_location_and_fluid_normal():
    point=np.array([[1.,2.,-1.]]);normal=np.array([[0.,0.,1.]]);phi=np.ones((1,6))*2
    result=pressure_matrix(point,normal,np.array([3.]),phi,np.zeros(3),1000.)
    expected=6000*np.array([0,0,1,2,-1,0])
    np.testing.assert_array_equal(result,np.broadcast_to(expected[:,None],(6,6)))


def test_telescope_preserves_cancellation_without_calling_remainder_source_error():
    coarse=np.zeros((6,6));fine=np.ones((6,6));quadrature=np.ones((6,6))*3
    total,output,remainder=telescoping(coarse,fine,quadrature,np.ones((6,6)),1e-10)
    np.testing.assert_array_equal(output,3);np.testing.assert_array_equal(remainder,-2)
    np.testing.assert_array_equal(total,output+remainder)


@pytest.mark.parametrize('order,degree',[(2,2),(4,4),(8,6)])
def test_duffy_gauss_exact_polynomial_moments(order,degree):
    from math import factorial
    bary,w=duffy_rule(order)
    assert np.all(bary>0) and np.all(w>0)
    np.testing.assert_allclose(bary.sum(axis=1),1.,atol=1e-14)
    for a in range(degree+1):
        for b in range(degree-a+1):
            for c in range(degree-a-b+1):
                expected=2*factorial(a)*factorial(b)*factorial(c)/factorial(a+b+c+2)
                actual=w@(bary[:,0]**a*bary[:,1]**b*bary[:,2]**c)
                assert actual==pytest.approx(expected,abs=1e-13)


def test_mirrored_targets_preserve_ownership_and_areas_without_mirroring_fields():
    t=np.array([[[0.,.2,-1.],[2.,.4,-1.],[0.,3.,-.5]]])
    triangles=np.r_[t,t[:,[1,0,2]]*np.array([1,-1,1])]
    points,owner,w,error=mirrored_gauss_targets(triangles,4)
    assert error==0 and len(points)==32
    np.testing.assert_array_equal(points[16:],points[:16]*np.array([1,-1,1]))
    np.testing.assert_allclose(np.bincount(owner,w),1.,atol=1e-14)
    owned_points(points,owner,triangles,1e-10)
    # A deliberately asymmetric scalar field must retain its odd contribution.
    from topwater.rankine import TrianglePanels
    panels=TrianglePanels(triangles,-1)
    potential=np.broadcast_to((1+.1*points[:,1])[:,None],(len(points),6))
    A=pressure_matrix(points,panels.normals[owner],w*panels.areas[owner],potential,np.zeros(3),1.)
    assert np.max(abs(A[1]))>1e-3
    bad=triangles.copy();bad[1,:,1]-=.1
    with pytest.raises(ValueError):mirrored_gauss_targets(bad,4)
