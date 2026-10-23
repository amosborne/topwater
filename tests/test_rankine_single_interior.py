"""Owned on-panel scalar integration; no invented boundary-gradient limits."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.integrate import quad
from scipy.spatial.transform import Rotation
from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels,panel_influence,single_layer_interior,single_layer_tangent


def independent_scalar(target,triangle):
    """Three positive fans; integrate radial coordinate analytically."""
    value=0.
    for a,b in zip(triangle-target,np.roll(triangle-target,-1,axis=0)):
        jacobian=np.linalg.norm(np.cross(a,b))
        radial,error=quad(lambda t:1/np.linalg.norm((1-t)*a+t*b),0.,1.,epsabs=1e-12,epsrel=1e-12)
        assert error<1e-10
        value+=jacobian*radial/(4*np.pi)
    return value


@pytest.mark.parametrize('sign',[1,-1])
def test_noncentroid_scalar_matches_independent_fan_integral_and_other_panels(sign):
    rotation=Rotation.from_rotvec([.31,-.24,.19]).as_matrix()
    triangle=np.array([[0.,0.,0.],[1.2,.1,0.],[.2,.9,0.]])@rotation.T+[.2,-.3,.4]
    panels=TrianglePanels(np.array([triangle,triangle+[2.,-.1,.3]]),sign)
    bary=np.array([[.17,.31,.52],[.61,.28,.11]])
    targets=bary@triangle
    result=single_layer_interior(targets,panels,np.zeros(2,dtype=int))
    assert_allclose(result[:,0],[independent_scalar(p,triangle) for p in targets],rtol=3e-12,atol=2e-14)
    other=panel_influence(targets,TrianglePanels(panels.triangles[1:],sign)).single_layer
    assert_allclose(result[:,1:],other,rtol=0,atol=0)
    tangent=single_layer_tangent(targets,panels,np.zeros(2,dtype=int))
    for row,point in enumerate(targets):
        expected=np.zeros(3)
        for edge in range(3):
            a,b=triangle[edge],triangle[(edge+1)%3]
            integral,_=quad(lambda t:np.linalg.norm(b-a)/(4*np.pi*np.linalg.norm(point-((1-t)*a+t*b))),
                            0.,1.,epsabs=1e-12,epsrel=1e-12)
            expected-=integral*panels.edge_outward[0,edge]
        assert_allclose(tangent[row,0],expected,rtol=3e-12,atol=3e-14)
    assert_allclose(np.einsum('mni,i->mn',tangent,panels.geometric_normals[0]),0.,atol=3e-15)
    direction=(triangle[1]-triangle[0])/np.linalg.norm(triangle[1]-triangle[0])
    errors=[]
    for step in [2e-4,1e-4]:
        derivative=(single_layer_interior(targets+step*direction,panels,[0,0])
                    -single_layer_interior(targets-step*direction,panels,[0,0]))/(2*step)
        errors.append(np.max(abs(derivative-tangent@direction)))
    assert errors[1]<.3*errors[0]+1e-10
    reverse=TrianglePanels(panels.triangles[::-1],sign)
    assert_allclose(single_layer_interior(targets,reverse,np.ones(2,dtype=int)),result[:,::-1],rtol=0,atol=0)
    with pytest.raises(ValueError,match='centroid'):
        panel_influence(targets,panels,self_indices=np.zeros(2,dtype=int))


def test_interior_owner_edges_plane_and_storage_guards():
    triangle=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    panels=TrianglePanels(np.array([triangle,triangle+[2.,0.,0.]]),1)
    with pytest.raises(ValueError,match='interior'):
        single_layer_interior([[.2,.3,0.]],panels,[1])
    for target in [[0.,0.,0.],[.5,0.,0.]]:
        with pytest.raises(ValueError,match='interior'):
            single_layer_interior([target],panels,[0])
        with pytest.raises(ValueError,match='interior'):
            single_layer_tangent([target],panels,[0])
    with pytest.raises(ValueError,match='plane'):
        single_layer_interior([[.2,.3,1e-4]],panels,[0])
    with pytest.raises(MemoryError,match='max_operator_bytes'):
        single_layer_interior([[.2,.3,0.]],panels,[0],controls=LaplaceControls(max_operator_bytes=1))
