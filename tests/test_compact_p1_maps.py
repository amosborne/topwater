"""Nodal scatter maps against the independently exercised field evaluator."""
import numpy as np
from numpy.testing import assert_allclose
from tests.experiments.free_surface.compact_p1_trace import CompactP1Trace
from tests.experiments.free_surface.compact_p1_maps import body_maps


def test_multicolumn_body_maps_preserve_full_source_node_order():
    x=np.linspace(-1.3,1.7,5);y=np.linspace(-.8,1.1,4)
    values=np.zeros((4,5,2));values[1:-1,1:-1]=np.array(
        [[[.4,-.3],[.1,.7],[-.5,.2]],[[.8,.4],[-.2,-.6],[.3,.9]]])
    geometry=CompactP1Trace(x,y,values[:,:,0])
    points=np.array([[.13,-.17,-.21],[-.43,.31,-.67],[2.3,1.4,-1.1]])
    normals=np.array([[.2,-.3,.6],[-.4,.1,-.7],[.5,-.8,.3]])
    p,n=body_maps(geometry,points,normals)
    actual_p=p@values.reshape(-1,2);actual_n=n@values.reshape(-1,2)
    for column in range(2):
        expected_p,expected_u=CompactP1Trace(x,y,values[:,:,column]).bulk(points)
        assert_allclose(actual_p[:,column],expected_p,rtol=3e-12,atol=4e-15)
        assert_allclose(actual_n[:,column],np.sum(expected_u*normals,axis=1),rtol=3e-12,atol=4e-15)
    # Independently differentiate map potential in the supplied normal direction.
    errors=[]
    for step in (2e-4,1e-4):
        plus=body_maps(geometry,points+step*normals,normals)[0]@values.reshape(-1,2)
        minus=body_maps(geometry,points-step*normals,normals)[0]@values.reshape(-1,2)
        errors.append(np.max(abs((plus-minus)/(2*step)-actual_n)))
    assert errors[1]<.3*errors[0]+1e-11
