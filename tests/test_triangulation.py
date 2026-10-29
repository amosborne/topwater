"""Shared triangulation preserves physical boundary vertices."""
import numpy as np


def test_ear_triangulation_retains_collinear_vertices_and_concavity():
    from topwater.geometry.triangulation import _ears
    v=np.array([[0.,0.,0.],[1,0,0],[2,0,0],[2,1,0],[1,1,0],[1,2,0],[0,2,0],[0,1,0]])
    triangles=np.array(_ears(range(len(v)),v,np.array([0.,0.,1.])))
    reversed_triangles=np.array(_ears(list(reversed(range(len(v)))),v,np.array([0.,0.,1.])))
    np.testing.assert_array_equal(triangles,reversed_triangles)
    assert len(triangles)==len(v)-2 and set(triangles.ravel())==set(range(len(v)))
    t=v[triangles];area=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])[:,2].sum()/2
    np.testing.assert_allclose(area,3.,atol=1e-12)
