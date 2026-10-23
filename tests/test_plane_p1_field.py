"""Original parent ownership and simultaneous real trace columns."""
import numpy as np
from numpy.testing import assert_allclose,assert_array_equal
import pytest
from tests.experiments.free_surface.compact_p1_trace import CompactP1Trace
from tests.experiments.free_surface.plane_p1_field import PlaneP1Field


def test_original_parent_lookup_and_multicolumn_centroid_trace():
    x=np.linspace(1.3,3.9,5);y=np.linspace(-2.7,-.9,4)
    values=np.zeros((len(y),len(x),2))
    values[1:-1,1:-1]=np.array([[[.3,-.7],[1.1,.2],[-.4,.8]],
        [[.5,.9],[-.2,.1],[.8,-.6]]])
    field=PlaneP1Field(x,y,values)
    targets=np.array([.19,.23,.58])@field.panels.triangles
    owners=field.point_owners(targets[:,:2])
    assert_array_equal(owners,np.arange(len(field.panels)))
    centres=field.panels.centroids
    ids=field.point_owners(centres[:,:2])
    actual=field.normal_flux(centres[:,:2],ids)
    for column in range(2):
        reference=CompactP1Trace(x,y,values[:,:,column]).surface_velocity(ids)[:,2]
        assert_allclose(actual[:,column],reference,rtol=2e-12,atol=3e-14)
    edge=np.mean(field.panels.triangles[0,:2,:2],axis=0)[None]
    with pytest.raises(ValueError,match='interior'):
        field.vector_potential(edge,field.point_owners(edge))
    with pytest.raises(ValueError,match='inside'):
        field.point_owners([[x[-1]+1.,y[1]]])
