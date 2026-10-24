"""Geometry invariants for the reference-only planar mesher."""
import numpy as np
import pytest
from .reference_mesh import facets,half_surface,moments,verify_surface


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


@pytest.mark.parametrize('defect',['missing','reversed','moved','overlap'])
def test_surface_proof_rejects_geometry_defects(defect):
    original=wet_box();changed=original.copy()
    if defect=='missing':changed=changed[:-1]
    if defect=='reversed':changed[0]=changed[0,[0,2,1]]
    if defect=='moved':changed[0,0,2]-=.001
    if defect=='overlap':changed=np.concatenate((changed,changed[:1]))
    with pytest.raises((AssertionError,ValueError)):
        verify_surface(original,changed)
