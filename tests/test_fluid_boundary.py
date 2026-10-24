import numpy as np
from scipy.spatial import ConvexHull
from topwater import HullMesh,RigidBody,InitialState
from topwater.fluid_boundary import closed_boundary
from topwater.surface import surface_integrals
from tests.validation.periodic_grid import annulus


def box_body():
    v=np.array([[x,y,z] for x in (-1.,1.) for y in (-1.,1.) for z in (-1.,1.)])
    h=ConvexHull(v);f=h.simplices.copy()
    for i,t in enumerate(v[f]):
        if np.dot(np.cross(t[1]-t[0],t[2]-t[0]),t.mean(axis=0))<0:f[i]=f[i,::-1]
    return RigidBody(HullMesh(v,f),4000.,[0,0,0],np.eye(3)*2000.)


def box_fluid(state=None,depth=4.):
    a=np.arange(8)*np.pi/4;r=1/np.maximum(abs(np.cos(a)),abs(np.sin(a)))
    v,f,_=annulus(a,r,4.,[0.,1.],np.zeros(2))
    return closed_boundary(box_body(),InitialState() if state is None else state,v,f,contact=np.arange(8),
        outer=np.arange(len(v)-8,len(v)),center_xy=np.zeros(2),depth=depth)


def test_closed_fluid_volume_and_hull_area():
    fluid=box_fluid();mesh=fluid.mesh
    volume,_=surface_integrals(mesh.vertices,mesh.faces)
    np.testing.assert_allclose(volume,8*16*np.sin(np.pi/4)/2*4-4,atol=1e-11)
    t=mesh.vertices[mesh.faces[fluid.tags=='body']]
    area=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1).sum()/2
    np.testing.assert_allclose(area,12.,atol=1e-12)
    edges=np.sort(np.concatenate((mesh.faces[:,[0,1]],mesh.faces[:,[1,2]],mesh.faces[:,[2,0]])),axis=1)
    assert np.all(np.unique(edges,axis=0,return_counts=True)[1]==2)


def test_canonical_calm_boundary_preserves_exact_displacement():
    from tests.reference_boat import reference_body,equilibrium
    from tests.validation.periodic_fields import world_vertices
    from tests.validation.periodic_grid import hull_planes,mean_angles,radius
    body=reference_body();state,loads=equilibrium();origin=state.position[:2]
    planes=hull_planes(world_vertices(body,state),body.mesh.faces)
    a=mean_angles(loads.wetted_triangles,origin,32);r=radius(planes,a,0.,origin)
    v,f,_=annulus(a,r,24.,[0.,.15,.4,1.,2.5,6.,12.],origin);n=len(a)
    fluid=closed_boundary(body,state,v,f,contact=np.arange(n),outer=np.arange(len(v)-n,len(v)),center_xy=origin,depth=24.)
    # Body-outward wet faces and the upward contact fan independently close
    # the displaced body. Their volume must retain the canonical equilibrium.
    wet=fluid.mesh.vertices[fluid.mesh.faces[fluid.tags=='body']][:,::-1]
    cap=np.array([[np.r_[origin,0.],v[i],v[(i+1)%n]] for i in range(n)])
    triangles=np.r_[wet,cap]
    volume=np.einsum('ni,ni->n',triangles[:,0],np.cross(triangles[:,1],triangles[:,2])).sum()/6
    np.testing.assert_allclose(volume,body.mass/1025.,rtol=1e-11,atol=1e-11)


def test_canonical_nonflat_contact_preserves_concave_wet_facets():
    from tests.reference_boat import reference_body,equilibrium
    from tests.validation.periodic_fields import world_vertices
    from tests.validation.periodic_grid import hull_planes,mean_angles,mapped_surface
    from tests.validation.mild import initial_body
    body=reference_body();state,loads=equilibrium();origin=state.position[:2]
    planes=hull_planes(world_vertices(body,state),body.mesh.faces)
    a=mean_angles(loads.wetted_triangles,origin,32)
    posed=initial_body(state,np.array([.01,-.005,.015,.008,-.0005,.0002]),np.zeros(6))
    actual=hull_planes(world_vertices(body,posed),body.mesh.faces)
    eta=np.broadcast_to(.02*np.cos(a)+.008*np.cos(2*a),(8,len(a)))
    target,v,f,checks=mapped_surface(a,eta,planes,actual,origin,24.,[0.,.15,.4,1.,2.5,6.,12.])
    n=len(target)
    fluid=closed_boundary(body,posed,v,f,contact=np.arange(n),outer=np.arange(len(v)-n,len(v)),center_xy=origin,depth=24.)
    assert fluid.diagnostics['panels_by_tag']['body']>0
    wet=fluid.mesh.vertices[fluid.mesh.faces[fluid.tags=='body']]
    # Every body triangle must remain on one actual physical support plane.
    residual=np.max(abs(wet@actual[:,:3].T-actual[:,3]),axis=1)
    assert np.all(np.min(residual,axis=1)<1e-10)


def test_ear_triangulation_retains_collinear_vertices_and_concavity():
    from topwater.fluid_boundary import _ears
    v=np.array([[0.,0.,0.],[1,0,0],[2,0,0],[2,1,0],[1,1,0],[1,2,0],[0,2,0],[0,1,0]])
    triangles=np.array(_ears(range(len(v)),v,np.array([0.,0.,1.])))
    reversed_triangles=np.array(_ears(list(reversed(range(len(v)))),v,np.array([0.,0.,1.])))
    np.testing.assert_array_equal(triangles,reversed_triangles)
    assert len(triangles)==len(v)-2 and set(triangles.ravel())==set(range(len(v)))
    t=v[triangles];area=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])[:,2].sum()/2
    np.testing.assert_allclose(area,3.,atol=1e-12)


def test_bottom_must_not_cut_hull():
    import pytest
    with pytest.raises(ValueError,match='bottom intersects'):box_fluid(depth=.5)
