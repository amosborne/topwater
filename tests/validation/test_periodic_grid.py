import numpy as np
import pytest
from .periodic_grid import hull_planes,radius,mean_angles,annulus,exterior_points,check_surface,mapped_surface


def box():
    from scipy.spatial import ConvexHull
    v=np.array([[x,y,z] for x in (-1.,1.) for y in (-1.,1.) for z in (-1.,1.)])
    h=ConvexHull(v);f=h.simplices.copy()
    for i,t in enumerate(v[f]):
        if np.dot(np.cross(t[1]-t[0],t[2]-t[0]),t.mean(axis=0))<0:f[i]=f[i,::-1]
    return v,f


def test_convex_sections_and_strict_reference_admissibility():
    v,f=box();p=hull_planes(v,f)
    np.testing.assert_allclose(radius(p,[0,np.pi/4],0.,np.zeros(2)),[1,np.sqrt(2)])
    exterior_points([[1.01,0,0],[0,0,-1.01]],p)
    for point in ([0,0,0],[1,0,0],[2,0,.001]):
        with pytest.raises(ValueError):exterior_points([point],p)
    with pytest.raises(ValueError):radius(p,[0],2.,np.zeros(2))


def test_contact_edges_must_preserve_hull_features():
    v,f=box();p=hull_planes(v,f)
    # Four axis rays lie on hull, but their connecting diamond cuts through it.
    a=np.arange(4)*np.pi/2;r=radius(p,a,0.,np.zeros(2))
    vertices,faces,_=annulus(a,r,4.,[0.,1.],np.zeros(2))
    with pytest.raises(ValueError,match='contact edge'):check_surface(vertices,faces,p,4)


def test_map_inserts_missing_hull_corners_and_has_no_interior_overlap():
    v,f=box();p=hull_planes(v,f);a=np.arange(8)*np.pi/4
    # Translate hull far enough to make its corner angles change; no pinned line.
    q=hull_planes(v+np.array([.1,.05,0]),f)
    eta=np.full((3,len(a)),.02)
    target,vertices,faces,checks=mapped_surface(a,eta,p,q,np.zeros(2),4.,[0.,1.])
    assert checks['inserted_angles']>=4 and checks['maximum_interior_overlap_m2']==0
    assert len(vertices)==3*len(target)
    np.testing.assert_allclose(vertices[:len(target),2],.02)


def test_nonconvex_geometry_is_rejected():
    v,f=box();v[0]=[0,0,0]
    with pytest.raises(ValueError,match='convex'):hull_planes(v,f)


def test_canonical_posed_hull_preserves_complete_contact_edges():
    from tests.reference_boat import reference_body,equilibrium
    from .periodic_fields import world_vertices
    from .mild import initial_body
    body=reference_body();state,loads=equilibrium()
    planes=hull_planes(world_vertices(body,state),body.mesh.faces);origin=state.position[:2]
    a=mean_angles(loads.wetted_triangles,origin,32)
    inner=radius(planes,a,0.,origin)
    v,f,_=annulus(a,inner,24.,[0.,.15,.4,1.,2.5,6.,12.],origin)
    check_surface(v,f,planes,len(a))
    posed=initial_body(state,np.array([.01,-.005,.015,.008,-.0005,.0002]),np.zeros(6))
    actual=hull_planes(world_vertices(body,posed),body.mesh.faces)
    eta=np.broadcast_to(.02*np.cos(a),(8,len(a)))
    target,mapped,faces,checks=mapped_surface(a,eta,planes,actual,origin,24.,[0.,.15,.4,1.,2.5,6.,12.])
    assert checks['maximum_interior_overlap_m2']==0
    assert checks['inserted_angles']>0
