"""Geometric conservation of the experimental fluid pocket construction."""
import numpy as np
import pytest
from .experiments.free_surface.pocket import grid, boundary


@pytest.mark.parametrize('n,amplitude,shift',[(6,0.,0),(12,.08,3),(24,0.,7)])
def test_pocket_boundary_encloses_its_tetrahedra(n,amplitude,shift):
    points,tets,near,_=grid(n,amplitude,shift)
    faces=boundary(points,tets[near])
    triangles=points[faces]
    surface_volume=np.einsum('ij,ij->',triangles[:,0],
        np.cross(triangles[:,1],triangles[:,2]))/6
    coordinates=points[tets[near]]
    cell_volume=np.abs(np.linalg.det(coordinates[:,1:]-coordinates[:,:1])).sum()/6
    assert surface_volume==pytest.approx(cell_volume,rel=1e-13)
    # A closed oriented surface has two oppositely directed copies of each edge.
    edges=np.stack((faces,np.roll(faces,-1,axis=1)),axis=-1).reshape(-1,2)
    _,inverse,counts=np.unique(np.sort(edges,axis=1),axis=0,
                               return_inverse=True,return_counts=True)
    assert np.all(counts==2)
    assert np.all(np.bincount(inverse,weights=np.sign(edges[:,1]-edges[:,0]))==0)
    if amplitude==0:
        box_volume=np.ptp(points[:,0])*np.ptp(points[:,1])*4
        assert cell_volume==pytest.approx(box_volume/9,rel=1e-13)
        all_coordinates=points[tets]
        total=np.abs(np.linalg.det(all_coordinates[:,1:]-all_coordinates[:,:1])).sum()/6
        assert total==pytest.approx(box_volume,rel=1e-13)
