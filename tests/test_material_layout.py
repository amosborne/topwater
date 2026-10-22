"""Numerical material subdivision must preserve the physical wet hull."""
from dataclasses import replace
import json
import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.spatial.transform import Rotation
from topwater import SimulationInputs
from topwater._clipping import clip_flat_water
from topwater.frames import rotation_matrix
from topwater.hydrostatics import equilibrate_heave
from .reference_hull import input_data
from .experiments.free_surface.material_layout import material_mesh


@pytest.fixture(scope='module')
def canoe():
    data=SimulationInputs.from_json(json.dumps(input_data(2)))
    rest,_=equilibrate_heave(data.body,data.initial_state)
    return data,rest


@pytest.mark.parametrize('transverse',[False,True])
@pytest.mark.parametrize('heave,roll',[(0.,0.),(.002,.01),(-.003,-.015)])
def test_material_layout_preserves_wet_integrals_and_shared_edges(canoe,transverse,heave,roll):
    data,rest=canoe
    state=replace(rest,position=rest.position+np.array([0,0,heave]),
        orientation=np.roll(Rotation.from_rotvec([roll,0.,0.]).as_quat(),1))
    rotation=rotation_matrix(state.orientation)
    def clip(vertices,faces,physical):
        world=(vertices-data.body.cg)@rotation.T+state.position
        return clip_flat_water(world,faces,physical,0.,True)
    original=clip(data.body.mesh.vertices,data.body.mesh.faces,data.body.mesh.physical_faces)
    vertices,faces=material_mesh(data.body,transverse)
    refined=clip(vertices,faces,np.arange(len(faces)))
    # Volume, first moment, pressure force/moment, and wet area stay physical.
    for expected,actual in zip(original[3:],refined[3:]):
        assert_allclose(actual,expected,rtol=2e-12,atol=2e-14)
    wet=refined[0]
    vertices,inverse=np.unique(wet.reshape(-1,3),axis=0,return_inverse=True)
    faces=inverse.reshape(-1,3)
    edges=np.stack((faces,np.roll(faces,-1,axis=1)),axis=-1).reshape(-1,2)
    unique,inverse,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    assert np.all((counts==1)|(counts==2))
    assert np.all(vertices[unique[counts==1],2]==0)
    orientation=np.bincount(inverse,weights=np.sign(edges[:,1]-edges[:,0]))
    assert np.all(orientation[counts==2]==0)
