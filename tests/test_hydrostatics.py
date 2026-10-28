"""Generic clipping/load controls on the canonical boat, without a fluid solve."""
from dataclasses import replace
import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.spatial.transform import Rotation
from topwater import HullMesh, InitialState
from topwater.hydrostatics import (
    FloodingError, LossOfContactError, NoEquilibriumError,
    equilibrate_heave, hydrostatics,
)
from .geometry_helpers import surface_integrals
from .reference_boat import DEPTH, GRAVITY, WATER_DENSITY, equilibrium, reference_body


def evaluated(body, state, **kwargs):
    return hydrostatics(body, state, density=WATER_DENSITY, gravity=GRAVITY, **kwargs)


def balanced(body, state=None, **kwargs):
    return equilibrate_heave(body, InitialState() if state is None else state,
        density=WATER_DENSITY, gravity=GRAVITY, **kwargs)


def heeled(rest):
    rotation=Rotation.from_rotvec([.08,.025,.4])*Rotation.from_quat(np.roll(rest.orientation,-1))
    return InitialState(rest.position+[1.2,-.7,.01],np.roll(rotation.as_quat(),1))


def test_pressure_and_all_six_world_load_derivatives_at_heeled_pose():
    body=reference_body();state=heeled(equilibrium()[0]);result=evaluated(body,state)
    expected=np.array([0.,0.,WATER_DENSITY*GRAVITY*result.displaced_volume])
    assert_allclose(result.force,expected,atol=1e-10)
    assert_allclose(result.moment,np.cross(result.center_of_buoyancy-state.position,expected),atol=1e-10)
    derivative=np.empty((6,6));step=1e-6
    rotation=Rotation.from_quat(np.roll(state.orientation,-1))
    for axis in range(6):
        loads=[]
        for sign in (-1,1):
            offset=np.eye(3)[axis%3]*step*sign
            moved=(replace(state,position=state.position+offset) if axis<3 else
                replace(state,orientation=np.roll((Rotation.from_rotvec(offset)*rotation).as_quat(),1)))
            value=evaluated(body,moved);loads.append(np.r_[value.force,value.moment])
        derivative[:,axis]=-(loads[1]-loads[0])/(2*step)
    assert_allclose(result.stiffness,derivative,atol=2e-5,rtol=1e-7)


def test_failure_boundaries_and_overloading():
    body=reference_body()
    with pytest.raises(FloodingError):
        evaluated(body,InitialState([0,0,body.cg[2]-DEPTH-.01]))
    with pytest.raises(LossOfContactError):
        evaluated(body,InitialState([0,0,10.]))
    with pytest.raises(NoEquilibriumError):
        balanced(replace(body,mass=1e6))


def test_closed_submerged_hull_volume_centroid_and_zero_heave_restoring():
    body=reference_body();mesh=replace(body.mesh,opening_rim=(),artificial_faces=())
    body=replace(body,mesh=mesh);state=InitialState([0,0,-2.])
    result=evaluated(body,state);volume,center=surface_integrals(mesh.vertices,mesh.faces)
    assert_allclose(result.displaced_volume,volume,rtol=1e-12)
    assert_allclose(result.center_of_buoyancy,center-body.cg+state.position,atol=1e-12)
    assert result.waterplane_area==0 and result.stiffness[2,2]==0
    assert result.rim_clearance is None


def test_equilibrium_is_independent_of_world_water_altitude():
    body=reference_body();rest,result=balanced(body,equilibrium()[0])
    level=1e9
    shifted,value=balanced(body,equilibrium()[0],water_level=level)
    assert_allclose(shifted.position[2]-level,rest.position[2],atol=np.spacing(level),rtol=0)
    # One representable world-height increment bounds the displacement error.
    force_allowance=abs(np.spacing(level))*result.stiffness[2,2]
    assert_allclose(value.force,result.force,atol=force_allowance,rtol=0)


def test_actual_vertex_waterline_preserves_displacement_derivative():
    body=reference_body();levels=np.unique(body.mesh.vertices[:,2]);level=levels[len(levels)//2]
    state=InitialState([0,0,body.cg[2]-level]);center=evaluated(body,state);step=1e-7
    below=evaluated(body,replace(state,position=state.position-[0,0,step]))
    above=evaluated(body,replace(state,position=state.position+[0,0,step]))
    assert_allclose((below.displaced_volume-above.displaced_volume)/(2*step),center.waterplane_area,rtol=1e-6)


@pytest.mark.parametrize('refinement',[0,1])
@pytest.mark.parametrize('pose',['equilibrium','heeled','vertex_row'])
def test_wet_shell_and_waterplane_join_by_exact_coordinates(refinement,pose):
    body=reference_body(refinement);rest=equilibrium(refinement)[0]
    if pose=='vertex_row':
        levels=np.unique(body.mesh.vertices[:,2])
        state=InitialState([0,0,body.cg[2]-levels[len(levels)//2]])
    else:
        state=rest if pose=='equilibrium' else heeled(rest)
    result=evaluated(body,state);lines=result.waterline_segments
    center=lines.reshape(-1,3).mean(axis=0)
    cap=np.stack((lines[:,1],lines[:,0],np.broadcast_to(center,(len(lines),3))),axis=1)
    triangles=np.concatenate((result.wetted_triangles,cap))
    vertices,indices=np.unique(triangles.reshape(-1,3),axis=0,return_inverse=True)
    closed=HullMesh(vertices,indices.reshape(-1,3));volume,centroid=surface_integrals(closed.vertices,closed.faces)
    assert_allclose(volume,result.displaced_volume,rtol=1e-12,atol=1e-14)
    assert_allclose(centroid,result.center_of_buoyancy,rtol=1e-12,atol=1e-13)
