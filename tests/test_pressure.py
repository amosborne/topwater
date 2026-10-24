import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from topwater import InitialState
from topwater.rankine import TrianglePanels
from topwater.pressure import PressureSystem,pressure_acceleration,body_normal_acceleration,hydrostatic_load
from .test_fluid_boundary import box_body,box_fluid


def test_material_normal_derivative_includes_rotation_and_hessian():
    state=InitialState(position=[.2,-.1,.3],linear_velocity=[.3,-.4,.2],angular_velocity=[.2,.1,-.3])
    r=np.array([.7,-.2,.5]);n=np.array([1.,2.,3.]);n/=np.linalg.norm(n)
    H=np.array([[.3,.1,.2],[.1,-.2,.05],[.2,.05,-.1]])
    U=np.array([.4,-.2,.1]);omega=state.angular_velocity
    x=state.position+r;u=U+H@x
    actual=body_normal_acceleration(x[None],n[None],state,u[None],H[None])[0]
    def residual(t):
        R=Rotation.from_rotvec(omega*t).as_matrix();rt=R@r
        point=state.position+state.linear_velocity*t+rt
        return np.dot(state.linear_velocity+np.cross(omega,rt)-U-H@point,R@n)
    step=1e-5;expected=(residual(step)-residual(-step))/(2*step)
    assert actual==pytest.approx(expected,abs=1e-10)
    with pytest.raises(ValueError):body_normal_acceleration(x[None],n[None],state,None,H[None])


def test_exact_hydrostatic_moment_matches_degree_two_triangle_quadrature():
    t=np.array([[[0.,0.,-1.],[2.,0.,-.3],[0.,3.,-.1]]]);n=np.array([[.2,-.3,.4]]);n/=np.linalg.norm(n)
    center=np.array([.4,-.2,.1]);area=np.linalg.norm(np.cross(t[0,1]-t[0,0],t[0,2]-t[0,0]))/2
    bary=np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]])
    x=bary@t[0];force=-1000*9.81*x[:,2,None]*n
    expected=area*np.r_[force.mean(axis=0),np.cross(x-center,force).mean(axis=0)]
    np.testing.assert_allclose(hydrostatic_load(t,n,center,1000.,9.81),expected,rtol=1e-13,atol=1e-10)


def test_owned_pressure_coupling_balances_all_six_nonzero_loads():
    fluid=box_fluid();p=TrianglePanels.from_mesh(fluid.mesh);body=fluid.tags=='body'
    owner=PressureSystem(fluid,fluid.tags=='surface',length_scale=4.)
    load=np.array([3.,-4.,5.,2.,-1.,.7]);count=len(p.areas)
    result=pressure_acceleration(owner,known_potential=np.zeros(count),
        known_normal_acceleration=np.zeros(count),body_speed_squared=np.zeros(body.sum()),
        density=1000.,gravity=9.81,applied_load=load)
    assert np.count_nonzero(abs(result.acceleration_world)>1e-8)==6
    assert result.added_mass.shape==(6,6)
    mass=np.diag([4000.]*3+[2000.]*3)
    np.testing.assert_allclose(mass@result.acceleration_world,result.fluid_load+result.external_load,atol=1e-10)
    # Reconstruct pressure force/moment independently of the owner's load matrix.
    t=p.triangles[body];area=p.areas[body];normal=p.normals[body]
    bary=np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]])
    points=np.einsum('qk,nkj->nqj',bary,t)
    pressure=-1000*(result.acceleration_potential[body,None]+9.81*points[:,:,2])
    force=pressure[:,:,None]*normal[:,None,:]
    independent=np.r_[np.sum(area[:,None]*force.mean(axis=1),axis=0),
        np.sum(area[:,None]*np.cross(points,force).mean(axis=1),axis=0)]
    np.testing.assert_allclose(result.fluid_load,independent,atol=1e-10)
    rest=pressure_acceleration(owner,known_potential=np.zeros(count),
        known_normal_acceleration=np.zeros(count),body_speed_squared=np.zeros(body.sum()),density=1000.,gravity=9.81)
    np.testing.assert_allclose(rest.acceleration_world,0.,atol=1e-10)
    with pytest.raises(ValueError):pressure_acceleration(owner,known_potential=np.zeros(count),
        known_normal_acceleration=np.zeros(count),body_speed_squared=None,density=1000.,gravity=9.81)
    with pytest.raises(ValueError):PressureSystem(fluid,np.ones(count,bool),length_scale=4.)


def test_constant_chi_cancels_uniform_drift_kinetic_pressure():
    U=np.array([.2,-.1,0.]);speed=U@U
    fluid=box_fluid(InitialState(linear_velocity=U))
    owner=PressureSystem(fluid,fluid.tags=='surface',length_scale=4.);n=len(owner.panels.areas)
    result=pressure_acceleration(owner,known_potential=np.full(n,-.5*speed),
        known_normal_acceleration=np.zeros(n),body_speed_squared=np.full(owner.body_mask.sum(),speed),density=1000.,gravity=9.81)
    np.testing.assert_allclose(result.acceleration_potential,-.5*speed,atol=1e-10)
    np.testing.assert_allclose(result.acceleration_world,0.,atol=1e-10)
