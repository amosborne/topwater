"""Periodic reference math and readiness controls; no BEM or runtime solve."""
import json
from pathlib import Path
import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from topwater import InitialState
from topwater.frames import rotation_matrix
from .mild import histories,spatial_envelope,initial_body,informative_signal,load_case


def test_periodic_sign_and_exact_derivatives():
    p=np.zeros((4,1,6),complex);p[:,:,0]=1+2j
    w=np.array([2.]);t=np.array([0.,.2,.7])
    q,v,a=histories(p,w,t)
    np.testing.assert_allclose(q[:,:,0],np.broadcast_to(np.cos(2*t)+2*np.sin(2*t),(4,3)))
    np.testing.assert_allclose(v[:,:,0],np.broadcast_to(-2*np.sin(2*t)+4*np.cos(2*t),(4,3)))
    np.testing.assert_allclose(a,-4*q)


def test_mesh_envelope_cannot_hide_cross_component_cancellation():
    p=np.zeros((4,2,6),complex)
    p[1,:,0]=[.0006,-.0006];p[2,:,0]=[.0008,-.0008];p[3,:,0]=[.0009,-.0009]
    r=spatial_envelope(p,np.array([1.,2.]),.005,np.ones(6)*.01,.1,np.ones(6)*1e-10)
    assert r['spatial'][0]==pytest.approx(.0004)
    assert r['sampling'][0]>0
    assert abs(histories(p,np.array([1.,2.]),np.array([0.]))[0][-1,0,0])<1e-15


def test_envelope_sum_and_sampling_share_reference_budget():
    p=np.zeros((4,2,6),complex)
    p[1,:,0]=[0,.0006];p[2,:,0]=[.0006,.0006];p[3,:,0]=[.0006,.0012]
    r=spatial_envelope(p,np.array([1.,2.]),.005,np.ones(6)*.01,.1,np.ones(6)*1e-10)
    assert r['spatial'][0]==pytest.approx(.0012)
    assert not r['gate']['passed']


def test_body_rate_is_derivative_of_exact_rotation_map():
    mean=InitialState(orientation=Rotation.from_euler('xyz',[.1,-.2,.3]).as_quat(scalar_first=True))
    q=np.array([.01,-.02,.03,.04,-.03,.02]);v=np.array([.1,.2,.3,.2,-.1,.4])
    state=initial_body(mean,q,v);step=1e-6
    plus=rotation_matrix(initial_body(mean,q+step*v,v).orientation)
    minus=rotation_matrix(initial_body(mean,q-step*v,v).orientation)
    expected=Rotation.from_matrix(minus.T@plus).as_rotvec()/(2*step)
    np.testing.assert_allclose(state.angular_velocity,expected,rtol=1e-8,atol=1e-10)
    np.testing.assert_array_equal(state.linear_velocity,v[:3])
    np.testing.assert_array_equal(state.position,mean.position+q[:3])


def test_near_tolerance_signal_does_not_pass_on_uncertainty():
    motion=np.zeros((2,6));motion[-1,[0,3]]=.011
    limits=np.ones(6)*.01
    assert informative_signal(motion,np.zeros(6),limits)[2]
    assert not informative_signal(motion,np.ones(6)*.001,limits)[2]


@pytest.mark.parametrize('change',['rest_start','rotation_chart','missing_level','relaxed_budget'])
def test_case_changes_fail_closed(tmp_path,change):
    c=json.loads(Path('validation/mild-periodic-case.json').read_text())
    if change=='rest_start':c['initialization']['rest_start_reference']=True
    if change=='rotation_chart':c['comparison']['rotation_error_chart']='Euler yaw pitch roll'
    if change=='missing_level':c['raw_reference']['required_levels']=[1,2,3]
    if change=='relaxed_budget':c['comparison']['reference_budget_fraction']=1.
    path=tmp_path/'case.json';path.write_text(json.dumps(c))
    with pytest.raises(ValueError):load_case(path)
