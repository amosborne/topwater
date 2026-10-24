import numpy as np
from .dynamics import ramp,integrate,radiation_kernel


def test_zero_force_stays_at_rest_with_coupling():
    t=np.linspace(0,1,101);M=np.eye(6);M[0,2]=M[2,0]=.2
    q,v,a=integrate(t,M,np.eye(6),np.eye(6),np.exp(-t[:,None,None])*np.eye(6),np.zeros((len(t),6)))
    assert not np.any(q) and not np.any(v) and not np.any(a)


def test_polynomial_acceleration_refines_at_second_order():
    errors=[]
    for N in (101,201):
        t=np.linspace(0,1,N);F=t[:,None]**2*np.arange(1,7)[None,:]
        q,v,a=integrate(t,np.eye(6),np.zeros((6,6)),np.zeros((6,6)),np.zeros((N,6,6)),F)
        np.testing.assert_allclose(a,F,atol=1e-14)
        errors.append(np.max(abs(q-t[:,None]**4*np.arange(1,7)[None,:]/12)))
    assert 3.9<errors[0]/errors[1]<4.1


def test_kernel_matches_analytic_exponential_damping_transform():
    w=np.linspace(0,20,10001);t=np.linspace(0,2,21)
    B=np.exp(-w[:,None,None])*np.eye(6)
    K=radiation_kernel(w,B,t)
    np.testing.assert_allclose(K,2/np.pi/(1+t[:,None,None]**2)*np.eye(6),atol=3e-7)


def test_ramp_has_exact_endpoints_and_zero_start_slope():
    np.testing.assert_array_equal(ramp(np.array([-1.,0.,10.,11.]),10),[0,0,1,1])
    assert ramp(np.array([1e-5]),10)[0]<1e-16
