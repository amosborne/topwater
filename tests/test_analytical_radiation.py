"""Independent integration oracle for the piecewise-linear B/omega model."""
import numpy as np
import pytest
from scipy.integrate import quad_vec
from topwater.solve.native.radiation import radiation_kernel


@pytest.mark.parametrize('dt,duration', [(1e-8, 3e-8), (.01, .04), (32.5, 130.)])
def test_analytical_kernel_matches_adaptive_quadrature_on_irregular_intervals(dt, duration):
    omega=np.array([.1,.37,1.4,2.9,6.3,12.])
    damping=np.random.default_rng(42).normal(size=(len(omega),6,6))
    time,actual=radiation_kernel(omega,damping,dt_s=dt,duration_s=duration)
    expected=[]
    for t in time:
        integral=np.zeros((6,6))
        for left,right,y0,y1 in zip(omega[:-1],omega[1:],damping[:-1]/omega[:-1,None,None],damping[1:]/omega[1:,None,None]):
            value,_=quad_vec(lambda w: w*(y0+(y1-y0)*(w-left)/(right-left))*np.cos(w*t),
                             left,right,epsabs=1e-10,epsrel=1e-12)
            integral+=value
        expected.append(2/np.pi*integral)
    np.testing.assert_allclose(actual,expected,atol=2e-10,rtol=2e-11)


def test_constant_damping_has_closed_form_cosine_integral():
    omega=np.array([.1,1.3,7.])
    # Constant B/omega gives B proportional to omega, the exact interpolant.
    time,kernel=radiation_kernel(omega,omega[:,None,None]*np.eye(6),dt_s=.01,duration_s=.1)
    a,b=omega[[0,-1]]
    exact=np.empty(len(time));exact[0]=(b*b-a*a)/np.pi
    t=time[1:]
    exact[1:]=2/np.pi*((b*np.sin(b*t)-a*np.sin(a*t))/t+(np.cos(b*t)-np.cos(a*t))/t**2)
    np.testing.assert_allclose(kernel,exact[:,None,None]*np.eye(6),atol=2e-11,rtol=1e-12)
