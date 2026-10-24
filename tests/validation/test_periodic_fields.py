import numpy as np
from .periodic_fields import incident,boundary_limits,envelope


def test_incident_elevation_phase_and_vertical_decay():
    w=2.;g=9.81;zeta=.02j;p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,0.,-1.]])
    phi=incident(p,w,0.,g,zeta)
    np.testing.assert_allclose(1j*w*phi/g,zeta*np.exp(w*w/g*(p[:,2]+1j*p[:,0])))
    assert abs(phi[2])<abs(phi[0])


def test_exterior_boundary_limit_preserves_complex_phase_and_reports_curvature():
    d=np.array([.02,.01,.005,.0025]);a=1+2j;b=-3+4j
    values=(a+b*d+5*d*d)[None,:,None]
    limit=boundary_limits(values)
    assert limit.shape==(1,3,1)
    np.testing.assert_allclose(limit[0,:,0],a-10*d[1:]**2)
    assert abs(limit[0,2,0]-a)<abs(limit[0,1,0]-a)


def test_field_mesh_envelope_cannot_cancel_wave_errors():
    v=np.zeros((4,2,1),complex);v[1,:,0]=[1,-1];v[2,:,0]=[1.5,-1.5];v[3,:,0]=[1.75,-1.75]
    u,trend=envelope(v,1e-9)
    np.testing.assert_allclose(u,[1.]);assert trend
