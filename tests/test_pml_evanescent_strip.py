"""Evanescent launch/coordinate/reflection controls; PDE runs follow review."""
import numpy as np

from .experiments.free_surface.pml_evanescent_strip import (
    make_evanescent_strip,exact_evanescent_field,solve_evanescent_case)
from .experiments.free_surface.pml_strip import stretching
from .experiments.free_surface.sponge_strip import vertical_modes, DEPTH, GRAVITY


def test_evanescent_source_bc_and_complex_wall_phase_are_exact():
    strip=make_evanescent_strip(2.,np.pi/4)
    _,mu,_=vertical_modes(strip.omega,1)
    alpha=np.sqrt(mu[0]**2+strip.ky**2)
    np.testing.assert_allclose(-2*alpha*strip.xlevels[0],1.,atol=2e-15)
    points=np.column_stack((np.full(5,strip.xlevels[0]),np.linspace(-DEPTH,0,5)))
    value,gradient,reflection=exact_evanescent_field(strip,points,8.)
    vertical=np.cos(mu[0]*(points[:,1]+DEPTH))
    # left outward derivative is -u_x (s=1), Lambda_evan=-alpha.
    np.testing.assert_allclose(-gradient[:,0]+alpha*value,2*alpha*vertical,atol=2e-14,rtol=2e-13)
    wall=points.copy();wall[:,0]=strip.xlevels[-1]
    assert np.max(abs(exact_evanescent_field(strip,wall,8.)[1][:,0]))<1e-13
    top=points.copy();top[:,0]=np.linspace(0,strip.xlevels[-1],5);top[:,1]=0.
    v,g,_=exact_evanescent_field(strip,top,8.)
    np.testing.assert_allclose(g[:,1],strip.omega**2/GRAVITY*v,atol=2e-13,rtol=2e-12)
    amplitude=np.exp(-2*alpha*strip.xlevels[-1])
    np.testing.assert_allclose(abs(reflection),amplitude,atol=1e-15,rtol=2e-13)
    for strength in (0.,4.,8.):
        actual=exact_evanescent_field(strip,points,strength)[2]
        np.testing.assert_allclose(abs(actual),amplitude,atol=1e-15,rtol=2e-13)


def test_evanescent_fullfield_solves_stretched_elliptic_equation():
    strip=make_evanescent_strip(4.,np.pi/4)
    _,mu,_=vertical_modes(strip.omega,1)
    point=np.array([[-.1,-3.],[.3,-.7],[.9,-2.],[1.3,-.2]])
    value,_,_=exact_evanescent_field(strip,point,8.)
    errors=[]
    for step in (2e-4,1e-4):
        flux=[]
        for sign in (-1,1):
            moved=point.copy();moved[:,0]+=sign*step
            _,g,_=exact_evanescent_field(strip,moved,8.)
            s=stretching(moved[:,0],strip.xlevels[-1],8.,strip.omega)[1]
            flux.append(g[:,0]/s)
        s=stretching(point[:,0],strip.xlevels[-1],8.,strip.omega)[1]
        residual=(flux[1]-flux[0])/(2*step)-s*(mu[0]**2+strip.ky**2)*value
        errors.append(np.max(abs(residual)))
    assert errors[-1]<2e-6
    assert errors[-1]<.4*errors[0]+1e-10


def test_selected_evanescent_field_and_wall_return_match_without_power_normalization():
    # Two small PDE solves only after source review.
    strip=make_evanescent_strip(2.,0.,spacing=np.pi/10,layers=8)
    outgoing=solve_evanescent_case(strip,0.,right='outgoing')
    wall=solve_evanescent_case(strip,4.,right='wall')
    for case in (outgoing,wall):
        assert case['relative_system_residual']<1e-9
        assert case['absolute_complex_reflection_error']<.002
        assert case['spurious_propagating_mode_amplitude_at_left']<.002
        assert case['fullfield_errors']['full_strip_potential_relative_l2']<.002
        assert case['fullfield_errors']['quiet_physical_velocity_relative_l2']<.002
        assert case['incident_physical_power_flux']==0.
        assert case['incident_power_normalization_used'] is False


def test_real_stretch_changes_evanescent_distance_without_incident_power_normalization():
    compact = 1.5943951023931955
    physical = 7.877580409572782
    real_terminal = 2*physical/compact-1
    sigma_terminal = 8*physical/compact
    for omega in (2., 4., 6.):
        strip = make_evanescent_strip(omega, 0., compact, spacing=2*np.pi/9, layers=4)
        _, roots, _ = vertical_modes(omega, 1)
        alpha = roots[0]
        point = np.array([[strip.xlevels[0], -.4], [.8, -2.], [compact, -.1]])
        value, derivative, reflected = exact_evanescent_field(strip, point,
            sigma_terminal, real_terminal=real_terminal)
        expected = np.exp(-2*alpha*physical-1j*alpha*8*physical/omega)
        np.testing.assert_allclose(reflected, expected, atol=2e-16, rtol=3e-14)
        np.testing.assert_allclose(abs(reflected), np.exp(-2*alpha*physical), atol=2e-16, rtol=3e-14)
        # Source and wall conormal identities remain exact after real stretch.
        vertical = np.cos(roots[0]*(point[0, 1]+DEPTH))
        np.testing.assert_allclose(-derivative[0, 0]+alpha*value[0],
                                   2*alpha*vertical, atol=2e-14, rtol=2e-13)
        assert abs(derivative[-1, 0]) < 1e-13
