"""Current exterior stiffness is checked against continuum wave dispersion."""
import numpy as np
from numpy.testing import assert_allclose
from .experiments.free_surface.exterior_dispersion import symbol,uniform_symbol,measure


def test_actual_periodic_tetra_assembly_matches_independent_stencil():
    for angle in (0.,np.pi/2,np.pi/4,-np.pi/4):
        for nz in (4,8,16):
            dx,dy=.52,.26
            kx,ky=2.7*np.cos(angle),2.7*np.sin(angle)
            actual=symbol(dx,dy,np.linspace(-4.,0.,nz+1),kx,ky)
            assert_allclose(actual,uniform_symbol(dx,dy,nz,kx,ky),rtol=2e-13,atol=2e-13)
    assert abs(symbol(.4,.2,np.linspace(-4.,0.,9),0.,0.))<1e-12


def test_periodic_wave_frequency_converges_to_continuum_in_tested_directions():
    for angle in (0.,np.pi/2,np.pi/4,-np.pi/4):
        errors=[abs(measure(n,2.,angle)['relative_frequency_error']) for n in (24,48,96)]
        # The exact finite-depth continuum solution supplies the reference;
        # three independent h levels expose the expected second-order error.
        assert .20<errors[1]/errors[0]<.30
        assert .20<errors[2]/errors[1]<.30


def test_horizontal_and_vertical_resolution_are_separate_error_sources():
    for angle in (0.,np.pi/2,np.pi/4,-np.pi/4):
        base=measure(48,6.,angle)['relative_frequency_error']
        horizontal=measure(48,6.,angle,horizontal_refinement=2)['relative_frequency_error']
        vertical=measure(48,6.,angle,vertical_refinement=2)['relative_frequency_error']
        both=measure(48,6.,angle,horizontal_refinement=2,vertical_refinement=2)['relative_frequency_error']
        assert 0<both<min(horizontal,vertical)<base
