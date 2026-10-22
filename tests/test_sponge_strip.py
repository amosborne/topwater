"""Independent port/integration controls and opt-in-by-selection strip PDE checks.

The source is reviewed before any PDE checks are executed. These are contained
body-free diagnostics, not a canoe or outgoing-boundary qualification suite.
"""
import numpy as np
from scipy.integrate import quad

from .experiments.free_surface.sponge_strip import (
    DEPTH, make_strip, port, surface_matrices, vertical_modes,
    volume_operator, depth_mode_amplitude, solve_case)


def test_high_mode_port_projection_is_resolved_by_independent_oscillatory_integral():
    strip = make_strip(6., np.pi/4, width=.8)
    actual = port(strip, 'left', 128)
    checked = port(strip, 'left', 128, quadrature_multiplier=2)
    assert max(actual['quadrature_orders']) > 80
    np.testing.assert_allclose(actual['projection'], checked['projection'], atol=2e-13, rtol=2e-11)
    np.testing.assert_allclose(actual['matrix'], checked['matrix'], atol=2e-11, rtol=2e-10)
    # The port basis reproduces these physical polynomials exactly. Independent
    # QUADPACK cosine-weighted integration handles oscillation analytically; it
    # does not reuse the source's Gauss nodes or per-segment sampling.
    height = (strip.points[actual['nodes'], 1]+DEPTH)/DEPTH
    for degree in range(5):
        moments = actual['projection'].T@(height**degree)
        for index in (0, 7, 63, 127):
            mu = actual['mu'][index]
            exact = quad(lambda t: (t/DEPTH)**degree, 0., DEPTH,
                         weight='cos', wvar=mu, epsabs=2e-13)[0]
            np.testing.assert_allclose(moments[index+1], exact, atol=3e-12, rtol=2e-10)
        k = actual['k']
        exact = quad(lambda z: ((z+DEPTH)/DEPTH)**degree
                     *np.cosh(k*(z+DEPTH))/np.cosh(k*DEPTH), -DEPTH, 0., epsabs=2e-13)[0]
        np.testing.assert_allclose(moments[0], exact, atol=2e-13, rtol=2e-12)
    assert np.linalg.eigvalsh(actual['matrix'].real).max() < 1e-12
    assert np.linalg.eigvalsh(actual['matrix'].imag).min() > -1e-12


def test_physical_sponge_is_psd_bounded_and_independent_of_quadrature():
    strip = make_strip(4., width=1.5943951023931955)
    free, mass, q8 = surface_matrices(strip, 4., order=8)
    free12, mass12, q12 = surface_matrices(strip, 4., order=12)
    np.testing.assert_array_equal(free, free12)
    np.testing.assert_allclose(mass.toarray(), mass12.toarray(), atol=0., rtol=0.)
    np.testing.assert_allclose(q8.toarray(), q12.toarray(), atol=2e-14, rtol=2e-12)
    assert np.linalg.eigvalsh(q12.toarray()).min() > -1e-13
    assert np.linalg.eigvalsh((4*mass-q12).toarray()).min() > -1e-13
    assert np.max(abs(q12@(np.ones(len(free))))) > 0


def test_actual_2d_volume_metric_and_independent_depth_sampling():
    strip = make_strip(4., width=.8, spacing=.6)
    operator = volume_operator(strip)
    direction = np.array([.7, -.4]);potential = strip.points@direction
    area = DEPTH*(strip.xlevels[-1]-strip.xlevels[0])
    np.testing.assert_allclose(potential@(operator@potential), area*(direction@direction), rtol=2e-12)
    # ky^2 mass must be present for the oblique full elliptic equation.
    oblique = make_strip(4., angle=np.pi/4, width=.8, spacing=.6)
    ones = np.ones(len(oblique.points))
    np.testing.assert_allclose(ones@(volume_operator(oblique)@ones),
        oblique.ky**2*DEPTH*(oblique.xlevels[-1]-oblique.xlevels[0]), rtol=2e-12)
    stations = np.linspace(.15, .75, 4)
    values = depth_mode_amplitude(strip, strip.points[:, 1]**4, stations)
    k, _, norm = vertical_modes(strip.omega, 1)
    exact = quad(lambda z: z**4*np.cosh(k*(z+DEPTH))/np.cosh(k*DEPTH),
                 -DEPTH, 0., epsabs=1e-12)[0]/norm[0]
    np.testing.assert_allclose(values, exact, atol=2e-11, rtol=2e-12)


def test_transparent_ports_and_reflecting_wall_have_physical_phase_and_power():
    # This test performs TWO small fresh PDE solves. Run only after source
    # review; the first study still reports64/128 and independent grid controls.
    strip = make_strip(4., angle=np.pi/4, width=1.5943951023931955,
                       spacing=np.pi/10, layers=8)
    through = solve_case(strip, 0., right='outgoing', modes=64)
    wall = solve_case(strip, 0., right='wall', modes=64)
    for result in (through, wall):
        assert result['relative_system_residual'] < 1e-9
        assert result['relative_flux_balance_error'] < 2e-9
        assert result['depth_mode_q12_q16_maximum_relative_difference'] < 1e-10
    assert through['port_reflection_amplitude'] < .002
    assert abs(through['transmission_amplitude']-1) < .002
    assert abs(through['transmission_phase_error_radians']) < .002
    np.testing.assert_allclose(wall['port_reflection_amplitude'], 1., atol=2e-9, rtol=0.)
    measured = complex(*wall['port_reflection_at_sponge_start'])
    expected = np.exp(2j*strip.kappa*strip.xlevels[-1])
    assert abs(measured-expected) < .002
    for fit in wall['quiet_region_fits']:
        assert abs(complex(*fit['reflection'])-expected) < .002
