"""Independent discrete-field checks on the actual reference-hull checkpoint."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.linalg import lu_solve

from topwater.constants import GRAVITY, WATER_DENSITY
from tests.reference_boat import simulation_inputs
from tests.experiments.free_surface.boat_radiation_time import Radiation, harmonic_maps, drive


@pytest.fixture(scope='module')
def system():
    return Radiation(simulation_inputs(), 33, 32., 1, .5)


def test_hidden_fill_reproduces_harmonic_polynomials(system):
    x, y = np.meshgrid(system.origin[0]+np.arange(system.size)*system.spacing,
                        system.origin[1]+np.arange(system.size)*system.spacing)
    exact = np.stack((1+.2*x-.1*y, x*y, x*x-y*y), axis=-1).reshape(-1, 3)
    full = system.extend(exact[system.exterior]).reshape(-1, 3)
    assert_allclose(full[system.hidden], exact[system.hidden], rtol=2e-13, atol=2e-13)


def test_oblique_fourier_maps_match_analytic_subsurface_field(system):
    m = np.array([2., 3.])
    k = 2*np.pi*m/system.length
    x, y = np.meshgrid(system.origin[0]+np.arange(system.size)*system.spacing,
                        system.origin[1]+np.arange(system.size)*system.spacing)
    phase = k[0]*x+k[1]*y+.37
    trace = np.stack((np.cos(phase), np.sin(phase)), axis=-1)
    spectrum = np.fft.fft2(trace, axes=(0, 1)).reshape(-1, 2)
    points, normals = system.panels.centroids, system.panels.normals
    argument = points[:, :2]@k+.37
    scale = np.exp(np.linalg.norm(k)*points[:, 2])
    value = scale[:, None]*np.column_stack((np.cos(argument), np.sin(argument)))
    derivative = scale[:, None]*np.column_stack((-np.sin(argument), np.cos(argument)))
    expected_normal = (normals[:, :2]@k)[:, None]*derivative
    expected_normal += normals[:, 2, None]*np.linalg.norm(k)*value
    p, n, _ = harmonic_maps(system.panels, system.size, system.spacing, system.origin)
    assert_allclose((p@spectrum).real, value, rtol=2e-12, atol=2e-13)
    assert_allclose((n@spectrum).real, expected_normal, rtol=2e-12, atol=2e-13)


def test_pressure_is_time_derivative_of_static_body_potential(system):
    # Nonzero exterior state, in the ramp; variable sponge contributes to psi_t.
    indices = system.exterior
    psi = np.column_stack((.03*np.sin(indices*.013), .02*np.cos(indices*.017)))
    eta = np.column_stack((.002*np.cos(indices*.007), .001*np.sin(indices*.011)))
    time, omega, ramp = 7.3, 2*np.pi*.1, 20.
    rates, actual, _, _ = system.rhs(time, np.stack((eta, psi)), omega, ramp)
    psi_rate = -GRAVITY*eta-system.sponge[:, None]*psi
    assert_allclose(rates[1], psi_rate, atol=0, rtol=0)
    p, n, _ = harmonic_maps(system.panels, system.size, system.spacing, system.origin)

    def static_moment(t, trace):
        spectrum = np.fft.fft2(system.extend(trace), axes=(0, 1)).reshape(-1, 2)
        speed, _ = drive(t, omega, ramp)
        density = lu_solve(system.factor, system.velocity*speed-(n@spectrum).real)
        potential = (p@spectrum).real+system.single@density
        return WATER_DENSITY*system.modes.T@(system.panels.areas[:, None]*potential)

    for step in [1e-4, 5e-5]:
        difference = -(static_moment(time+step, psi+step*psi_rate)
                       -static_moment(time-step, psi-step*psi_rate))/(2*step)
        assert_allclose(actual, difference, rtol=2e-8, atol=2e-8)
