"""Independent discrete-field checks on the actual reference-hull checkpoint."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.linalg import lu_factor, lu_solve

from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.bem import assemble_dirichlet_halfspace
from tests.reference_boat import simulation_inputs
from tests.experiments.free_surface.boat_radiation_time import Radiation, harmonic_maps, drive
from tests.experiments.free_surface.boat_spatial_work import source_matrix


@pytest.fixture(scope='module')
def system():
    pytest.importorskip('cupy')
    return Radiation(simulation_inputs(), 33, 32., 1, .5)


def cpu_extend(system, exterior):
    """Independent NumPy reference only; runtime owns a single GPU path."""
    full = np.zeros((system.size*system.size, exterior.shape[1]))
    full[system.exterior] = exterior
    full[system.hidden] = system.cp.asnumpy(system.fill_weights)@full[system.fill_boundary]
    return full.reshape(system.size, system.size, -1)


def test_hidden_fill_reproduces_harmonic_polynomials(system):
    x, y = np.meshgrid(system.origin[0]+np.arange(system.size)*system.spacing,
                        system.origin[1]+np.arange(system.size)*system.spacing)
    exact = np.stack((1+.2*x-.1*y, x*y, x*x-y*y), axis=-1).reshape(-1, 3)
    full = system.cp.asnumpy(system.extend(exact[system.exterior])).reshape(-1, 3)
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
    psi_rate = -GRAVITY*eta-system.cp.asnumpy(system.sponge)[:, None]*psi
    assert_allclose(system.cp.asnumpy(rates[1]), psi_rate, atol=2e-17, rtol=2e-15)
    p, n, _ = harmonic_maps(system.panels, system.size, system.spacing, system.origin)
    factor = lu_factor(source_matrix(system.panels))
    single = assemble_dirichlet_halfspace(system.panels).single_layer
    velocity = system.modes@system.cp.asnumpy(system.body_velocity)

    def static_moment(t, trace):
        spectrum = np.fft.fft2(cpu_extend(system, trace), axes=(0, 1)).reshape(-1, 2)
        speed, _ = drive(t, omega, ramp)
        density = lu_solve(factor, velocity*speed-(n@spectrum).real)
        potential = (p@spectrum).real+single@density
        return WATER_DENSITY*system.modes.T@(system.panels.areas[:, None]*potential)

    for step in [1e-4, 5e-5]:
        difference = -(static_moment(time+step, psi+step*psi_rate)
                       -static_moment(time-step, psi-step*psi_rate))/(2*step)
        assert_allclose(system.cp.asnumpy(actual), difference, rtol=2e-8, atol=2e-8)


def test_gpu_rhs_matches_independent_cpu_field_and_pressure(system):
    # Both driven columns, nonzero state and variable sponge inside the ramp.
    i = system.exterior
    eta = np.column_stack((.002*np.cos(i*.007), .001*np.sin(i*.011)))
    psi = np.column_stack((.03*np.sin(i*.013), .02*np.cos(i*.017)))
    time, omega, ramp = 7.3, 2*np.pi*.1, 20.
    actual = system.rhs(time, np.stack((eta, psi)), omega, ramp)
    p, n, wave = harmonic_maps(system.panels, system.size, system.spacing, system.origin)
    factor = lu_factor(source_matrix(system.panels))
    single = assemble_dirichlet_halfspace(system.panels).single_layer
    velocity = system.modes@system.cp.asnumpy(system.body_velocity)
    speed, acceleration = drive(time, omega, ramp)
    spectrum = np.fft.fft2(cpu_extend(system, psi), axes=(0, 1)).reshape(-1, 2)
    density = lu_solve(factor, velocity*speed-(n@spectrum).real)
    flux = np.fft.ifft2(spectrum.reshape(system.size, system.size, 2)*wave[:, :, None],
                         axes=(0, 1)).real.reshape(-1, 2)[system.exterior]
    flux += system.cp.asnumpy(system.flux_map)@density
    sponge = system.cp.asnumpy(system.sponge)[:, None]
    rates = np.stack((flux-sponge*eta, -GRAVITY*eta-sponge*psi))
    rate_spectrum = np.fft.fft2(cpu_extend(system, rates[1]), axes=(0, 1)).reshape(-1, 2)
    density_rate = lu_solve(factor, velocity*acceleration-(n@rate_spectrum).real)
    potential_rate = (p@rate_spectrum).real+single@density_rate
    weighted = system.panels.areas[:, None]
    force = -WATER_DENSITY*system.modes.T@(weighted*potential_rate)
    potential = (p@spectrum).real+single@density
    energy = .5*WATER_DENSITY*(system.spacing**2*np.sum(GRAVITY*eta**2+psi*flux, axis=0)
                              +np.sum(weighted*potential*velocity*speed, axis=0))
    power = -WATER_DENSITY*system.spacing**2*np.sum(sponge*(GRAVITY*eta**2+psi*flux), axis=0)
    for gpu, cpu in zip(actual, (rates, force, energy, power), strict=True):
        assert_allclose(system.cp.asnumpy(gpu), cpu, rtol=2e-11, atol=2e-11)
