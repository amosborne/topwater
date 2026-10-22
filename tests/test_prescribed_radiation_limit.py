"""Prescribed nonlinear reactions must recover finite-frequency force phase."""
from importlib.util import find_spec

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

if find_spec('cupy') is None:
    pytest.skip('optional CUDA environment required', allow_module_level=True)

from topwater.hydrostatics import hydrostatics
from .experiments.free_surface.canonical_canoe import encode, physical_rate
from .experiments.free_surface.linear_canoe_radiation import solve
from .experiments.free_surface.moving_canoe import make_domain, initial_state, moved_orientation
from .experiments.free_surface.prescribed_canoe import reaction
from .experiments.free_surface.variational_canoe import rhs


@pytest.fixture(scope='module')
def domain():
    return make_domain(12, True)


def test_prescribed_reaction_recovers_all_six_complex_radiation_forces(domain):
    omega = 4.
    size = len(domain.free_nodes)
    with threadpool_limits(1):
        report, solved = solve(domain, [omega], sponge_strength=4.)
        response = report['cases'][0]
        added_mass = np.asarray(response['added_mass'])
        damping = np.asarray(response['radiation_damping'])
        restoring = hydrostatics(domain.body, domain.rest).stiffness
        # Mixed phases expose sign/convention errors in both added mass and
        # damping, including all cross-mode force and moment responses.
        displacement = np.array([.4+.2j, -.3+.5j, .7-.4j,
                                 .1+.3j, -.2-.1j, .3+.2j])
        velocity = -1j*omega*displacement
        acceleration = -omega**2*displacement
        potential = solved['potential_per_velocity'][0]@velocity
        height = solved['height_per_velocity'][0]@velocity
        expected_force = ((omega**2*added_mass+1j*omega*damping-restoring)@displacement).real
        expected_fluid = np.r_[(-1j*omega*height).real,
                                (-1j*omega*potential[domain.free_nodes]).real]
        assert np.all(np.abs(expected_force) > .01)
        force_errors = []
        fluid_errors = []
        for amplitude in (1e-4, 5e-5):
            forces = []
            fluid_rates = []
            for sign in (1, -1):
                scale = sign*amplitude
                state = initial_state(domain, False)
                state[:3] += scale*displacement[:3].real
                state[3:7] = moved_orientation(state[3:7], scale*displacement[3:].real, 1.)
                state[7:13] = scale*velocity.real
                state[13:13+size] = scale*height.real
                state[13+size:] = scale*potential[domain.free_nodes].real
                pi = encode(domain, state)[-size:]
                actual = reaction(domain, state[:13], state[13:13+size], pi,
                    scale*acceleration.real, damping=solved['damping'])
                forces.append(actual['hydrodynamic_load'])
                fluid_rates.append(np.r_[actual['eta_rate'], actual['psi_rate']])
            centered_force = (forces[0]-forces[1])/(2*amplitude)
            centered_fluid = (fluid_rates[0]-fluid_rates[1])/(2*amplitude)
            force_errors.append(np.abs((centered_force-expected_force)/expected_force))
            fluid_errors.append(np.linalg.norm(centered_fluid-expected_fluid)
                                /np.linalg.norm(expected_fluid))
    # Componentwise force checks prevent a large heave force from hiding an
    # incorrect small rotational response. Centering removes the DC buoyancy.
    assert np.max(force_errors[1]) < 2e-5, force_errors
    assert np.all(force_errors[1] < .4*force_errors[0]), force_errors
    assert fluid_errors[1] < 2e-5, fluid_errors
    assert fluid_errors[1] < .4*fluid_errors[0], fluid_errors


def test_prescribing_actual_free_acceleration_requires_no_controller(domain):
    size = len(domain.free_nodes)
    state = initial_state(domain, scale=10.)
    with threadpool_limits(1):
        canonical = encode(domain, state)
        free = rhs(domain, canonical)
        decoded = physical_rate(domain, canonical, free)
        actual = reaction(domain, state[:13], state[13:13+size], canonical[-size:],
                          decoded[7:13])
    np.testing.assert_allclose(actual['controller_load'], 0., atol=2e-6)
    np.testing.assert_allclose(actual['eta_rate'], free[6:6+size], atol=2e-12, rtol=2e-10)
    np.testing.assert_allclose(actual['pi_rate'], free[-size:], atol=2e-11, rtol=2e-10)
    np.testing.assert_allclose(actual['psi_rate'], decoded[13+size:], atol=2e-8, rtol=2e-7)
