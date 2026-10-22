"""Instantaneous physical limits and nonlinear prescribed-reaction checks."""
from importlib.util import find_spec

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

if find_spec('cupy') is None:
    pytest.skip('optional CUDA environment required', allow_module_level=True)

from topwater.hydrostatics import hydrostatics
from .experiments.free_surface.canonical_canoe import encode
from .experiments.free_surface.canoe_radiation import sponge
from .experiments.free_surface.moving_canoe import make_domain, initial_state, field
from .experiments.free_surface.prescribed_canoe import reaction


@pytest.fixture(scope='module')
def canoe():
    with threadpool_limits(1):
        domain = make_domain(12, True)
        state = initial_state(domain, scale=10.)
        canonical = encode(domain, state)
    return domain, state, canonical[-len(domain.free_nodes):]


def test_all_six_instantaneous_reactions_are_buoyancy_minus_added_mass(canoe):
    domain, _, _ = canoe
    rest = initial_state(domain, False)
    size = len(domain.free_nodes)
    eta = pi = np.zeros(size)
    wet = hydrostatics(domain.body, domain.rest)
    buoyancy = np.r_[wet.force, wet.moment]
    with threadpool_limits(1):
        solved = field(domain, rest[:3], rest[3:7], eta, np.zeros(size), np.zeros(6))
        for acceleration in np.eye(6):
            actual = reaction(domain, rest[:13], eta, pi, acceleration)
            expected = buoyancy-solved['added_mass']@acceleration
            np.testing.assert_allclose(actual['hydrodynamic_load'], expected, atol=2e-8, rtol=2e-10)
            np.testing.assert_allclose(actual['eta_rate'], 0., atol=1e-11)


def test_nonlinear_reaction_and_fluid_momentum_derivatives_are_step_resolved(canoe):
    domain, state, pi = canoe
    size = len(domain.free_nodes)
    acceleration = np.array([.03, -.02, .05, .1, -.07, .04])
    damping, _ = sponge(domain.reference_points[domain.free_nodes], 4.)
    with threadpool_limits(1):
        coarse = reaction(domain, state[:13], state[13:13+size], pi, acceleration,
                          damping=damping, step=2e-4)
        fine = reaction(domain, state[:13], state[13:13+size], pi, acceleration,
                        damping=damping, step=1e-4)
    # Nonflat waves, all six body velocities/accelerations, moving contact and
    # the weak sponge participate. Pi_dot is independently recovered through
    # full physical-state encoding, not the local Ceta directional formula.
    np.testing.assert_allclose(fine['hydrodynamic_load'], coarse['hydrodynamic_load'], atol=2e-6, rtol=2e-7)
    assert np.max(np.abs(fine['fluid_momentum_pushforward_error'])) < 2e-7
    assert np.linalg.norm(fine['fluid_momentum_pushforward_error']) < (
        .4*np.linalg.norm(coarse['fluid_momentum_pushforward_error'])+1e-9)


def test_undamped_reaction_preserves_constant_potential_gauge(canoe):
    domain, state, pi = canoe
    size = len(domain.free_nodes)
    shifted = state.copy()
    shifted[13+size:] += .3
    acceleration = np.array([.03, -.02, .05, .1, -.07, .04])
    with threadpool_limits(1):
        shifted_pi = encode(domain, shifted)[-size:]
        original = reaction(domain, state[:13], state[13:13+size], pi, acceleration)
        actual = reaction(domain, shifted[:13], shifted[13:13+size], shifted_pi, acceleration)
    np.testing.assert_allclose(actual['hydrodynamic_load'], original['hydrodynamic_load'], atol=2e-6, rtol=2e-7)
    np.testing.assert_allclose(actual['eta_rate'], original['eta_rate'], atol=1e-10, rtol=1e-8)
    np.testing.assert_allclose(actual['psi_rate'], original['psi_rate'], atol=1e-8, rtol=1e-7)
