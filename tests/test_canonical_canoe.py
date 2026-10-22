"""Physical identities for the independent canonical canoe reference."""
import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from .experiments.free_surface.canonical_canoe import body_map, encode, decode
from .experiments.free_surface.moving_canoe import make_domain, initial_state, field, surface
from .experiments.free_surface.audit_canoe_energy import energy_terms, gravity_moments


def test_exponential_body_chart_rates_are_world_angular_velocity():
    rotvec = np.array([.2, -.3, .1])
    rate = np.array([.4, .1, -.2])
    rest = Rotation.from_rotvec([-.4, .2, .3])
    rotation = (Rotation.from_rotvec(rotvec)*rest).as_matrix()
    step = 1e-6
    plus = (Rotation.from_rotvec(rotvec+step*rate)*rest).as_matrix()
    minus = (Rotation.from_rotvec(rotvec-step*rate)*rest).as_matrix()
    spin = ((plus-minus)/(2*step))@rotation.T
    expected = np.array([spin[2, 1], spin[0, 2], spin[1, 0]])
    np.testing.assert_allclose(body_map(rotvec)[3:, 3:]@rate, expected, atol=1e-10)


@pytest.fixture(scope='module')
def canoe():
    domain = make_domain(12, gpu=False)
    state = initial_state(domain, scale=10.)
    with threadpool_limits(1):
        canonical = encode(domain, state)
        recovered, hamiltonian, qdot = decode(domain, canonical)
    return domain, state, canonical, recovered, hamiltonian, qdot


def test_canonical_momentum_round_trip_and_physical_energy(canoe):
    domain, state, _, recovered, hamiltonian, _ = canoe
    np.testing.assert_allclose(recovered, state, atol=3e-13, rtol=3e-12)
    # Independently solve the actual velocity boundary problem; this checks the
    # effective-mass Hamiltonian against physical fluid/body energy, not itself.
    with threadpool_limits(1):
        physical, _, _ = energy_terms(domain, state,
            gravity_moments(domain, domain.reference_points))
    np.testing.assert_allclose(hamiltonian, physical.sum(), atol=2e-12, rtol=2e-11)


def test_canonical_kinematics_match_actual_moving_surface(canoe):
    domain, state, canonical, _, _, qdot = canoe
    size = len(domain.free_nodes)
    with threadpool_limits(1):
        current = field(domain, state[:3], state[3:7], state[13:13+size],
                        state[13+size:], state[7:13])
        height_rate, _, _ = surface(domain, current, state[13+size:], state[7:13])
        # A momentum-direction FD tests dH/dp independently of the analytic
        # qdot formula, including body/free-surface cross coupling.
        direction = np.zeros_like(canonical)
        direction[len(canonical)//2:] = np.cos(np.arange(len(canonical)//2))
        step = 1e-5
        plus = decode(domain, canonical+step*direction)[1]
        minus = decode(domain, canonical-step*direction)[1]
    np.testing.assert_allclose(qdot[6:], height_rate, atol=3e-12, rtol=3e-11)
    np.testing.assert_allclose(body_map(canonical[3:6])@qdot[:6], state[7:13], atol=3e-13)
    np.testing.assert_allclose((plus-minus)/(2*step),
        qdot@direction[len(canonical)//2:], atol=3e-9, rtol=3e-8)
