"""Independent physical energy budget under nonlinear prescribed motion."""
from importlib.util import find_spec
import json

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

if find_spec('cupy') is None:
    pytest.skip('optional CUDA environment required', allow_module_level=True)

from .experiments.free_surface.audit_canoe_energy import energy_terms, gravity_moments
from .experiments.free_surface.canonical_canoe import encode
from .experiments.free_surface.canoe_radiation import sponge
from .experiments.free_surface.moving_canoe import make_domain, initial_state
from .experiments.free_surface.prescribed_canoe import reaction
from .experiments.free_surface.variational_canoe import rhs


@pytest.fixture(scope='module')
def canoe():
    with threadpool_limits(1):
        domain = make_domain(12, True)
        state = initial_state(domain, scale=10.)
        canonical = encode(domain, state)
        free = rhs(domain, canonical)
    return domain, state, canonical, free


@pytest.mark.parametrize('damped', [False, True], ids=['undamped', 'weak-far-sponge'])
def test_controller_and_declared_sponge_balance_independent_physical_work(canoe, damped):
    domain, state, canonical, free = canoe
    size = len(domain.free_nodes)
    damping = sponge(domain.reference_points[domain.free_nodes], 4.)[0] if damped else None
    acceleration = np.array([.03, -.02, .05, .1, -.07, .04])
    reference = gravity_moments(domain, domain.reference_points)
    with threadpool_limits(1):
        prescribed = reaction(domain, state[:13], state[13:13+size], canonical[-size:],
                              acceleration, damping=damping)
        eta_source = prescribed['eta_rate']-free[6:6+size]
        pi_source = prescribed['pi_rate']-free[-size:]
        controller_power = state[7:13]@prescribed['controller_load']
        # Hq=-pdot_free and Hp=qdot_free. Sources affect eta/Pi and body
        # momentum; this prediction is compared with physical energy below.
        sponge_power = -free[-size:]@eta_source+free[6:6+size]@pi_source
        expected = controller_power+sponge_power
        measured = []
        for step in (2e-4, 1e-4):
            energies = []
            for sign in (1, -1):
                moved = prescribed['physical_state']+sign*step*prescribed['physical_rate']
                moved[3:7] /= np.linalg.norm(moved[3:7])
                # Fresh imposed-velocity fluid solves and exact volume gravity
                # moments are independent of the canonical work contraction.
                energies.append(energy_terms(domain, moved, reference)[0])
            measured.append(float(np.sum(energies[0]-energies[1])/(2*step)))
    errors = np.abs(np.asarray(measured)-expected)
    # Two parts per million of physical stage power, or excitation energy over
    # one second, rather than the multi-megajoule constant tank gravity datum.
    power_scale = max(abs(controller_power)+abs(sponge_power), abs(prescribed['hamiltonian'])/1.)
    tolerance = 2e-6*power_scale
    print(json.dumps(dict(damped=damped, controller_power_watts=float(controller_power),
        sponge_power_watts=float(sponge_power), expected_power_watts=float(expected),
        physical_energy_rates_watts=measured, absolute_errors_watts=errors.tolist(),
        tolerance_watts=float(tolerance))))
    assert errors[1] < tolerance
    # Resolve the derivative under step refinement, allowing the already-small
    # residual to reach floating-point/geometry-difference noise.
    assert abs(measured[1]-measured[0]) < tolerance/2
    assert errors[1] < .4*errors[0]+.1*tolerance
