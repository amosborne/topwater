"""Actual canoe linear-limit checks against nonlinear variational stages."""
from importlib.util import find_spec

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

if find_spec('cupy') is None:
    pytest.skip('optional CUDA environment required',allow_module_level=True)

from .experiments.free_surface.moving_canoe import make_domain,initial_state,moved_orientation
from .experiments.free_surface.canonical_canoe import encode,physical_rate
from .experiments.free_surface.variational_canoe import rhs
from .experiments.free_surface.linear_canoe_radiation import solve


def test_six_mode_response_is_reciprocal_and_passive():
    with threadpool_limits(1):
        report,_=solve(make_domain(12,True),[2.,4.,6.])
    for case in report['cases']:
        assert case['relative_system_residual']<1e-10
        assert case['reciprocity_relative']<1e-10
        assert case['passivity_identity_relative']<1e-10
        assert case['minimum_damping_eigenvalue']>-1e-10


def test_frequency_equations_are_the_actual_nonlinear_fluid_first_variation():
    with threadpool_limits(1):
        domain=make_domain(12,True)
        size=len(domain.free_nodes)
        omega=4.
        _,solved=solve(domain,[omega],sponge_strength=0.)
        # All six coordinates contribute, including finite rotational rates.
        displacement=np.array([.4+.2j,-.3+.5j,.7-.4j,.1+.3j,-.2-.1j,.3+.2j])
        velocity=-1j*omega*displacement
        potential=solved['potential_per_velocity'][0]@velocity
        eta=solved['height_per_velocity'][0]@velocity
        expected=np.r_[(-1j*omega*eta).real,(-1j*omega*potential[domain.free_nodes]).real]
        errors=[]
        for amplitude in (1e-4,5e-5):
            derivatives=[]
            for sign in (1,-1):
                scale=sign*amplitude
                state=initial_state(domain,False)
                state[:3]+=scale*displacement[:3].real
                state[3:7]=moved_orientation(state[3:7],scale*displacement[3:].real,1.)
                state[7:13]=scale*velocity.real
                state[13:13+size]=scale*eta.real
                state[13+size:]=scale*potential[domain.free_nodes].real
                canonical=encode(domain,state)
                derivatives.append(physical_rate(domain,canonical,rhs(domain,canonical))[13:])
            actual=(derivatives[0]-derivatives[1])/(2*amplitude)
            errors.append(np.linalg.norm(actual-expected)/np.linalg.norm(expected))
    assert errors[1]<2e-5,errors
    assert errors[1]<.4*errors[0],errors
