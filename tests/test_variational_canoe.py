"""GPU checks of actual nonlinear variational forces, gauge and physical work."""
from importlib.util import find_spec
import numpy as np
import pytest
from threadpoolctl import threadpool_limits

if find_spec('cupy') is None:
    pytest.skip('optional CUDA environment required',allow_module_level=True)

from .experiments.free_surface.moving_canoe import make_domain,initial_state
from .experiments.free_surface.canonical_canoe import encode,decode,physical_rate
from .experiments.free_surface.variational_canoe import rhs
from .experiments.free_surface.audit_canoe_energy import energy_terms,gravity_moments


@pytest.fixture(scope='module')
def canoe():
    with threadpool_limits(1):
        domain=make_domain(12,True)
        physical=initial_state(domain,scale=10.)
        canonical=encode(domain,physical)
        derivative=rhs(domain,canonical)
        decoded_rate=physical_rate(domain,canonical,derivative)
    return domain,physical,canonical,derivative,decoded_rate


def test_complete_canoe_forces_match_fixed_momentum_hamiltonian_variations(canoe):
    domain,_,canonical,derivative,_=canoe
    size=len(canonical)//2
    directions=list(np.eye(size)[:6])
    xy=domain.reference_points[domain.free_nodes,:2]
    for values in (np.cos(.5*xy[:,0])*np.cos(xy[:,1]),np.sin(.7*xy[:,0]+.4*xy[:,1])):
        values-=values.mean()
        directions.append(np.r_[np.zeros(6),values/np.max(np.abs(values))])
    estimates=[]
    with threadpool_limits(1):
        for step in (2e-5,1e-5):
            row=[]
            for direction in directions:
                offset=np.r_[step*direction,np.zeros(size)]
                row.append((decode(domain,canonical+offset)[1]-decode(domain,canonical-offset)[1])/(2*step))
            estimates.append(np.array(row))
    independent=(4*estimates[1]-estimates[0])/3
    force=-derivative[size:]
    np.testing.assert_allclose(independent,[force@direction for direction in directions],atol=8e-8,rtol=3e-8)


def test_constant_potential_gauge_preserves_physical_canoe_rates(canoe):
    domain,physical,canonical,_,decoded_rate=canoe
    shifted=physical.copy()
    shifted[13+len(domain.free_nodes):]+=.3
    with threadpool_limits(1):
        encoded=encode(domain,shifted)
        shifted_rate=physical_rate(domain,encoded,rhs(domain,encoded))
        energy=decode(domain,encoded)[1]
        original_energy=decode(domain,canonical)[1]
    np.testing.assert_allclose(shifted_rate,decoded_rate,atol=2e-8,rtol=1e-7)
    np.testing.assert_allclose(energy,original_energy,atol=2e-10,rtol=0)


def test_decoded_nonlinear_motion_balances_independent_physical_energy(canoe):
    domain,physical,_,_,decoded_rate=canoe
    reference=gravity_moments(domain,domain.reference_points)
    step=1e-4
    with threadpool_limits(1):
        plus=energy_terms(domain,physical+step*decoded_rate,reference)
        minus=energy_terms(domain,physical-step*decoded_rate,reference)
    # Physical velocity-boundary solves and separate volume moments evaluate
    # this energy; algebraic Hamiltonian cancellation cannot pass the check.
    assert abs(np.sum(plus[0]-minus[0])/(2*step))<1e-6
    assert abs((plus[1]-minus[1])/(2*step))<1e-8
