"""Checks of the complete nonlinear fluid/pressure/six-DOF stage."""
import numpy as np
import pytest
from threadpoolctl import threadpool_limits
from .experiments.free_surface.moving_canoe import make_domain,initial_state,rhs


@pytest.fixture(scope='module')
def domain():
    return make_domain(12,gpu=False)


def test_coupled_canoe_retains_independent_hydrostatic_equilibrium(domain):
    # Equilibrium came from independently clipped hydrostatics; the dynamic
    # stage must balance total Bernoulli pressure and body weight itself.
    diagnostics={}
    with threadpool_limits(1):
        value=rhs(domain,initial_state(domain,False),diagnostics=diagnostics)
    assert np.max(np.abs(value))<1e-10
    np.testing.assert_allclose(diagnostics['fluid_volume'],
        32*np.pi**2-domain.body.mass/1000,rtol=0,atol=1e-12)


def test_nonlinear_six_dof_pressure_derivative_is_step_resolved(domain):
    state=initial_state(domain)
    with threadpool_limits(1):
        coarse=rhs(domain,state,2e-4)
        fine=rhs(domain,state,1e-4)
    # Nonflat water, nonzero potential, six-axis body motion and moving contact
    # all participate in this derivative, including its current added mass.
    acceleration=fine[7:13]
    assert np.all(np.abs(acceleration)>1e-7)
    assert np.linalg.norm(coarse[7:13]-acceleration)/np.linalg.norm(acceleration)<1e-7
    assert np.all(np.isfinite(fine))
