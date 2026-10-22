"""Analytic ODE controls for rejecting invalid numerical geometry trials."""
import numpy as np
import pytest
from scipy.integrate import RK23

from topwater.hydrostatics import FloodingError, LossOfContactError
from .experiments.free_surface.canoe_contact import ContactTopologyError
from .experiments.free_surface.geometry_rk23 import GeometryRK23
from .experiments.free_surface.p4_geometry import InvalidGeometryError


def test_positive_domain_decay_rejects_oversized_trial_and_recovers_exact_solution():
    holder, rejected_accepted_states = {}, []

    def decay(time, state):
        if state[0] <= 0:
            solver = holder['solver']
            rejected_accepted_states.append((solver.t, solver.y.copy(), solver.f.copy()))
            raise InvalidGeometryError('positive-domain trial left admissible geometry')
        return -state

    solver = GeometryRK23(decay, 0., [1.], 3., first_step=3., max_step=3.,
                          rtol=1e-8, atol=1e-11)
    holder['solver'] = solver
    while solver.status == 'running':
        solver.step()
        assert solver.y[0] > 0
    assert solver.status == 'finished'
    assert solver.invalid_geometry_rejections > 0
    assert solver.invalid_geometry_rejections == len(rejected_accepted_states)
    # First oversized stage is y(0)+1.5*f(0)=-.5. It must not advance
    # the accepted state/derivative, or contaminate the retry's stage buffer.
    time, state, derivative = rejected_accepted_states[0]
    assert time == 0.
    np.testing.assert_array_equal(state, [1.])
    np.testing.assert_array_equal(derivative, [-1.])
    assert solver.y[0] == pytest.approx(np.exp(-3.), rel=5e-8)


@pytest.mark.parametrize('error_type', [ValueError, RuntimeError, FloodingError,
                                      LossOfContactError, ContactTopologyError])
def test_unrelated_and_physical_failures_propagate_without_retry(error_type):
    def rate(time, state):
        if time > 0:
            raise error_type('physical or unrelated failure')
        return -state

    solver = GeometryRK23(rate, 0., [1.], 1., first_step=1.)
    original_fun = solver.fun
    with pytest.raises(error_type, match='physical or unrelated failure'):
        solver.step()
    assert solver.fun is original_fun
    assert solver.invalid_geometry_rejections == 0
    assert solver.t == 0.
    np.testing.assert_array_equal(solver.y, [1.])
    np.testing.assert_array_equal(solver.f, [-1.])


def test_no_invalid_path_is_bitwise_plain_rk23_including_dense_output():
    def oscillator(time, state):
        return np.array([state[1], -state[0]])

    options = dict(first_step=.3, max_step=.7, rtol=1e-8, atol=1e-11)
    plain = RK23(oscillator, 0., [1., 0.], 2., **options)
    guarded = GeometryRK23(oscillator, 0., [1., 0.], 2., **options)
    while plain.status == 'running':
        assert plain.step() == guarded.step()
        assert plain.status == guarded.status
        assert plain.t == guarded.t
        assert plain.h_abs == guarded.h_abs
        assert plain.nfev == guarded.nfev
        np.testing.assert_array_equal(plain.y, guarded.y)
        np.testing.assert_array_equal(plain.f, guarded.f)
        np.testing.assert_array_equal(plain.K, guarded.K)
        midpoint = (plain.t_old+plain.t)/2
        np.testing.assert_array_equal(plain.dense_output()(midpoint), guarded.dense_output()(midpoint))
    assert guarded.invalid_geometry_rejections == 0


def test_persistent_invalid_trial_stops_at_representable_minimum_without_advancing():
    def rate(time, state):
        if time != 1.:
            raise InvalidGeometryError('no admissible positive trial')
        return -state

    minimum = 10*(np.nextafter(1., np.inf)-1.)
    solver = GeometryRK23(rate, 1., [1.], 2., first_step=minimum)
    message = solver.step()
    assert solver.status == 'failed'
    assert 'minimum representable RK23 step' in message
    assert solver.invalid_geometry_rejections == 1
    assert solver.t == 1.
    np.testing.assert_array_equal(solver.y, [1.])
    np.testing.assert_array_equal(solver.f, [-1.])


def test_geometry_backoff_uses_failed_stage_after_internal_accuracy_rejection():
    times = []

    def decay(time, state):
        times.append(time)
        # Initial f plus three stages complete the first oversized trial.
        # Its normal accuracy rejection shrinks h internally; inject a geometry
        # failure at the FIRST stage of the smaller trial, not at the old h.
        if len(times) == 5:
            raise InvalidGeometryError('injected failure after accuracy rejection')
        return -state

    solver = GeometryRK23(decay, 0., [1.], 1., first_step=.8, rtol=1e-8, atol=1e-11)
    original_fun = solver.fun
    solver.step()
    assert solver.fun is original_fun
    assert solver.invalid_geometry_rejections == 1
    assert times[4] < times[1]  # normal rejection already shrank the trial
    assert times[5] <= .5*times[4]  # retry cannot use half the stale outer h
    while solver.status == 'running':
        solver.step()
    assert solver.status == 'finished'
    assert solver.y[0] == pytest.approx(np.exp(-1.), rel=5e-8)
