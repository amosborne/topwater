"""Analytic passive causal impedance controls, with no BEM or fluid solve."""
import numpy as np
import pytest

from .capytaine_time_reference import radiation_kernel, radiation_force, reconstruct_coefficients


def coefficients():
    # Dense coupled passive memory, with a distinct noncommuting inertial map.
    lower = np.tril(np.ones((6, 6))*.12)+np.diag(np.arange(1., 7.)/4)
    coupling = lower@lower.T
    added = np.diag(np.arange(10., 16.))
    return 2., coupling, added


def test_known_causal_impedance_recovers_all_six_time_force_columns():
    alpha, coupling, added = coefficients()
    beta = .7
    errors = []
    for step in (.02, .01):
        times = np.linspace(0., 3., round(3/step)+1)
        kernel = np.exp(-alpha*times)[:, None, None]*coupling
        # Six separate input columns prevent compensation in one mixed twist.
        column_errors = []
        for column in range(6):
            axis = np.eye(6)[column]
            velocity = (1-np.exp(-beta*times))[:, None]*axis
            acceleration = (beta*np.exp(-beta*times))[:, None]*axis
            actual = radiation_force(times, velocity, acceleration, added, kernel,
                                     added_mass_infinite_depth=4.)
            memory = ((1-np.exp(-alpha*times))/alpha
                      -(np.exp(-beta*times)-np.exp(-alpha*times))/(alpha-beta))
            expected = -acceleration@added.T-memory[:, None]*(coupling@axis)
            np.testing.assert_array_equal(actual[0], expected[0])
            column_errors.append(np.max(abs(actual-expected)))
        errors.append(max(column_errors))
    assert errors[1] < .27*errors[0]
    assert errors[1] < 1e-4


def test_finite_band_cosine_to_causal_force_preserves_startup_and_cross_coupling():
    alpha, coupling, added = coefficients()
    times = np.linspace(0., 2., 401)
    omega = np.linspace(0., 200., 4001)
    damping = alpha/(alpha**2+omega[:, None, None]**2)*coupling
    kernel = radiation_kernel(omega, damping, times)
    beta = .7
    direction = np.array([1., -.3, .6, -.2, .8, -.5])
    velocity = (1-np.exp(-beta*times))[:, None]*direction
    acceleration = (beta*np.exp(-beta*times))[:, None]*direction
    actual = radiation_force(times, velocity, acceleration, added, kernel,
                             added_mass_infinite_depth=4.)
    memory = ((1-np.exp(-alpha*times))/alpha
              -(np.exp(-beta*times)-np.exp(-alpha*times))/(alpha-beta))
    expected = -acceleration@added.T-memory[:, None]*(coupling@direction)
    # Known omitted Lorentz tail bounds the kernel at t0. There is no hidden
    # taper or inferred A_inf; the force comparison includes finite-band error.
    tail = 1-2/np.pi*np.arctan(omega[-1]/alpha)
    np.testing.assert_allclose(kernel[0], (1-tail)*coupling, atol=2e-5, rtol=0.)
    assert np.max(abs(actual-expected)) < 2e-4


def test_exponential_memory_back_transform_has_correct_added_mass_sign():
    alpha, coupling, added = coefficients()
    times = np.linspace(0., 10., 10001)
    omega = np.array([.5, 2., 4., 6.])
    kernel = np.exp(-alpha*times)[:, None, None]*coupling
    recovered_a, recovered_b = reconstruct_coefficients(times, kernel, added, omega)
    expected_a = added-coupling/(alpha**2+omega[:, None, None]**2)
    expected_b = alpha*coupling/(alpha**2+omega[:, None, None]**2)
    np.testing.assert_allclose(recovered_a, expected_a, atol=3e-6, rtol=0.)
    np.testing.assert_allclose(recovered_b, expected_b, atol=6e-6, rtol=0.)


def test_causality_keeps_future_motion_out_of_earlier_force_and_requires_depth4():
    alpha, coupling, added = coefficients()
    times = np.linspace(0., 1., 101)
    kernel = np.exp(-alpha*times)[:, None, None]*coupling
    velocity = np.zeros((len(times), 6)); acceleration = velocity.copy()
    baseline = radiation_force(times, velocity, acceleration, added, kernel,
                               added_mass_infinite_depth=4.)
    velocity[60:, 4] = np.linspace(0., 1., len(times)-60)
    moved = radiation_force(times, velocity, acceleration, added, kernel,
                            added_mass_infinite_depth=4.)
    np.testing.assert_array_equal(moved[:60], baseline[:60])
    with pytest.raises(ValueError, match='depth4'):
        radiation_force(times, velocity, acceleration, added, kernel,
                         added_mass_infinite_depth=np.inf)
