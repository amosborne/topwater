"""Independent physical invariants of the checkout-only numerical experiment."""

import numpy as np
import pytest
from threadpoolctl import threadpool_limits
from .experiments.free_surface.verify import run
from .experiments.free_surface.evolution import run as nonlinear_check


def test_closed_tank_energy_and_flux():
    with threadpool_limits(1):
        coarse = run((4, 2, 2), order=6)
        refined_quadrature = run((4, 2, 2), order=10)
    # Constant potential has no physical velocity; the independently known
    # volume fixes the energy of all three world-coordinate linear fields.
    for result in (coarse, refined_quadrature):
        assert result["constant_flux_max"] < 1e-11
        assert np.max(np.abs(result["net_top_flux"])) < 1e-11
        eigenvalues = result["lowest_eigenvalues"]
        assert abs(eigenvalues[0]) < 1e-11
        assert eigenvalues[1] == pytest.approx(0.25 * np.tanh(1.0), rel=0.05)
        assert np.min(eigenvalues[1:]) > 0
    assert (
        refined_quadrature["linear_energy_relative_error"]
        < coarse["linear_energy_relative_error"] / 4
    )


def test_nonlinear_finite_depth_wave_terms_refine():
    # Independent nonzero second-order standing-wave terms, not a polynomial
    # used to build the numerical basis. Whole-field quadrature includes edges;
    # the separate interior-node check cannot qualify the continuous field.
    with threadpool_limits(1):
        coarse = nonlinear_check((8, 4, 4), 0.02, 0.25)
        fine = nonlinear_check((16, 8, 8), 0.02, 0.25)
        finest = nonlinear_check((32, 16, 16), 0.02, 0.25)
    for metric in (
        "quadratic_kinematic_relative_rms",
        "quadratic_dynamic_relative_rms",
    ):
        assert fine[metric] < coarse[metric] / 3
        assert fine[metric] < 0.01
    assert np.max(np.abs(fine["integrated_elevation_rates"])) < 1e-11
    for metric in (
        "continuous_quadratic_kinematic_relative_rms",
        "continuous_quadratic_dynamic_relative_rms",
    ):
        assert fine[metric] < coarse[metric] / 2
        assert finest[metric] < fine[metric] / 2
        assert finest[metric] < 0.01
    assert np.max(np.abs(finest["integrated_elevation_rates"])) < 1e-11
