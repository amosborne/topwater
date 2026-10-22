"""Analytic and discrete controls for the contained matched-coordinate strip.

PDE tests execute only after source review. No canoe or nonlinear boundary is
validated by this module. Existing sponge-strip tests own the port controls.
"""
import numpy as np

from .experiments.free_surface.sponge_strip import make_strip, volume_operator, surface_matrices
from .experiments.free_surface.pml_strip import stretching, exact_field, pml_matrices, solve_pml_case


def test_coordinate_derivative_fullfield_bc_and_stretched_laplace_equation():
    strip = make_strip(4., np.pi/4, 1.5943951023931955)
    x = np.array([-.5, .2, .6, 1.3, 2.])
    coordinate, stretch, _ = stretching(x, strip.xlevels[-1], 8., strip.omega)
    for step in (2e-5, 1e-5):
        plus = stretching(x+step, strip.xlevels[-1], 8., strip.omega)[0]
        minus = stretching(x-step, strip.xlevels[-1], 8., strip.omega)[0]
        np.testing.assert_allclose((plus-minus)/(2*step), stretch, atol=2e-9, rtol=2e-9)
    point = np.column_stack((x[:-1], np.array([-3., -.8, -2., -.1])))
    value, _, reflection = exact_field(strip, point, 8.)
    errors = []
    for step in (2e-4, 1e-4):
        derivative = []
        for sign in (-1, 1):
            moved = point.copy();moved[:, 0] += sign*step
            _, gradient, _ = exact_field(strip, moved, 8.)
            s = stretching(moved[:, 0], strip.xlevels[-1], 8., strip.omega)[1]
            derivative.append(gradient[:, 0]/s)
        s = stretching(point[:, 0], strip.xlevels[-1], 8., strip.omega)[1]
        residual = (derivative[1]-derivative[0])/(2*step)+s*(strip.k**2-strip.ky**2)*value
        errors.append(np.max(abs(residual)))
    assert errors[-1] < 2e-6
    assert errors[-1] < .4*errors[0]+1e-10
    top = point.copy();top[:, 1] = 0.
    u, grad, _ = exact_field(strip, top, 8.)
    np.testing.assert_allclose(grad[:, 1], strip.omega**2/9.81*u, atol=1e-13, rtol=2e-12)
    bed = point.copy();bed[:, 1] = -4.
    np.testing.assert_array_equal(exact_field(strip, bed, 8.)[1][:, 1], np.zeros(len(bed)))
    wall = point.copy();wall[:, 0] = strip.xlevels[-1]
    assert np.max(abs(exact_field(strip, wall, 8.)[1][:, 0])) < 1e-13
    expected = np.exp(2j*strip.kappa*strip.xlevels[-1]
                       -strip.kappa*8*strip.xlevels[-1]/strip.omega)
    np.testing.assert_allclose(reflection, expected, atol=1e-16, rtol=2e-14)


def test_zero_stretch_recovers_unstretched_elliptic_operator_and_mass():
    strip = make_strip(4., np.pi/4, .8, spacing=.6)
    energy, free, ms = pml_matrices(strip, 0., volume_order=12)
    reference = volume_operator(strip)
    free0, mass, _ = surface_matrices(strip, 0.)
    np.testing.assert_array_equal(free, free0)
    np.testing.assert_allclose(energy.toarray(), reference.toarray(), atol=2e-13, rtol=2e-12)
    np.testing.assert_allclose(ms.toarray(), mass.toarray(), atol=2e-14, rtol=2e-12)
    # Nonzero stretch is COMPLEX SYMMETRIC, not a real-SPD physical-Q matrix.
    stretched, _, ms = pml_matrices(strip, 4., volume_order=12)
    assert np.linalg.norm((stretched-stretched.T).data) < 1e-12
    assert np.max(abs(stretched.imag.data)) > 0
    assert np.max(abs(ms.imag.data)) > 0


def test_scalar_augmented_pml_give_same_complex_reflection_and_analytic_field():
    # Two small PDE solves only AFTER root/source review. The eventual study
    # still needs coefficient-quadrature and grid controls at fixed settings.
    strip = make_strip(4., np.pi/4, 1.5943951023931955, spacing=np.pi/10, layers=8)
    augmented = solve_pml_case(strip, 4., volume_order=12, formulation='augmented')
    scalar = solve_pml_case(strip, 4., volume_order=12, formulation='scalar')
    for result in (augmented, scalar):
        assert result['relative_system_residual'] < 1e-9
        assert result['dynamic_eta_identity_relative_error'] < 1e-10
        assert result['discrete_modal_flux_balance_error'] < 1e-9
        assert result['independent_modal_flux_balance_error'] < 1e-7
        assert result['absolute_complex_reflection_error'] < .002
        assert result['fullfield_errors']['full_strip_potential_relative_l2'] < .002
        assert result['fullfield_errors']['quiet_physical_velocity_relative_l2'] < .002
    np.testing.assert_allclose(augmented['reflection_at_layer_start'], scalar['reflection_at_layer_start'],
                               atol=2e-10, rtol=2e-9)
