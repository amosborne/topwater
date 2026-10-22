"""Algebra controls for the offline actual-canoe surface-space diagnostic."""
import numpy as np
from scipy.sparse import csr_matrix

from tests.experiments.free_surface.p4_reduced_spectrum import projected_columns
from tests.experiments.free_surface.p4_surface_basis import rectangle_cosines


def test_restricted_schur_matches_stationary_full_field_energy():
    # SPD full energy, mixed boundary traces; compare to the stationary lifted
    # field's complete energy, independently of a Schur-matrix construction.
    generator = np.random.default_rng(5)
    values = generator.normal(size=(7, 7))
    full = values.T@values+np.eye(7)
    kii, kif, kff = full[:3, :3], full[:3, 3:], full[3:, 3:]
    basis = generator.normal(size=(4, 2))
    action, right, extension = projected_columns(csr_matrix(kif), csr_matrix(kff), basis,
        lambda rhs: np.linalg.solve(kii, rhs))
    lifted = np.vstack((-extension, basis))
    np.testing.assert_allclose(kii@extension, right, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(basis.T@action, lifted.T@full@lifted, rtol=1e-13, atol=1e-13)
    for column in range(2):
        individual, _, _ = projected_columns(csr_matrix(kif), csr_matrix(kff), basis[:, column:column+1],
            lambda rhs: np.linalg.solve(kii, rhs))
        np.testing.assert_allclose(individual[:, 0], action[:, column], rtol=1e-13, atol=1e-13)


def test_actual_rectangle_counts_nesting_and_constant_gauge_trace():
    xy = np.array([[-2*np.pi, -np.pi], [2*np.pi, np.pi], [.17, -.23], [-1., .5]])
    basis, pairs = rectangle_cosines(xy, .9, enrichment=2.)
    assert len(pairs) == 326
    np.testing.assert_array_equal(basis[:, 0], np.ones(len(xy)))
    for enrichment, count in ((1., 85), (1.5, 189), (2., 326)):
        small, modes = rectangle_cosines(xy, .9, enrichment=enrichment)
        assert len(modes) == count
        np.testing.assert_array_equal(modes, pairs[:count])
        np.testing.assert_array_equal(small, basis[:, :count])
