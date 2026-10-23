"""Independent cut-domain polynomial moments and hidden-fill verification."""
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from tests.experiments.free_surface.compact_p1_trace import CompactP1Trace
from tests.experiments.free_surface.physical_p1_patch import PhysicalP1Patch


def box(x0, y0, x1, y1):
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])


def polynomial_moments(patch):
    p, a = patch.cut_triangles, patch.cut_areas
    # Integrate coordinate polynomials by their exact triangle formulas.
    first = np.sum(a[:, None]*p.mean(axis=1), axis=0)
    second = np.sum(a[:, None]/12*(np.sum(p*p, axis=1)+np.sum(p, axis=1)**2), axis=0)
    return a.sum(), first, second


def test_straight_cut_exact_area_polynomials_and_single_hat_mass():
    hx, hy = .31, .23
    x, y = np.array([-hx, 0, hx]), np.array([-hy, 0, hy])
    patch = PhysicalP1Patch(x, y, box(0, -2, 2, 2))
    area, first, second = polynomial_moments(patch)
    assert_allclose(area, 2*hx*hy, rtol=2e-15)
    assert_allclose(first, [-hy*hx**2, 0], atol=2e-17)
    assert_allclose(second, [2*hy*hx**3/3, 2*hx*hy**3/3], rtol=3e-15)
    assert_array_equal(patch.active_nodes, [4])
    # Original single hat has total squared integral hx*hy/2; inversion
    # symmetry makes a half-plane cut exactly half that value.
    assert_allclose(patch.mass.toarray(), [[hx*hy/4]], rtol=2e-15)
    assert_allclose(patch.physical_node_moments.sum(), area, rtol=2e-15)
    assert_array_equal(patch.faces, CompactP1Trace(x, y, np.zeros((3, 3))).faces)
    assert_allclose(patch.basis_values.sum(axis=2), 1., rtol=2e-15)
    assert_allclose(patch.cut_areas.sum()+patch.hidden_areas.sum(), 4*hx*hy, rtol=2e-15)
    assert_allclose(np.einsum('tvi,tij->tvj', patch.hidden_basis_values,
                             patch.nodes[patch.faces[patch.hidden_parent_id]]), patch.hidden_triangles, atol=5e-17)
    assert_allclose(np.einsum('tvi,tij->tvj', patch.basis_values,
                             patch.nodes[patch.faces[patch.parent_id]]), patch.cut_triangles, atol=5e-17)


def test_corner_cut_and_fully_hidden_rows():
    hx, hy = .31, .23
    x, y = [-hx, 0, hx], [-hy, 0, hy]
    patch = PhysicalP1Patch(x, y, box(0, 0, 2, 2))
    area, first, second = polynomial_moments(patch)
    assert_allclose(area, 3*hx*hy, rtol=2e-15)
    assert_allclose(first, [-hy*hx**2/2, -hx*hy**2/2], atol=2e-17)
    assert_allclose(second, [hx**3*hy, hx*hy**3], rtol=3e-15)
    # Two of six equal-area hat parents are removed by this quadrant.
    assert_allclose(patch.mass.toarray(), [[hx*hy/3]], rtol=2e-15)
    hidden = PhysicalP1Patch(x, y, box(-2, -2, 2, 2))
    assert hidden.mass.shape == (0, 0)
    assert_array_equal(hidden.covered_nodes, [4])
    assert hidden.E.shape == (9, 0)
    assert len(hidden.cut_areas) == 0
    assert_array_equal(hidden.expand(np.zeros(0)), np.zeros(9))
    outside = PhysicalP1Patch(x, y, box(2, 2, 3, 3))
    assert_allclose(outside.mass.toarray(), [[hx*hy/2]], rtol=2e-15)


def test_fill_constant_affine_condensation_and_active_basis_order():
    axis = np.linspace(-2, 2, 9)
    patch = PhysicalP1Patch(axis, axis, box(-1.1, -1.1, 1.1, 1.1))
    assert len(patch.covered_nodes) > 1
    assert np.all(patch.physical_node_moments[patch.covered_nodes] == 0)
    coefficient = 2+.3*patch.nodes[patch.active_nodes, 0]-.2*patch.nodes[patch.active_nodes, 1]
    expanded = patch.expand(coefficient)
    expected = 2+.3*patch.nodes[patch.covered_nodes, 0]-.2*patch.nodes[patch.covered_nodes, 1]
    assert_allclose(expanded[patch.covered_nodes], expected, rtol=3e-15)
    assert_allclose(patch.expand(np.ones(len(coefficient)))[patch.covered_nodes], 1., rtol=3e-15)
    perimeter = np.setdiff1d(np.arange(len(patch.nodes)), patch.interior_nodes)
    assert_array_equal(expanded[perimeter], 0)
    maps = np.vstack((np.ones(len(patch.nodes)), patch.nodes.T))
    assert_allclose(patch.condense_maps(maps) @ coefficient, maps @ expanded, rtol=3e-15)
    # Independent degree-two target quadrature of the interpolated field,
    # followed by the same check under an explicitly permuted active basis.
    field = np.einsum('tvi,ti->tv', patch.basis_values,
                     expanded[patch.faces[patch.parent_id]])
    exact_integral = np.sum(patch.cut_areas/12*(np.sum(field**2, axis=1)+np.sum(field, axis=1)**2))
    assert_allclose(coefficient @ (patch.mass @ coefficient), exact_integral, rtol=4e-15)
    permutation = np.arange(len(coefficient))[::-1]
    permuted_mass = patch.mass[permutation][:, permutation]
    assert_allclose(coefficient[permutation] @ (permuted_mass @ coefficient[permutation]), exact_integral, rtol=4e-15)
    assert np.linalg.eigvalsh(patch.mass.toarray())[0] > 0
