"""Independent affine moments and divergence tests on actual hat parents."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from tests.experiments.free_surface.cut_hat_geometry import cut_hat


def rectangle(x0, y0, x1, y1):
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])


def check_complement(result, hx, hy):
    inside, outside = result['interior'], result['exterior']
    assert np.all(inside['areas'] > 0) and np.all(outside['areas'] > 0)
    area = inside['areas'].sum()+outside['areas'].sum()
    mass = inside['hat_moments'].sum()+outside['hat_moments'].sum()
    assert_allclose(area, 3*hx*hy, rtol=3e-14)
    assert_allclose(mass, hx*hy, rtol=3e-14)
    for parent in range(6):
        original = result['parents'][parent]
        expected_area = abs(np.linalg.det(original[1:]-original[:1]))/2
        pieces = [(p, p['parent_id'] == parent) for p in (inside, outside)]
        actual_area = sum(p['areas'][mask].sum() for p, mask in pieces)
        actual_first = sum((p['areas'][mask, None]*p['triangles'][mask].mean(axis=1)).sum(axis=0)
                           for p, mask in pieces)
        assert_allclose(actual_area, expected_area, rtol=3e-14)
        assert_allclose(actual_first, expected_area*original.mean(axis=0), atol=2e-17)
    # Independent integration-by-parts check for the cut hat itself. Outer
    # patch boundary contributes zero because the complete hat trace vanishes.
    gradient_integral = (outside['areas'][:, None]*outside['gradients']).sum(axis=0)
    line = result['waterline']
    lengths = np.linalg.norm(line['segments'][:, 1]-line['segments'][:, 0], axis=1)
    boundary_integral = (lengths[:, None]*line['values'].mean(axis=1)[:, None]
                         *line['hole_outward_normals']).sum(axis=0)
    assert_allclose(gradient_integral+boundary_integral, 0., atol=2e-15)


def test_straight_and_corner_cut_affine_hat_moments():
    hx, hy = .31, .23
    straight = cut_hat([0, 0], hx, hy, rectangle(0, -2, 2, 2))
    corner = cut_hat([0, 0], hx, hy, rectangle(0, 0, 2, 2))
    check_complement(straight, hx, hy)
    check_complement(corner, hx, hy)
    assert_allclose(straight['physical_hat_mass'], hx*hy/2, rtol=2e-15)
    assert_allclose(corner['interior']['hat_moments'].sum(), hx*hy/3, rtol=2e-15)
    assert_allclose(straight['waterline']['hole_outward_normals'],
                    np.tile([-1., 0.], (len(straight['waterline']['segments']), 1)))


def test_small_hole_integral_and_translation():
    hx, hy, a = .31, .23, .25
    footprint = rectangle(-a*hx, -a*hy, a*hx, a*hy)
    result = cut_hat([0, 0], hx, hy, footprint)
    check_complement(result, hx, hy)
    # Integrate the original piecewise-affine hat over four rectangular
    # quadrants: same-sign quadrants use 1-max(u,v), opposite use 1-u-v.
    expected = hx*hy*(4*a*a-10*a**3/3)
    assert_allclose(result['interior']['hat_moments'].sum(), expected, rtol=3e-15)
    # A separate generic cut avoids making float-translation parity depend
    # on a symmetric footprint corner lying exactly on a parent diagonal.
    footprint = rectangle(-.21*hx, -.17*hy, .29*hx, .36*hy)
    result = cut_hat([0, 0], hx, hy, footprint)
    shift = np.array([.137, -.081])
    moved = cut_hat(shift, hx, hy, footprint+shift)
    assert_allclose(moved['physical_hat_mass'], result['physical_hat_mass'], rtol=3e-15)
    for kind in ('interior', 'exterior'):
        assert_allclose(moved[kind]['triangles']-shift, result[kind]['triangles'], atol=8e-17)
        assert_allclose(moved[kind]['values'], result[kind]['values'], atol=8e-16)


def test_fully_hidden_and_exterior_rows_and_rejected_geometry():
    hx, hy = .31, .23
    hidden = cut_hat([0, 0], hx, hy, rectangle(-2, -2, 2, 2))
    exterior = cut_hat([0, 0], hx, hy, rectangle(2, 2, 3, 3))
    for result in (hidden, exterior):
        check_complement(result, hx, hy)
        assert len(result['waterline']['segments']) == 0
    assert hidden['physical_hat_mass'] == 0
    assert_allclose(exterior['physical_hat_mass'], hx*hy, rtol=2e-15)
    with pytest.raises(ValueError, match='convex'):
        cut_hat([0, 0], hx, hy, rectangle(-2, -2, 2, 2)[::-1])
    with pytest.raises(ValueError):
        cut_hat([0, 0], 0, hy, rectangle(-2, -2, 2, 2))
