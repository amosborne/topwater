"""Analytic modes and independent triangle-plane checks on rectangular graphs."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from tests.experiments.free_surface.graph_geometry_reference import tensor_fourier, structured_P1


def test_oblique_fourier_values_and_derivatives_at_irregular_tensor_queries():
    lengths, origin = np.array([7., 11.]), np.array([-.4, 1.3])
    nx, ny = 25, 19
    gx, gy = np.meshgrid(origin[0]+np.arange(nx)*lengths[0]/nx,
                         origin[1]+np.arange(ny)*lengths[1]/ny)
    k, ell = 2*np.pi*np.array([3., -2.])/lengths, 2*np.pi*np.array([-2., 1.])/lengths
    def exact(x, y):
        p, q = k[0]*x+k[1]*y+.31, ell[0]*x+ell[1]*y-.27
        height = .7+.12*np.cos(p)+.07*np.sin(q)
        gradient = -.12*np.sin(p)[..., None]*k+.07*np.cos(q)[..., None]*ell
        return height, gradient
    eta = exact(gx, gy)[0]
    x, y = np.array([-1.2, .3, 1.9, 7.8]), np.array([-.6, 2.1, 4.7])
    actual = tensor_fourier(eta, x, y, lengths, origin)
    expected = exact(*np.meshgrid(x, y))
    for a, b in zip(actual, expected, strict=True):
        assert_allclose(a, b, atol=3e-14, rtol=3e-13)


def test_affine_roof_and_quadratic_nodal_triangle_planes():
    nx, ny = np.array([-2., -.7, 1.2, 3.]), np.array([-.9, .4, 2.6])
    gx, gy = np.meshgrid(nx, ny)
    x, y = np.array([-2., -1.3, -.2, 2.1, 3.]), np.array([-.9, -.3, .8, 2.6])
    xx, yy = np.meshgrid(x, y)
    value, gradient = structured_P1(nx, ny, 1+.2*gx-.3*gy, x, y)
    assert_allclose(value, 1+.2*xx-.3*yy, atol=5e-16)
    assert_allclose(gradient, np.broadcast_to([.2, -.3], gradient.shape), atol=2e-16)
    heights = gx**2+.3*gx*gy-.4*gy**2
    value, gradient = structured_P1(nx, ny, heights, x, y)
    # Independently solve the physical affine plane of the containing triangle.
    for j, qy in enumerate(y):
        for i, qx in enumerate(x):
            a = min(np.searchsorted(nx, qx, side='right')-1, len(nx)-2)
            b = min(np.searchsorted(ny, qy, side='right')-1, len(ny)-2)
            u, v = (qx-nx[a])/(nx[a+1]-nx[a]), (qy-ny[b])/(ny[b+1]-ny[b])
            ids = [(b, a), (b+1, a+1), (b, a+1) if v <= u else (b+1, a)]
            matrix = np.array([[1., nx[c], ny[r]] for r, c in ids])
            plane = np.linalg.solve(matrix, [heights[r, c] for r, c in ids])
            assert_allclose(value[j, i], plane@[1., qx, qy], atol=3e-15)
            assert_allclose(gradient[j, i], plane[1:], atol=3e-15)
    with pytest.raises(ValueError, match='outside'):
        structured_P1(nx, ny, heights, [nx[0]-1e-12], [ny[0]])
