"""Independent exterior target integration and signed-lag indexing controls."""
from types import SimpleNamespace
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.testing import assert_allclose
import pytest
from tests.laplace_reference import triangle_samples
from tests.experiments.free_surface.cut_hat_geometry import cut_hat
from tests.experiments.free_surface.plane_p1_field import PlaneP1Field
from tests.experiments.free_surface.physical_p1_patch import PhysicalP1Patch
from tests.experiments.free_surface.physical_weak_flux import PhysicalWeakFlux


def box(a, b, c, d):
    return np.array([[a, b], [c, b], [c, d], [a, d]])


def reference_kernel(shape, hx, hy, order=8):
    """Independent CPU target work; no GPU lag shader/FFT implementation."""
    ny, nx = shape
    result = np.empty((2*ny-1, 2*nx-1))
    geometry = cut_hat([0., 0.], hx, hy, box(10, 10, 11, 11))
    target = np.concatenate((geometry['parents'], np.zeros((6, 3, 1))), axis=2)
    points, weights = triangle_samples(target, order)
    areas = np.full(6, hx*hy/2)
    for dy in range(-ny+1, ny):
        for dx in range(-nx+1, nx):
            ix = np.arange(min(-1, dx-1), max(1, dx+1)+1)
            iy = np.arange(min(-1, dy-1), max(1, dy+1)+1)
            values = np.zeros((len(iy), len(ix), 1))
            values[np.flatnonzero(iy == dy)[0], np.flatnonzero(ix == dx)[0], 0] = 1
            field = PlaneP1Field(ix*hx, iy*hy, values)
            p = points[:, :, :2].reshape(-1, 2)
            A = field.vector_potential(p, field.point_owners(p))[:, :, 0].reshape(6, -1, 2)
            result[-dy+ny-1, -dx+nx-1] = 2*np.einsum('t,q,ti,tqi->',
                areas, weights, geometry['parent_gradients'], A)
    return SimpleNamespace(kernel=result, spacing=(hx, hy), shape=shape)


def exterior_row(field, center, hx, hy, footprint, order):
    """Direct physical area+boundary work, not full-minus-hidden correction."""
    geometry = cut_hat(center, hx, hy, footprint)
    pieces = geometry['exterior']
    triangles = np.concatenate((pieces['triangles'], np.zeros((*pieces['triangles'].shape[:2], 1))), axis=2)
    points, weights = triangle_samples(triangles, order)
    p = points[:, :, :2].reshape(-1, 2)
    A = field.vector_potential(p, field.point_owners(p))[:, :, 0].reshape(len(triangles), -1, 2)
    area = 2*np.einsum('t,q,ti,tqi->', pieces['areas'], weights, pieces['gradients'], A)
    line = geometry['waterline']
    g, w = leggauss(order)
    t, w = (g+1)/2, w/2
    boundary = 0.
    for segment, normal, values in zip(line['segments'], line['hole_outward_normals'], line['values']):
        p = segment[:1]+t[:, None]*(segment[1:]-segment[:1])
        A = field.vector_potential(p, field.point_owners(p))[:, :, 0]
        boundary += 2*np.linalg.norm(segment[1]-segment[0])*np.einsum('q,q,qi,i->',
            w, (1-t)*values[0]+t*values[1], A, normal)
    return area+boundary


def test_physical_rows_against_independent_exterior_work_and_scaling():
    hx, hy = .3, .23
    x, y = np.arange(-2, 3)*hx, np.arange(-2, 2)*hy
    footprint = box(.11, .07, 2, 2)
    patch = PhysicalP1Patch(x, y, footprint)
    owner = reference_kernel((2, 3), hx, hy)
    operator = PhysicalWeakFlux(patch, footprint, owner, order=8)
    coefficients = np.array([.2, -.1, .7, -.3, .4, .1])
    full = patch.expand(coefficients).reshape(len(y), len(x), 1)
    field = PlaneP1Field(x, y, full)
    expected = np.array([exterior_row(field, point, hx, hy, footprint, 16)
                         for point in patch.nodes[patch.active_nodes]])
    # Independent quadrature partitions differ: this is a numerical agreement
    # control, not an exact-integral certificate or a production tolerance.
    assert_allclose(operator.apply(coefficients), expected, atol=3e-4, rtol=3e-3)
    scale = 1.7
    enlarged = PhysicalP1Patch(scale*x, scale*y, scale*footprint)
    scaled_owner = SimpleNamespace(kernel=scale*owner.kernel, spacing=(scale*hx, scale*hy), shape=owner.shape)
    scaled = PhysicalWeakFlux(enlarged, scale*footprint, scaled_owner, order=8)
    assert_allclose(scaled.matrix, scale*operator.matrix, atol=2e-13, rtol=2e-12)
    assert operator.diagnostics['source_pairs'] > 0
    assert not operator.matrix.flags.writeable


def test_exterior_signed_toeplitz_order_and_budget_rejection():
    x, y = np.arange(5)*.3, np.arange(4)*.23
    footprint = box(4, 4, 5, 5)
    patch = PhysicalP1Patch(x, y, footprint)
    # Deliberately nonradial/asymmetric values expose reversal/index errors;
    # this test does not claim the synthetic kernel is a physical operator.
    kernel = np.arange(15, dtype=float).reshape(3, 5)
    owner = SimpleNamespace(kernel=kernel, spacing=(.3, .23), shape=(2, 3))
    operator = PhysicalWeakFlux(patch, footprint, owner)
    expected = np.array([[kernel[(i//3-j//3)+1, (i%3-j%3)+2] for j in range(6)] for i in range(6)])
    assert_allclose(operator.matrix, expected, atol=0, rtol=0)
    assert operator.diagnostics['source_pairs'] == 0
    with pytest.raises(MemoryError):
        PhysicalWeakFlux(patch, footprint, owner, maximum_bytes=1)
    cut = box(.41, .31, 2, 2)
    physical = PhysicalP1Patch(x, y, cut)
    with pytest.raises(RuntimeError, match='budget'):
        PhysicalWeakFlux(physical, cut, owner, maximum_pairs=1)
    with pytest.raises(ValueError, match='footprint'):
        PhysicalWeakFlux(physical, cut+np.array([.01, 0]), owner)
