"""Composite target integration; physical panels and unknowns stay unchanged."""
import numpy as np
from numpy.polynomial.legendre import leggauss

from tests.laplace_reference import triangle_samples


def composite_samples(triangles, order, levels):
    """Red subdivide integration regions, returning original-area weights."""
    if levels not in (1, 2):
        raise ValueError('bounded composite probe requires one or two levels')
    triangles = np.asarray(triangles)
    children = triangles[:, None]
    for _ in range(levels):
        a, b, c = np.moveaxis(children, 2, 0)
        ab, bc, ca = (a+b)/2, (b+c)/2, (c+a)/2
        children = np.stack((np.stack((a, ab, ca), axis=2),
                             np.stack((ab, b, bc), axis=2),
                             np.stack((ca, bc, c), axis=2),
                             np.stack((ab, bc, ca), axis=2)), axis=2)
        children = children.reshape(len(triangles), -1, 3, 3)
    count = children.shape[1]
    points, weights = triangle_samples(children.reshape(-1, 3, 3), order)
    return points.reshape(len(triangles), -1, 3), np.tile(weights/count, count)


def edge_samples(triangles, order, power):
    """Average three edge-clustered rules; same rule for every source in a row.

    s=t**power, x=s*A+(1-s)*((1-v)*B+v*C). The normalized area
    Jacobian is 2*(1-s)*power*t**(power-1). Order six integrates all
    affine target moments exactly for powers two and four.
    """
    if power not in (2, 4):
        raise ValueError('bounded edge probe requires power two or four')
    nodes, gauss = leggauss(order)
    unit, gauss = (nodes+1)/2, gauss/2
    t, v = np.meshgrid(unit, unit, indexing='ij')
    s = t**power
    barycentric = np.column_stack((s.ravel(), ((1-s)*(1-v)).ravel(), ((1-s)*v).ravel()))
    weights = (2*(1-s)*power*t**(power-1)*np.outer(gauss, gauss)).ravel()
    barycentric = np.concatenate([np.roll(barycentric, shift, axis=1) for shift in range(3)])
    return barycentric@triangles, np.tile(weights/3, 3)


def prepare_composite_energy(faces, panels, order, levels):
    points, weights = composite_samples(panels.triangles, order, levels)
    return _prepare_sampled_energy(faces, panels, order, points, weights)


def prepare_edge_energy(faces, panels, order, power):
    points, weights = edge_samples(panels.triangles, order, power)
    return _prepare_sampled_energy(faces, panels, order, points, weights)


def _prepare_sampled_energy(faces, panels, order, points, weights):
    """Same pair kernel and single self jump, with common points per target row.

    This bounded, single-process experiment restores the kernel even on failure.
    Its returned assembly is explicitly marked nonuniform by prepare_energy.
    """
    from . import cuda_operator as gpu
    device_points = gpu.cp.asarray(np.ascontiguousarray(points))
    device_weights = gpu.cp.asarray(weights)
    original = gpu._pairs

    def dispatch(grid, block, args):
        args = list(args)
        args[10:12] = device_points, device_weights
        args[14] = np.int32(len(weights))
        original(grid, block, tuple(args))

    gpu._pairs = dispatch
    try:
        return gpu.prepare_energy(faces, panels, order=order)
    finally:
        gpu._pairs = original
