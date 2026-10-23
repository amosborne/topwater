"""Shared odd-source work maps for two actual fixed-hull verification callers."""
from time import perf_counter
import numpy as np
from tests.laplace_reference import triangle_samples
from .compact_p1_maps import body_maps
from .panel_gradient_reference import checked_influence


def source_flux(patch, panels, order, *, maximum_pairs, maximum_bytes):
    """Direct exterior F=2 int Ni*partial_z(S), without hidden subtraction."""
    started = perf_counter()
    count = len(panels)
    pairs = len(patch.cut_triangles)*order*order*count
    if pairs > maximum_pairs or 24*len(patch.nodes)*count > maximum_bytes:
        raise RuntimeError('odd source-flux assembly preflight exceeds work/storage budget')
    result = np.zeros((len(patch.nodes), count))
    bary, weights = triangle_samples(np.eye(3)[None], order)
    bary = bary[0]
    work = 0
    for first in range(0, len(patch.cut_triangles), 16):
        last = min(first+16, len(patch.cut_triangles))
        xy = np.einsum('qj,tji->tqi', bary, patch.cut_triangles[first:last])
        basis = np.einsum('qj,tji->tqi', bary, patch.basis_values[first:last])
        measure = (patch.cut_areas[first:last, None, None]*weights[None, :, None]*basis).reshape(-1, 3)
        points = np.column_stack((xy.reshape(-1, 2), np.zeros(len(xy)*len(weights))))
        nodes = np.repeat(patch.faces[patch.parent_id[first:last]], len(weights), axis=0)
        for start in range(0, len(points), 128):
            stop = min(start+128, len(points))
            work += (stop-start)*count
            if work > maximum_pairs:
                raise RuntimeError('odd F source-pair budget exhausted')
            # Gradient source-edge guards are preserved; scalar plane S is
            # not a fallback for a genuinely singular gradient query.
            flux = 2*checked_influence(points[start:stop], panels).gradient_single_layer[:, :, 2]
            for local in range(3):
                np.add.at(result, nodes[start:stop, local], measure[start:stop, local, None]*flux)
    return result[patch.active_nodes], dict(source_pairs=work, assembly_seconds=perf_counter()-started)


def harmonic_field(patch, psi, points, normals):
    """H psi at fixed old probe points; caller may reuse across body partitions."""
    full = patch.expand(psi)
    p, n = body_maps(patch, points, normals)
    return p@full, n@full


def source_field(panels, sigma, points, normals):
    """Odd body-source potential and fluid-side normal at strict bulk targets."""
    potential = np.empty(len(points), dtype=np.result_type(sigma, float))
    normal = np.empty_like(potential)
    for first in range(0, len(points), 128):
        last = min(first+128, len(points))
        real = checked_influence(points[first:last], panels)
        image = checked_influence(points[first:last]*[1., 1., -1.], panels)
        potential[first:last] = (real.single_layer-image.single_layer)@sigma
        gradient = real.gradient_single_layer-image.gradient_single_layer*[1., 1., -1.]
        normal[first:last] = np.einsum('mni,mi->mn', gradient, normals[first:last])@sigma
    return potential, normal
