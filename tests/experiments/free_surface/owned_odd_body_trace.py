"""True fluid-side odd-source body trace at explicit owned face interiors."""
import numpy as np
from topwater.rankine import TrianglePanels, single_layer_interior
from .panel_gradient_reference import checked_influence


def owned_odd_trace(panels, sigma, points, owners):
    """Return R sigma and its fluid-outward normal, without normal offsets.

    Only the flat self-panel normal PV is replaced by its analytic value zero,
    with exactly one +1/2 domain jump. Every other pair keeps bulk guards.
    Scalar S uses the existing strict owned-interior API, never edge snapping.
    """
    points, owners, sigma = np.asarray(points), np.asarray(owners), np.asarray(sigma)
    if (not isinstance(panels, TrianglePanels) or owners.ndim != 1
            or points.shape != (len(owners), 3) or owners.dtype.kind not in 'iu'
            or np.any(owners < 0) or np.any(owners >= len(panels))
            or sigma.shape != (len(panels),) or not np.isfinite(sigma).all()
            or np.iscomplexobj(points) or not np.isfinite(points).all() or np.any(points[:, 2] >= 0)):
        raise ValueError('finite owned body-face interiors below the plane and one source field required')
    potential = np.empty(len(points), dtype=np.result_type(sigma, float))
    normal = np.empty_like(potential)
    pairs = 0
    for owner in np.unique(owners):
        indices = np.flatnonzero(owners == owner)
        mask = np.arange(len(panels)) != owner
        other = TrianglePanels(panels.triangles[mask], panels.normal_sign[mask]) if mask.any() else None
        for first in range(0, len(indices), 128):
            ids = indices[first:first+128]
            target = points[ids]
            scalar = single_layer_interior(target, panels, np.full(len(ids), owner))
            value = scalar@sigma
            derivative = np.full(len(ids), .5*sigma[owner], dtype=normal.dtype)
            if other is not None:
                influence = checked_influence(target, other)
                derivative += np.einsum('mni,i->mn', influence.gradient_single_layer,
                                        panels.normals[owner])@sigma[mask]
            image = checked_influence(target*[1., 1., -1.], panels)
            value -= image.single_layer@sigma
            derivative -= np.einsum('mni,i->mn', image.gradient_single_layer,
                                    panels.normals[owner]*[1., 1., -1.])@sigma
            potential[ids], normal[ids] = value, derivative
            pairs += len(ids)*(3*len(panels)-1)
    return potential, normal, dict(source_pairs=pairs, analytic_self_normal_jumps=len(points))
