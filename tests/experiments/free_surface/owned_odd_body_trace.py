"""True fluid-side odd-source body trace at explicit owned face interiors."""
import numpy as np
from topwater.rankine import TrianglePanels, single_layer_interior
from .panel_gradient_reference import checked_influence


def self_scalar_trace(panels, owners, barycentric):
    """Owned scalar self S in an intrinsic face chart, for ALL owned points."""
    owners, barycentric = np.asarray(owners), np.asarray(barycentric)
    if (owners.ndim != 1 or owners.dtype.kind not in 'iu' or np.any(owners < 0)
            or np.any(owners >= len(panels)) or barycentric.shape != (len(owners), 3)
            or np.iscomplexobj(barycentric) or not np.isfinite(barycentric).all()
            or np.any(barycentric <= 0)
            or np.any(abs(barycentric.sum(axis=1)-1) > 64*np.finfo(float).eps)):
        raise ValueError('positive original barycentric coordinates and valid owners required')
    values = np.empty(len(owners))
    metric_error = 0.
    relative_metric_error = 0.
    for owner in np.unique(owners):
        ids = np.flatnonzero(owners == owner)
        ab = panels.triangles[owner, 1]-panels.triangles[owner, 0]
        ac = panels.triangles[owner, 2]-panels.triangles[owner, 0]
        length = np.linalg.norm(ab)
        local = np.array([[0., 0., 0.], [length, 0., 0.],
            [np.dot(ac, ab)/length, np.linalg.norm(np.cross(ab, ac))/length, 0.]])
        chart = TrianglePanels(local[None], 1)
        error = float(abs(chart.edge_lengths[0]-panels.edge_lengths[owner]).max())
        if error > 128*np.finfo(float).eps*panels.diameters[owner]:
            raise ValueError('intrinsic source chart does not preserve original edge metric')
        metric_error = max(metric_error, error)
        relative_metric_error = max(relative_metric_error, error/panels.diameters[owner])
        for first in range(0, len(ids), 128):
            rows = ids[first:first+128]
            target = barycentric[rows]@local
            values[rows] = single_layer_interior(target, chart, np.zeros(len(rows), dtype=int))[:, 0]
    return values, dict(chart_edge_metric_max_absolute=metric_error,
                        chart_edge_metric_max_relative_diameter=relative_metric_error)


def owned_odd_maps(panels, points, owners, barycentric):
    """Return R/J maps and their true fluid-outward normal, without offsets.

    Only the flat self-panel normal PV is replaced by its analytic value zero,
    with exactly one +1/2 domain jump. Every other pair keeps bulk guards.
    Self scalar S uses the strict owned-interior API in its intrinsic chart.
    Original positive barycentrics are mandatory; no world point is projected.
    """
    points, owners = np.asarray(points), np.asarray(owners)
    if (not isinstance(panels, TrianglePanels) or owners.ndim != 1
            or points.shape != (len(owners), 3) or owners.dtype.kind not in 'iu'
            or np.any(owners < 0) or np.any(owners >= len(panels))
            or np.iscomplexobj(points) or not np.isfinite(points).all() or np.any(points[:, 2] >= 0)):
        raise ValueError('finite owned body-face interiors below the plane required')
    self_value, chart_diagnostics = self_scalar_trace(panels, owners, barycentric)
    barycentric = np.asarray(barycentric)
    expected = np.einsum('mj,mjk->mk', barycentric, panels.triangles[owners])
    coordinate_allowance = 64*np.finfo(float).eps*np.einsum('mj,mjk->mk',
        barycentric, abs(panels.triangles[owners]))
    if np.any(abs(points-expected) > coordinate_allowance):
        raise ValueError('world target does not correspond to its declared original barycentrics')
    potential = np.zeros((len(points), len(panels)))
    normal = np.zeros_like(potential)
    pairs = 0
    for owner in np.unique(owners):
        indices = np.flatnonzero(owners == owner)
        mask = np.arange(len(panels)) != owner
        other = TrianglePanels(panels.triangles[mask], panels.normal_sign[mask]) if mask.any() else None
        for first in range(0, len(indices), 128):
            ids = indices[first:first+128]
            target = points[ids]
            potential[ids, owner] = self_value[ids]
            normal[ids, owner] = .5
            if other is not None:
                influence = checked_influence(target, other)
                potential[np.ix_(ids, mask)] = influence.single_layer
                normal[np.ix_(ids, mask)] = np.einsum('mni,i->mn', influence.gradient_single_layer,
                                                     panels.normals[owner])
            image = checked_influence(target*[1., 1., -1.], panels)
            potential[ids] -= image.single_layer
            normal[ids] -= np.einsum('mni,i->mn', image.gradient_single_layer,
                                    panels.normals[owner]*[1., 1., -1.])
            pairs += len(ids)*2*len(panels)
    return potential, normal, dict(source_pairs=pairs, analytic_self_normal_jumps=len(points),
        **chart_diagnostics, maximum_world_barycentric_difference=float(abs(points-expected).max(initial=0)))


def owned_odd_trace(panels, sigma, points, owners, barycentric):
    """Apply the single owned map implementation to one finite source field."""
    sigma = np.asarray(sigma)
    if sigma.shape != (len(panels),) or not np.isfinite(sigma).all():
        raise ValueError('one finite source field required')
    potential, normal, diagnostics = owned_odd_maps(panels, points, owners, barycentric)
    return potential@sigma, normal@sigma, diagnostics
