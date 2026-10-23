"""Local analytic tangent of the existing guarded canonical-boat chart.

Columns are world translations, left/world rotations about the body CG, then
fixed-parent heights. No finite-step validity or roundoff bound is claimed.
Actual states still pass evaluate_template; separating aliases/incidences reject.
"""
import json

import numpy as np

from .boat_geometry_template import InvalidTemplate, _inputs, _plane, evaluate_template


def _plane_tangent(label, reference, world, roof, world_jacobian, roof_jacobian):
    kind, ids = label
    count = world_jacobian.shape[-1]
    if kind in ('bottom', 'background_edge'):
        return np.zeros((4, count))
    if kind == 'hull':
        selected = reference['source_body_faces'][ids[0]]
        points, tangent = world[selected], world_jacobian[selected]
    elif kind == 'roof':
        selected = reference['parents'][ids]
        points, tangent = roof[selected], roof_jacobian[selected]
    else:
        raise ValueError('unknown canonical boat plane definition')
    edges = points[1:]-points[:1]
    delta = tangent[1:]-tangent[:1]
    raw = np.cross(edges[0], edges[1])
    length = np.linalg.norm(raw)
    if not np.isfinite(length) or length == 0:
        raise InvalidTemplate('degenerate tangent plane')
    normal = raw/length
    draw = (np.cross(delta[0].T, edges[1])
            + np.cross(edges[0], delta[1].T)).T
    dn = (draw-normal[:, None]*(normal@draw)[None])/length
    return np.vstack((dn, points[0]@dn+normal@tangent[0]))


def _centroid_tangent(points, jacobian, loop, normal):
    """Differentiate the signed-area quotient used by _centroid."""
    if len(loop) == 3:
        # A triangle centroid is its vertex average. Differentiating the
        # equivalent area quotient introduces needless cancellation on slivers.
        return jacobian[loop].mean(axis=0)
    origin, do = points[loop[0]], jacobian[loop[0]]
    a, b = points[loop[1:-1]], points[loop[2:]]
    da, db = jacobian[loop[1:-1]], jacobian[loop[2:]]
    axis = np.argmax(abs(normal))
    weights = np.cross(a-origin, b-origin)[:, axis]
    dw = (np.cross((da-do).transpose(0, 2, 1), (b-origin)[:, None])
          + np.cross((a-origin)[:, None], (db-do).transpose(0, 2, 1)))[:, :, axis]
    measure = weights.sum()
    if not np.isfinite(measure) or measure == 0:
        raise InvalidTemplate('degenerate tangent polygon')
    numerator = np.einsum('f,fi->i', weights, origin+a+b)
    dnumerator = (np.einsum('fk,fi->ik', dw, origin+a+b)
                  + np.einsum('f,fik->ik', weights, do+da+db))
    center = numerator/(3*measure)
    return (dnumerator-3*center[:, None]*dw.sum(axis=0))/(3*measure)


def chart_tangent(reference, *, pose, parent_eta, parent_psi):
    """Return guarded current arrays and local geometry/trace Jacobians.

    point_jacobian has shape (N,3,6+Nparent); parent_map_jacobian has shape
    (Nfree,Nparent,6+Nparent). Alias and relevant-incidence tangent checks are
    engineering floating-point checks, not interval certificates.
    """
    current, base_diagnostics = evaluate_template(reference, pose=pose,
        parent_eta=parent_eta, parent_psi=parent_psi)
    _, world, roof = _inputs(reference, pose, parent_eta)
    count = 6+len(roof)
    world_jacobian = np.zeros((*world.shape, count))
    world_jacobian[:, :, :3] = np.eye(3)
    world_jacobian[:, :, 3:6] = np.cross(np.eye(3)[:, None],
        (world-pose.position)[None]).transpose(1, 2, 0)
    roof_jacobian = np.zeros((*roof.shape, count))
    roof_jacobian[:, 2, 6:] = np.eye(len(roof))
    planes, derivatives = [], []
    eps = np.finfo(float).eps
    aliases = json.loads(str(reference['template_plane_aliases']))
    for definitions in aliases:
        values = [_plane(label, reference, world, roof)*orientation
                  for label, orientation in definitions]
        tangents = [_plane_tangent(label, reference, world, roof,
            world_jacobian, roof_jacobian)*orientation for label, orientation in definitions]
        first = tangents[0]
        for other in tangents[1:]:
            budget = 512*eps*np.maximum(1., np.maximum(abs(first), abs(other)))
            if np.any(abs(other-first) > budget):
                raise InvalidTemplate('physical-plane alias derivatives separated')
        planes.append(values[0]); derivatives.append(first)
    planes, derivatives = np.asarray(planes), np.asarray(derivatives)
    triples = reference['template_owner_triples']
    matrices, right = planes[triples, :3], planes[triples, 3]
    owners = np.linalg.solve(matrices, right[..., None])[..., 0]
    rhs = derivatives[triples, 3]-np.einsum('mijk,mj->mik', derivatives[triples, :3], owners)
    owner_jacobian = np.linalg.solve(matrices, rhs)
    node_to_owner = reference['template_node_to_owner']
    selected = node_to_owner >= 0
    jacobian = np.zeros((*current['points'].shape, count))
    jacobian[selected] = owner_jacobian[node_to_owner[selected]]
    vertices, offsets = reference['template_fan_vertices'], reference['template_fan_offsets']
    for fan, (start, stop) in enumerate(zip(offsets[:-1], offsets[1:])):
        loop = vertices[start:stop]
        normal = planes[reference['template_fan_planes'][fan], :3]
        jacobian[reference['template_fan_centers'][fan]] = _centroid_tangent(
            current['points'], jacobian, loop, normal)
    witnesses = np.vstack((owners, current['points'][reference['template_fan_centers']]))
    witness_jacobian = np.concatenate((owner_jacobian,
        jacobian[reference['template_fan_centers']]))
    residual = (np.einsum('pi,wik->wpk', planes[:, :3], witness_jacobian)
                + np.einsum('pik,wi->wpk', derivatives[:, :3], witnesses)
                - derivatives[None, :, 3])
    incident = ((reference['template_reference_signs'] == 0)
                & reference['template_relevant_planes'])
    # Scale by the arithmetic terms, not a user-selected geometry tolerance.
    magnitude = (np.einsum('pi,wik->wpk', abs(planes[:, :3]), abs(witness_jacobian))
                 + np.einsum('pik,wi->wpk', abs(derivatives[:, :3]), abs(witnesses))
                 + abs(derivatives[None, :, 3]))
    budget = 2048*eps*np.maximum(1., magnitude)
    if np.any(abs(residual[incident]) > budget[incident]):
        raise InvalidTemplate('relevant incident-plane derivatives separated')
    if not np.isfinite(jacobian).all():
        raise InvalidTemplate('nonfinite chart tangent')
    mapping = np.zeros((len(current['free_nodes']), len(roof), count))
    for row, node in enumerate(current['free_nodes']):
        support = current['free_parent_indices'][row]
        support = support[support >= 0]
        xy = roof[support, :2]
        moving = jacobian[node, :2]
        if len(support) == 1:
            local = np.zeros((1, count))
        elif len(support) == 2:
            edge = xy[1]-xy[0]
            dt = edge@moving/(edge@edge)
            local = np.stack((-dt, dt))
        else:
            duv = np.linalg.solve((xy[1:]-xy[0]).T, moving)
            local = np.vstack((-duv.sum(axis=0), duv))
        mapping[row, support] = local
    return current, dict(point_jacobian=jacobian, parent_map_jacobian=mapping,
        maximum_incident_derivative_residual=float(np.max(abs(residual[incident]), initial=0.)),
        base_diagnostics=base_diagnostics,
        direction_order=('world_translation[3]', 'left_world_rotation[3]', 'parent_eta'))
