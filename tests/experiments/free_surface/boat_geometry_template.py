"""Guarded CURRENT geometry on one canonical-boat arrangement chart.

No operator is cached. Uncertain arithmetic or changed topology raises
InvalidTemplate; the caller must construct a fresh exact reference instead.
Roundoff radii are engineering guards, not a formal interval proof; acceptance
must be checked against exact-reference geometry and current physical fields.
"""
from itertools import combinations
from hashlib import sha256
import json
import numpy as np

from topwater import HullMesh
from topwater.frames import rotation_matrix


class InvalidTemplate(ValueError):
    pass


class InvalidFan(InvalidTemplate):
    """Physical chart guards passed, but its chosen triangulation is invalid."""


def chart_key(reference):
    """Identity of immutable node recipes/ownership, independent of pose/eta/psi."""
    names = sorted(name for name in reference if name.startswith('template_'))
    names += ['source_body_vertices', 'source_body_faces', 'source_body_cg', 'parents',
              'body_faces', 'free_faces', 'wall_faces', 'bottom_faces', 'free_nodes',
              'body_source_group', 'free_source_parent', 'free_parent_indices']
    digest = sha256()
    for name in [*names, 'parent_xy']:
        value = (reference['parent_nodes'][:, :2] if name == 'parent_xy' else reference[name])
        value = np.ascontiguousarray(value)
        digest.update(name.encode()); digest.update(value.dtype.str.encode())
        digest.update(str(value.shape).encode()); digest.update(value.tobytes())
    return np.array(digest.hexdigest())


def _plane(label, arrays, world, roof):
    kind, ids = label
    if kind == 'hull':
        p = world[arrays['source_body_faces'][ids[0]]]
    elif kind == 'roof':
        p = roof[arrays['parents'][ids]]
    elif kind == 'bottom':
        p = roof[arrays['parents'][ids]].copy()[::-1]
        p[:, 2] = -2.
    elif kind == 'background_edge':
        a, b = roof[list(ids)]
        normal = np.array([b[1]-a[1], a[0]-b[0], 0.])
        normal /= np.linalg.norm(normal)
        return np.r_[normal, normal@a]
    else:
        raise ValueError('unknown canonical boat plane definition')
    normal = np.cross(p[1]-p[0], p[2]-p[0])
    normal /= np.linalg.norm(normal)
    return np.r_[normal, normal@p[0]]


def _inputs(arrays, pose, eta):
    eta = np.asarray(eta, dtype=float)
    if eta.shape != (len(arrays['parent_nodes']),) or not np.isfinite(eta).all():
        raise ValueError('one finite eta per fixed parent node required')
    rotation = rotation_matrix(pose.orientation)
    world = (arrays['source_body_vertices']-arrays['source_body_cg'])@rotation.T+pose.position
    roof = np.column_stack((arrays['parent_nodes'][:, :2], eta))
    return rotation, world, roof


def export_template(registry, meta, arrays, critical_vertices):
    """Cold exact witness export, including DRY hull/prism vertices."""
    from topwater import InitialState
    pose = InitialState(arrays['pose_position'], arrays['pose_orientation'])
    _, world, roof = _inputs(arrays, pose, arrays['parent_nodes'][:, 2])
    canonical = np.array(registry.planes, dtype=float)
    canonical /= np.linalg.norm(canonical[:, :3], axis=1)[:, None]
    aliases = []
    for plane, labels in enumerate(registry.labels):
        definitions = []
        seen = set()
        for label in labels:
            key = json.dumps(label, default=int)
            if key not in seen:
                raw = _plane(label, arrays, world, roof)
                definitions.append((label, 1 if raw[:3]@canonical[plane, :3] > 0 else -1))
                seen.add(key)
        aliases.append(definitions)
    active = meta['active']
    owners = sorted(set(critical_vertices)|{int(i) for i in active if i < len(registry.points)})
    index = {owner: row for row, owner in enumerate(owners)}
    node_to_owner = np.array([index[int(i)] if i < len(registry.points) else -1 for i in active])
    inverse = {int(i): node for node, i in enumerate(active)}
    fan_vertices, offsets, centers, fan_planes, fan_signs = [], [0], [], [], []
    for polygon, plane, sign, _, _ in meta['polygons']:
        fan_vertices.extend(inverse[i] for i in polygon)
        offsets.append(len(fan_vertices)); fan_planes.append(plane); fan_signs.append(sign)
    # Each polygon's fan is contiguous in the reference face list.
    face_start = np.r_[0, np.cumsum(np.diff(offsets))[:-1]]
    centers = arrays['faces'][face_start, 0]
    exact_points = [registry.points[i] for i in owners]+[meta['points'][active[i]] for i in centers]
    signs = registry.exact_signs(exact_points)
    if registry.relevant is None:
        raise ValueError('template export requires actual clipping-constraint provenance')
    relevant = np.zeros(signs.shape, dtype=bool)
    for row, owner in enumerate(owners):
        relevant[row, list(registry.relevant[owner])] = True
    for fan, plane in enumerate(fan_planes):
        # Boundary vertices preserve each real clipping piece. A merged body
        # fan need not stay on one artificial column; only its own facet plane
        # and strictly oriented interior fan are constraints of the center.
        relevant[len(owners)+fan, plane] = True
    triples = []
    for point in owners:
        choices = list(combinations(sorted(registry.incidence[point]&registry.relevant[point]), 3))
        determinants = np.abs(np.linalg.det(canonical[np.asarray(choices), :3]))
        if not len(choices) or max(determinants) == 0:
            raise ValueError('reference witness has no independent owner triple')
        triples.append(choices[int(np.argmax(determinants))])
    result = dict(template_plane_aliases=np.array(json.dumps(aliases, default=int)),
        template_faces=arrays['faces'].copy(),
        template_owner_triples=np.asarray(triples, dtype=int),
        template_owner_determinant_sign=np.sign(np.linalg.det(canonical[np.asarray(triples), :3])),
        template_node_to_owner=node_to_owner,
        template_reference_signs=signs,
        template_relevant_planes=relevant,
        template_fan_vertices=np.asarray(fan_vertices), template_fan_offsets=np.asarray(offsets),
        template_fan_centers=centers, template_fan_planes=np.asarray(fan_planes),
        template_fan_orientation=np.asarray(fan_signs))
    result['node_chart_key'] = chart_key(dict(arrays, **result))
    return result


def _centroid(points, loop, normal):
    origin = points[loop[0]]
    a, b = points[loop[1:-1]], points[loop[2:]]
    axis = np.argmax(abs(normal))
    weights = np.cross(a-origin, b-origin)[:, axis]
    measure = weights.sum()
    if not np.isfinite(measure) or measure == 0:
        raise InvalidTemplate('degenerate current polygon')
    return (weights[:, None]*(origin+a+b)).sum(axis=0)/(3*measure)


def evaluate_template(reference, *, pose, parent_eta, parent_psi):
    """Reconstruct this chart or reject it. No cold refresh happens implicitly."""
    if (len(reference['template_node_to_owner']) != len(reference['points'])
            or not np.array_equal(reference['faces'], reference['template_faces'])):
        raise InvalidTemplate('reference template was reindexed or refined')
    key = chart_key(reference)
    if 'node_chart_key' in reference and reference['node_chart_key'] != key:
        raise InvalidTemplate('reference chart recipes or ownership changed')
    rotation, world, roof = _inputs(reference, pose, parent_eta)
    psi = np.asarray(parent_psi, dtype=float)
    if psi.shape != (len(roof),) or not np.isfinite(psi).all():
        raise ValueError('one finite psi per parent node required')
    scale = max(1., np.max(abs(world)), np.max(abs(roof)))
    eps = np.finfo(float).eps
    plane_error = 64*eps*scale
    aliases = json.loads(str(reference['template_plane_aliases']))
    planes = []
    for definitions in aliases:
        current = [_plane(label, reference, world, roof)*orientation for label, orientation in definitions]
        first = current[0]
        if any(np.max(abs(other[:3]-first[:3])) > 64*eps or abs(other[3]-first[3]) > plane_error
               for other in current[1:]):
            raise InvalidTemplate('reference physical-plane aliases separated')
        planes.append(first)
    planes = np.asarray(planes)
    triples = reference['template_owner_triples']
    matrices, right = planes[triples, :3], planes[triples, 3]
    if np.any(np.sign(np.linalg.det(matrices)) != reference['template_owner_determinant_sign']):
        raise InvalidTemplate('owner-plane orientation changed')
    try:
        owner_points = np.linalg.solve(matrices, right[..., None])[..., 0]
        inverse = np.linalg.inv(matrices)
    except np.linalg.LinAlgError as error:
        raise InvalidTemplate('current owner planes are dependent') from error
    condition = np.linalg.cond(matrices)
    residual = np.max(abs(np.einsum('mij,mj->mi', matrices, owner_points)-right), axis=1)
    radius = np.linalg.norm(inverse, ord=np.inf, axis=(1, 2))*(residual+plane_error)
    if not np.isfinite(owner_points).all() or not np.isfinite(radius).all():
        raise InvalidTemplate('nonfinite owner intersection')
    points = np.full(reference['points'].shape, np.nan)
    node_to_owner = reference['template_node_to_owner']
    selected = node_to_owner >= 0
    points[selected] = owner_points[node_to_owner[selected]]
    node_radius = np.zeros(len(points))
    node_radius[selected] = radius[node_to_owner[selected]]
    fan_radii, invalid_fans = [], []
    vertices, offsets = reference['template_fan_vertices'], reference['template_fan_offsets']
    for fan, (start, stop) in enumerate(zip(offsets[:-1], offsets[1:])):
        loop = vertices[start:stop]
        normal = planes[reference['template_fan_planes'][fan], :3]
        center = _centroid(points, loop, normal)
        node = reference['template_fan_centers'][fan]
        points[node] = center
        n = np.cross(points[loop]-center, points[np.roll(loop, -1)]-center)
        sectors = n@normal
        if not np.all(sectors*reference['template_fan_orientation'][fan] > 0):
            invalid_fans.append(fan)
        fan_radius = 8*max(node_radius[loop])+plane_error
        node_radius[node] = fan_radius
        fan_radii.append(fan_radius)
    witnesses = np.vstack((owner_points, points[reference['template_fan_centers']]))
    errors = np.r_[radius, fan_radii]
    values = witnesses@planes[:, :3].T-planes[:, 3]
    bounds = errors[:, None]*np.abs(planes[:, :3]).sum(axis=1)+plane_error
    signs = reference['template_reference_signs']
    relevant = reference['template_relevant_planes']
    incident = (signs == 0)&relevant
    if np.any(abs(values[incident]) > bounds[incident]):
        raise InvalidTemplate('reference plane incidence split')
    nonincident = (signs != 0)&relevant
    if np.any((signs*values)[nonincident] <= bounds[nonincident]):
        raise InvalidTemplate('changed or uncertain nonincident plane sign')
    # Conservative sufficient guard: otherwise let the exact reference decide.
    rim = world[reference['source_body_rim']]
    if np.min(rim[:, 2])-np.max(roof[:, 2]) <= plane_error:
        raise InvalidTemplate('cannot certify rim clearance from global height bound')
    if (np.min(world[:, 2]) <= -2+plane_error
            or np.any(world[:, :2].min(axis=0) <= roof[:, :2].min(axis=0)+plane_error)
            or np.any(world[:, :2].max(axis=0) >= roof[:, :2].max(axis=0)-plane_error)):
        raise InvalidTemplate('body reaches the bed or fixed outer walls')
    if invalid_fans:
        # Distinguish a failed numerical fan from a failed physical clipping
        # arrangement only AFTER every owner/birth/sign/domain guard passes.
        raise InvalidFan('current polygon lost its interior oriented fan: '+str(invalid_fans))
    try:
        HullMesh(points, reference['faces'])
    except ValueError as error:
        raise InvalidTemplate('current closed mesh failed its geometric guard') from error
    result = {key: value for key, value in reference.items() if not key.startswith('template_')}
    result['node_chart_key'] = key
    result.update(points=points, parent_nodes=roof, parent_psi=psi.copy(), body_vertices=world,
                  pose_position=pose.position.copy(), pose_orientation=pose.orientation.copy())
    material = np.linalg.solve(rotation, (points-pose.position).T).T+reference['source_body_cg']
    result['body_material_vertices'] = material[result['faces'][result['body_faces']]]
    result.pop('body_vertex_source_triangle', None)
    result.pop('body_vertex_source_barycentric', None)
    # Existing support IDs are semantic chart ownership, but weights are CURRENT.
    indices = reference['free_parent_indices']
    weights = np.zeros(indices.shape)
    for row, node in enumerate(reference['free_nodes']):
        support = indices[row][indices[row] >= 0]
        xy = roof[support, :2]
        target = points[node, :2]
        if len(support) == 1:
            local = np.ones(1)
        elif len(support) == 2:
            edge = xy[1]-xy[0]
            t = (target-xy[0])@edge/(edge@edge)
            local = np.array([1-t, t])
        else:
            uv = np.linalg.solve((xy[1:]-xy[0]).T, target-xy[0])
            local = np.r_[1-uv.sum(), uv]
        coordinate_budget = 8*node_radius[node]+plane_error
        if (np.max(abs(local@xy-target)) > coordinate_budget
                or abs(local@roof[support, 2]-points[node, 2]) > coordinate_budget
                or np.min(local) < -64*eps):
            raise InvalidTemplate('current free trace left its fixed-parent support')
        weights[row, :len(support)] = local
    result['free_parent_weights'] = weights
    triangles = points[result['faces']]
    fluid_volume = np.einsum('fi,fi->f', triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2])).sum()/6
    p = roof[result['parents']]
    area = np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0])[:, 2]/2
    tank_volume = area@(p[:, :, 2].mean(axis=1)+2)
    if not len(result['body_faces']) or tank_volume-fluid_volume <= 64*eps*max(1., tank_volume):
        raise InvalidTemplate('cannot certify nonzero displaced body volume')
    return result, dict(maximum_owner_condition=float(max(condition)),
        maximum_intersection_uncertainty_m=float(max(radius)),
        rim_clearance_lower_bound_m=float(np.min(rim[:, 2])-np.max(roof[:, 2])),
        displaced_volume_m3=float(tank_volume-fluid_volume), witnesses=len(witnesses))
