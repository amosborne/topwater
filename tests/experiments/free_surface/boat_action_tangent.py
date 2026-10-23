"""Exact small action-object tangents on one currently owned boat chart.

Columns are world CG translations, left/world rotations, then parent heights.
No operator is differentiated here; the current vertex energy gradient is an
explicit input. The finite-difference action path remains the reference.
"""
import numpy as np

from tests.reference_boat import GRAVITY, WATER_DENSITY
from .boat_parent_action import parent_mass, factor_parent_mass, port_transport_covector
from .boat_world_field import refine_boundary


def refine_tangent(arrays, point_jacobian, map_jacobian, levels, *, local_body_refinement=None):
    """Prolong both tangents through the existing sorted body-edge midpoints.

    The midpoint coordinates are checked bitwise against the refiner owner;
    neither the red/green face split nor its validity rules are reimplemented.
    """
    current, dx = arrays, np.asarray(point_jacobian)
    ds = np.zeros((len(current['points']), len(current['parent_nodes']), dx.shape[-1]))
    ds[current['free_nodes']] = map_jacobian
    def refine(current, dx, ds, selected=None):
        faces = current['faces'][current['body_faces'] if selected is None else selected]
        edges = np.unique(np.sort(np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]])), axis=1), axis=0)
        refined = refine_boundary(current, region='body', marked_faces=selected)
        expected = np.vstack((current['points'], current['points'][edges].mean(axis=1)))
        if not np.array_equal(refined['points'], expected):
            raise ValueError('body refinement no longer follows the owned midpoint prolongation')
        dx = np.concatenate((dx, dx[edges].mean(axis=1)))
        ds = np.concatenate((ds, ds[edges].mean(axis=1)))
        return refined, dx, ds
    for _ in range(levels):
        current, dx, ds = refine(current, dx, ds)
    local = local_body_refinement or {}
    for _ in range(local.get('levels', 0)):
        selected = current['body_faces'][np.any(current['faces'][current['body_faces']] == local['contact_node'], axis=1)]
        current, dx, ds = refine(current, dx, ds, selected)
    return current, dx, ds[current['free_nodes']]


def object_derivatives(arrays, point_jacobian, map_jacobian):
    """Differentiate B, C and fluid gravity polynomials at current nodes."""
    points, faces = arrays['points'], arrays['faces']
    nparent = len(arrays['parent_nodes'])
    ndirection = 6+nparent
    dx = np.asarray(point_jacobian)
    ds = np.asarray(map_jacobian)
    if (dx.shape != (*points.shape, ndirection) or ds.shape != (len(arrays['free_nodes']), nparent, ndirection)
            or not np.isfinite(dx).all() or not np.isfinite(ds).all()):
        raise ValueError('finite current point/map Jacobians with all owned directions required')
    tri = points[faces]
    dtri = np.moveaxis(dx[faces], -1, 0)
    e1, e2 = tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]
    da = .5*(np.cross(dtri[:, :, 1]-dtri[:, :, 0], e2)
              +np.cross(e1, dtri[:, :, 2]-dtri[:, :, 0]))
    area_vector = .5*np.cross(e1, e2)
    # Exact P1 rigid-motion load uses oriented area vectors: no normal/area
    # division is needed in this derivative.
    body = arrays['body_faces']
    dcg = np.zeros((ndirection, 3)); dcg[:3] = np.eye(3)
    lever = tri[body]-arrays['pose_position']
    dlever = dtri[:, body]-dcg[:, None, None, :]
    drot = (np.cross(dlever, area_vector[None, body, None, :])
            +np.cross(lever[None], da[:, body, None, :]))
    dtranslation = np.broadcast_to(da[:, body, None, :], drot.shape)
    dmodes = np.concatenate((dtranslation, drot), axis=-1)
    db = np.zeros((len(points), 6, ndirection))
    for corner in range(3):
        values = (dmodes.sum(axis=2)+dmodes[:, :, corner])/12
        np.add.at(db, faces[body, corner], np.moveaxis(values, 0, -1))
    mapping = np.zeros((len(points), nparent))
    for column in range(3):
        ids = arrays['free_parent_indices'][:, column]
        valid = ids >= 0
        mapping[arrays['free_nodes'][valid], ids[valid]] += arrays['free_parent_weights'][valid, column]
    full_ds = np.zeros((ndirection, len(points), nparent))
    full_ds[:, arrays['free_nodes']] = np.moveaxis(ds, -1, 0)
    free = arrays['free_faces']
    trace, dtrace = mapping[faces[free]], full_ds[:, faces[free]]
    unit = (np.ones((3, 3))+np.eye(3))/12
    area, darea = area_vector[free, 2], da[:, free, 2]
    if np.any(area <= 0):
        raise ValueError('invalid projected free triangle')
    dmass = np.einsum('df,fci,ce,fej->dij', darea, trace, unit, trace, optimize=True)
    left = np.einsum('f,dfci,ce,fej->dij', area, dtrace, unit, trace, optimize=True)
    dmass += left+left.swapaxes(1, 2)
    z, dz = tri[:, :, 2], dtri[:, :, :, 2]
    average_squared = (z.sum(axis=1)**2+(z*z).sum(axis=1))/12
    daverage = (z.sum(axis=1)*dz.sum(axis=2)+(z[None]*dz).sum(axis=2))/6
    dgravity = WATER_DENSITY*GRAVITY/2*np.sum(da[:, :, 2]*average_squared+area_vector[:, 2]*daverage, axis=1)
    return dict(body_load_jacobian=db, parent_mass_jacobian=dmass, gravity_jacobian=dgravity)


def action_terms(arrays, field, gradient, values, twist, time, point_jacobian, map_jacobian):
    """All parent rates and body pose loads, without shape FD evaluations."""
    if gradient.shape != arrays['points'].shape or not np.isfinite(gradient).all():
        raise ValueError('finite current vertex energy gradient required')
    derivatives = object_derivatives(arrays, point_jacobian, map_jacobian)
    mass = parent_mass(arrays)
    solve, condition = factor_parent_mass(mass)
    eta_rate = solve(field['parent_map'].T@field['flux'])
    eta, psi = arrays['parent_nodes'][:, 2], arrays['parent_psi']
    covector = port_transport_covector(arrays, field['potential'], values, time)
    kinetic = -np.einsum('ni,nid->d', gradient, point_jacobian)
    body_load = WATER_DENSITY*np.einsum('n,nkd,k->d', field['potential'], derivatives['body_load_jacobian'], twist)
    moving_trace = -WATER_DENSITY*np.einsum('n,npd,p->d', field['flux'], map_jacobian, psi)
    incident_port = -WATER_DENSITY*np.einsum('ni,nid->d', covector, point_jacobian)
    dc = derivatives['parent_mass_jacobian']
    cdot = np.einsum('d,dij->ij', np.r_[twist, eta_rate], dc)
    one_form = WATER_DENSITY*np.einsum('i,dij,j->d', psi, dc, eta_rate)
    parent_terms = dict(kinetic=kinetic[6:], body_load=body_load[6:], moving_trace=moving_trace[6:],
        incident_port=incident_port[6:], gravity=-WATER_DENSITY*GRAVITY*(mass@eta),
        mass_transport=one_form[6:]-WATER_DENSITY*(cdot.T@psi))
    weak = sum(parent_terms.values())
    body_terms = dict(kinetic=kinetic[:6], boundary_load=body_load[:6], moving_trace=moving_trace[:6],
        incident_port=incident_port[:6], gravity=-derivatives['gravity_jacobian'][:6], surface_one_form=one_form[:6])
    body_pose = sum(body_terms.values())
    return dict(eta_rate=eta_rate, psi_rate=solve(weak/WATER_DENSITY), parent_weak_load=weak,
        body_pose_load=body_pose, parent_terms=parent_terms, body_terms=body_terms,
        parent_mass_condition=condition, parent_mass=mass, mass_time_derivative=cdot, object_derivatives=derivatives)
