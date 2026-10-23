"""Restricted parent loads with prescribed open-port virtual work.

Diagnostic candidate, not enabled in the moving solver. The supplied vertex
gradient must belong to this current geometry, potential, density and energy.
Only geometry/load/trace stencils are evaluated here: no fluid operator/solve.
"""
from collections import Counter

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cho_factor, cho_solve
from scipy.spatial.transform import Rotation

from topwater import InitialState
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.reference_boat import GRAVITY, WATER_DENSITY
from .boat_geometry_template import InvalidTemplate
from .boat_world_field import boundary_data
from .incident_baseline_benchmark import fused_cpu


def parent_mass(arrays):
    """Integrate fixed world parent basis products over exposed projected faces."""
    mapping = np.zeros((len(arrays['points']), len(arrays['parent_nodes'])))
    for column in range(3):
        ids = arrays['free_parent_indices'][:, column]
        valid = ids >= 0
        mapping[arrays['free_nodes'][valid], ids[valid]] += arrays['free_parent_weights'][valid, column]
    faces = arrays['faces'][arrays['free_faces']]
    triangles = arrays['points'][faces, :2]
    area = np.linalg.det(triangles[:, 1:]-triangles[:, :1])/2
    if np.any(area <= 0):
        raise ValueError('invalid projected free triangle')
    trace = mapping[faces]
    unit = (np.ones((3, 3))+np.eye(3))/12
    return np.einsum('f,fci,cd,fdj->ij', area, trace, unit, trace, optimize=True)


def port_transport_covector(arrays, potential, values, time):
    """Current integral of q*grad_Gamma(phi)*N_i on wall AND bottom ports."""
    points = arrays['points']
    faces = arrays['faces'][np.r_[arrays['wall_faces'], arrays['bottom_faces']]]
    panels = TrianglePanels(points[faces], 1)
    samples, weights = triangle_samples(panels.triangles, 3)
    shape = triangle_samples(np.eye(3)[None], 3)[0][0]
    velocity = fused_cpu(samples.reshape(-1, 3), values, time)[:, 6:].reshape(samples.shape)
    flux = np.einsum('fqi,fi->fq', velocity, panels.normals)
    opposite = np.roll(panels.triangles, -2, axis=1)-np.roll(panels.triangles, -1, axis=1)
    basis_gradient = np.cross(panels.geometric_normals[:, None], opposite)/(2*panels.areas[:, None, None])
    gradient = np.einsum('fci,fc->fi', basis_gradient, potential[faces])
    covector = np.zeros_like(points)
    for corner in range(3):
        load = panels.areas*np.einsum('q,q,fq->f', weights, shape[:, corner], flux)
        np.add.at(covector, faces[:, corner], load[:, None]*gradient)
    return covector


def port_seam_load(arrays, field, values, time, order=12):
    """Independent line integral -rho*int_perimeter psi*q*R_j ds_horizontal."""
    points, faces = arrays['points'], arrays['faces']
    free_edges = Counter(tuple(sorted((int(a), int(b))))
        for face in faces[arrays['free_faces']] for a, b in zip(face, np.roll(face, -1)))
    wall_ids = arrays['wall_faces']
    panels = TrianglePanels(points[faces[wall_ids]], 1)
    edges, normals = [], []
    for face, normal in zip(faces[wall_ids], panels.normals, strict=True):
        for a, b in zip(face, np.roll(face, -1)):
            edge = tuple(sorted((int(a), int(b))))
            if free_edges[edge] == 1:
                edges.append(edge); normals.append(normal)
    if not edges or len(set(edges)) != len(edges):
        raise ValueError('missing or multiply owned free/wall seam')
    edges, normals = np.array(edges), np.array(normals)
    nodes, weights = leggauss(order)
    line_shape = np.column_stack(((1-nodes)/2, (1+nodes)/2))
    weights = weights/2
    samples = np.einsum('qc,eci->eqi', line_shape, points[edges])
    velocity = fused_cpu(samples.reshape(-1, 3), values, time)[:, 6:].reshape(samples.shape)
    flux = np.einsum('eqi,ei->eq', velocity, normals)
    mapping = np.zeros((len(points), len(arrays['parent_nodes'])))
    mapping[arrays['free_nodes']] = field['parent_map']
    trace = np.einsum('qc,eci->eqi', line_shape, mapping[edges])
    psi = np.einsum('qc,ec->eq', line_shape, field['potential'][edges])
    length = np.linalg.norm(points[edges[:, 1], :2]-points[edges[:, 0], :2], axis=1)
    return -WATER_DENSITY*np.einsum('e,q,eq,eq,eqi->i', length, weights, psi, flux, trace)


def action_rates(arrays, field, gradient, values, twist, time, geometry):
    """Candidate all-parent rates from one CURRENT fixed-nodal energy gradient.

    geometry(pose, eta) must replay the same owned/refined chart at fixed psi.
    Both h/h2 results are returned; this does not choose a physical tolerance.
    """
    points, fixed = arrays['points'], arrays['free_nodes']
    if gradient.shape != points.shape or not np.isfinite(gradient).all():
        raise ValueError('finite current vertex gradient required')
    eta, psi = arrays['parent_nodes'][:, 2], arrays['parent_psi']
    pose = InitialState(arrays['pose_position'], arrays['pose_orientation'])
    mass = parent_mass(arrays)
    scale = np.sqrt(np.diag(mass))
    if not np.isfinite(scale).all() or np.any(scale <= 0):
        raise ValueError('unsupported parent coefficient')
    equilibrated = mass/scale[:, None]/scale[None, :]
    condition = np.linalg.cond(equilibrated, 1)
    if not np.isfinite(condition) or condition*np.finfo(float).eps > 1e-8:
        raise ValueError('parent mass exceeds roundoff budget')
    factor = cho_factor(equilibrated)
    def solve(load):
        return cho_solve(factor, load/scale)/scale
    eta_rate = solve(field['parent_map'].T@field['flux'])
    geometry_evaluations = 0
    def checked(current_pose, height):
        nonlocal geometry_evaluations
        geometry_evaluations += 1
        current = geometry(current_pose, height)
        for key in ('faces', 'free_nodes', 'free_faces', 'body_faces', 'wall_faces', 'bottom_faces', 'node_chart_key'):
            if not np.array_equal(current[key], arrays[key]):
                raise ValueError('parent derivative leaves the owned chart/basis')
        return current
    for halving in range(13):
        dt = 1e-5/2**halving
        try:
            timed = []
            for offset in (dt, -dt, dt/2, -dt/2):
                orientation = (Rotation.from_rotvec(offset*twist[3:])
                    *Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
                timed.append(checked(InitialState(pose.position+offset*twist[:3], orientation), eta+offset*eta_rate))
        except InvalidTemplate:
            continue
        break
    else:
        raise ValueError('no valid owned mass-transport stencil')
    time_mass = [parent_mass(current) for current in timed]
    cdot = [(time_mass[0]-time_mass[1])/(2*dt), (time_mass[2]-time_mass[3])/dt]
    covector = port_transport_covector(arrays, field['potential'], values, time)
    count = len(eta)
    names = ('kinetic', 'body_load', 'incident_load', 'moving_trace', 'mass_transport', 'port_work', 'gravity')
    terms = {name: np.zeros((2, count)) for name in names}
    steps = []
    for index in range(count):
        direction = np.eye(count)[index]
        for halving in range(13):
            h = 1e-5/2**halving
            try:
                moved = [checked(pose, eta+offset*direction) for offset in (h, -h, h/2, -h/2)]
            except InvalidTemplate:
                continue
            break
        else:
            raise ValueError('no valid owned parent-height stencil')
        steps.append(h)
        loads = [boundary_data(current, values, twist, time=time) for current in moved]
        masses = [parent_mass(current) for current in moved]
        for level in range(2):
            step = h/2**level
            plus, minus = moved[2*level:2*level+2]
            bp, bm = loads[2*level:2*level+2]
            moving = (plus['points']-minus['points'])/(2*step)
            dbody = ((bp['body_load']-bm['body_load'])@twist)/(2*step)
            dincident = (bp['right']-bm['right'])/(2*step)-dbody
            dtrace = ((bp['parent_map']-bm['parent_map'])@psi)/(2*step)
            dc = (masses[2*level]-masses[2*level+1])/(2*step)
            terms['kinetic'][level, index] = -np.sum(gradient*moving)
            terms['body_load'][level, index] = WATER_DENSITY*(field['potential']@dbody)
            terms['incident_load'][level, index] = WATER_DENSITY*(field['potential']@dincident)
            terms['moving_trace'][level, index] = -WATER_DENSITY*(field['flux']@dtrace)
            terms['gravity'][level, index] = -WATER_DENSITY*GRAVITY*(mass@eta)[index]
            terms['mass_transport'][level, index] = WATER_DENSITY*(psi@dc@eta_rate-(cdot[level].T@psi)[index])
            terms['port_work'][level, index] = -WATER_DENSITY*(field['potential']@dincident+np.sum(covector*moving))
    raw = sum(value for name, value in terms.items() if name != 'port_work')
    corrected = raw+terms['port_work']
    return dict(eta_rate=eta_rate, psi_rate=np.stack([solve(load/WATER_DENSITY) for load in corrected]),
        parent_mass=mass, mass_time_derivative=np.stack(cdot), raw_action_weak=raw,
        corrected_action_weak=corrected, terms=terms, height_steps_m=np.array(steps),
        transport_step_seconds=dt, geometry_evaluations=geometry_evaluations,
        parent_mass_condition=float(condition), port_seam_weak=port_seam_load(arrays, field, values, time),
        port_seam_weak_order6=port_seam_load(arrays, field, values, time, order=6))
