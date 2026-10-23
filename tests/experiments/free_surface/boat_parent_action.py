"""Projected parent mass and incident-port work for the restricted action."""
from collections import Counter

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cho_factor, cho_solve

from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from topwater.constants import WATER_DENSITY
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


def factor_parent_mass(mass):
    """One equilibrated parent solve and its roundoff-budget check."""
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
    return solve, float(condition)


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
