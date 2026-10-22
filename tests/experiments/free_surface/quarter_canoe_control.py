"""Coarse exact-mirror parity falsifier; not a Capytaine qualification.

The full operator and full six-mode body loads are assembled independently.
The mirrored quadrature is part of this explicitly mirrored discretization.
Neither the old full-tank TetGen mesh nor nonlinear motion is parity reduced.
"""
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from . import p4_tetra
from .p4_canoe import loads
from .pml_corner import coordinates, coefficients
from .pocket import boundary
from .quarter_canoe import PARITY, mirror


@lru_cache(None)
def rule(order, swapped=False):
    if isinstance(order, bool) or not isinstance(order, int) or order < 5:
        raise ValueError('integer quadrature order at least five required')
    x, w = leggauss(order); x = (x+1)/2; w = w/2
    u, v, t = np.meshgrid(x, x, x, indexing='ij')
    wu, wv, wt = np.meshgrid(w, w, w, indexing='ij')
    xyz = np.column_stack((u.ravel(), ((1-u)*v).ravel(), ((1-u)*(1-v)*t).ravel()))
    bary = np.column_stack((1-xyz.sum(axis=1), xyz))
    if swapped:
        bary[:, :2] = bary[:, 1::-1]
    _, d = p4_tetra.basis(bary)
    return bary, d[:, :, 1:]-d[:, :, :1], (wu*wv*wt*(1-u)**2*(1-v)).ravel()


def scalar_operator(mesh, omega, onset, width, strength, *, order=12, swaps=None, surface_faces=None):
    """Physical complex-stretch weak form K_s - omega² M_s/g.

    Ordinary transpose gives the complex-symmetric form. No conjugation,
    physical sponge Q, symmetrization, or boundary row correction is applied.
    This deliberately straightforward assembler is for coarse controls only.
    """
    if len(mesh.points) > 20000:
        raise ValueError('coarse control exceeds its 20000-node assembly cap')
    if not np.isfinite(omega) or omega <= 0:
        raise ValueError('positive finite frequency required')
    if swaps is None:
        swaps = np.zeros(len(mesh.tets), bool)
    if np.shape(swaps) != (len(mesh.tets),):
        raise ValueError('one quadrature permutation per tetrahedron required')
    matrices = []
    for tet, swapped in zip(mesh.tets, swaps):
        p = mesh.points[tet]; edges = p[1:]-p[:1]
        determinant = np.linalg.det(edges)
        if determinant <= 0:
            raise ValueError('positive physical orientation required')
        bary, derivative, weights = rule(order, bool(swapped))
        gradient = derivative@np.linalg.inv(edges).T
        _, stretch = coordinates(bary@p, onset, width, strength, omega)
        tensor, _ = coefficients(stretch)
        matrices.append(np.einsum('qai,qbi,qi,q->ab', gradient, gradient, tensor, weights)*determinant)
    nodes = mesh.elements
    matrix = coo_matrix((np.asarray(matrices).ravel(),
        (np.repeat(nodes, 35, axis=1).ravel(), np.tile(nodes, (1, 35)).ravel())),
        shape=(len(mesh.points),)*2).tocsr()
    faces = boundary(mesh.points, mesh.tets)
    top = faces[np.all(mesh.points[faces, 2] == 0, axis=1)]
    if surface_faces is not None:
        identities = lambda f: {tuple(sorted(map(int, row))) for row in f}
        if len(surface_faces) != len(top) or identities(surface_faces) != identities(top):
            raise ValueError('mirrored surface quadrature changed the physical top facets')
        # Ordering is inherited from the quarter face, including reflections.
        # Thus even when a ramp breakpoint cuts a face, its integration rule is
        # the exact physical reflection of the quarter rule. Area is unsigned;
        # no surface normal enters this scalar Robin integral.
        top = np.asarray(surface_faces)
    panels = TrianglePanels(mesh.points[top], 1)
    samples, weights = triangle_samples(panels.triangles, order)
    bary, _ = triangle_samples(np.eye(3)[None], order)
    values, _ = p4_tetra.basis(bary[0])
    _, stretch = coordinates(samples, onset, width, strength, omega)
    _, measure = coefficients(stretch)
    local = np.einsum('f,q,qa,qb,fq->fab', panels.areas, weights, values, values, measure)
    nodes = p4_tetra.face_nodes(mesh, top)
    mass = coo_matrix((local.ravel(), (np.repeat(nodes, 15, axis=1).ravel(),
        np.tile(nodes, (1, 15)).ravel())), shape=matrix.shape).tocsr()
    return matrix-omega**2/9.81*mass


def prolongation(quarter, full, maps, signs, parity):
    """Exact parity injection, excluding odd symmetry-plane trace nodes."""
    fixed = np.zeros(len(quarter.points), bool)
    for axis in (0, 1):
        if parity[axis] == -1:
            fixed |= quarter.points[:, axis] == 0
    unknown = np.flatnonzero(~fixed)
    inverse = np.full(len(quarter.points), -1, int); inverse[unknown] = np.arange(len(unknown))
    owner = np.full(len(full.points), -1, int); weight = np.zeros(len(full.points), int)
    for mapping, reflection in zip(maps, signs):
        factor = int(np.prod([parity[a] if reflection[a] < 0 else 1 for a in (0, 1)]))
        for q in unknown:
            f = mapping[q]
            if owner[f] >= 0 and (owner[f] != inverse[q] or weight[f] != factor):
                raise ValueError('inconsistent parity at a shared reflection-plane node')
            owner[f], weight[f] = inverse[q], factor
    rows = np.flatnonzero(owner >= 0)
    return coo_matrix((weight[rows], (rows, owner[rows])),
                      shape=(len(full.points), len(unknown))).tocsr(), unknown


def compare(quarter, body_faces, center, omega, onset, width, strength, *, order=12):
    """All six raw responses and four exact restriction identities.

    Independent full-domain load integration catches an erroneous mode parity,
    torque center or fourfold work factor; it is not copied from quarter B.
    """
    full, maps, swaps, signs = mirror(quarter)
    aq = scalar_operator(quarter, omega, onset, width, strength, order=order)
    qfaces = boundary(quarter.points, quarter.tets)
    qtop = qfaces[np.all(quarter.points[qfaces, 2] == 0, axis=1)]
    mirrored_top = np.concatenate([mapping[qtop] for mapping in maps])
    af = scalar_operator(full, omega, onset, width, strength, order=order, swaps=swaps,
                         surface_faces=mirrored_top)
    bq, _ = loads(quarter, body_faces, center)
    # Fresh outward boundary winding determines reflected body triangles.
    body_keys = {tuple(sorted(map(int, face))) for mapping in maps for face in mapping[body_faces]}
    full_faces = boundary(full.points, full.tets)
    full_body = np.array([face for face in full_faces if tuple(sorted(map(int, face))) in body_keys])
    if len(full_body) != 4*len(body_faces):
        raise ValueError('reflected body facet coverage changed')
    bf, _ = loads(full, full_body, center)
    full_phi = splu(af.tocsc()).solve(bf.astype(complex))
    rf = bf.T@full_phi
    rq = np.zeros((6, 6), complex); diagnostics = []
    for parity in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
        modes = np.flatnonzero(np.all(PARITY == parity, axis=1))
        t, unknown = prolongation(quarter, full, maps, signs, parity)
        local = aq[unknown][:, unknown]
        rhs = bq[unknown][:, modes]
        phi = splu(local.tocsc()).solve(rhs.astype(complex))
        rq[np.ix_(modes, modes)] = 4*rhs.T@phi
        difference = (t.T@af@t-4*local).tocsr()
        matrix_norm = np.sqrt(np.sum(abs(local.data)**2))
        load_error = np.linalg.norm(t.T@bf[:, modes]-4*rhs)/max(np.linalg.norm(4*rhs), np.finfo(float).tiny)
        reconstructed = t@phi
        diagnostics.append(dict(parity=list(parity), modes=modes.tolist(),
            matrix_restriction_relative=float(np.sqrt(np.sum(abs(difference.data)**2))/(4*matrix_norm)),
            load_restriction_relative=float(load_error),
            full_solution_relative=float(np.linalg.norm(reconstructed-full_phi[:, modes])/np.linalg.norm(full_phi[:, modes]))))
    return dict(qualified=False, omega=omega, quarter_nodes=len(quarter.points), full_nodes=len(full.points),
        raw_response_quarter_real=rq.real.tolist(), raw_response_quarter_imag=rq.imag.tolist(),
        raw_response_full_real=rf.real.tolist(), raw_response_full_imag=rf.imag.tolist(),
        full_matrix_relative_difference=float(np.linalg.norm(rq-rf)/np.linalg.norm(rf)),
        full_residual_relative=float(np.linalg.norm(af@full_phi-bf)/np.linalg.norm(bf)),
        parity_controls=diagnostics,
        scope='Coarse exact-mirrored discrete correctness control only. Structurally zero '
              'cross-parity entries are derived, not independently computed quarter entries. '
              'No propagation, mesh convergence, nonlinear or Capytaine qualification.')
