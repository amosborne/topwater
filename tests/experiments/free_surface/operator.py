"""Experimental symmetric P1/P0 Laplace energy on a complete closed boundary.

S = W + H.T V^-1 H, H = K + .5*mixed_mass. W uses surface curls.
Uniform outer quadrature preserves the closed-boundary constant identity.
This is not a tolerance-controlled Galerkin assembler: see README.md.
"""
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.sparse import coo_matrix
from tests.laplace_reference import triangle_samples
from .kernel import assemble_pairs, symmetrize_single_layer


def assemble_energy(faces, panels, order=3):
    points, weights = triangle_samples(panels.triangles, order)
    nf, nv = len(panels), int(faces.max())+1
    opposite = (np.roll(panels.triangles, -2, axis=1)
                - np.roll(panels.triangles, -1, axis=1))
    gradient = (np.cross(panels.geometric_normals[:, None], opposite)
                / (2*panels.areas[:, None, None]))
    v, h = assemble_pairs(
        panels.triangles, panels.geometric_normals, panels.edge_outward,
        panels.tangents, panels.edge_lengths, panels.areas, panels.centroids,
        panels.normal_sign, faces, gradient, points, weights, nv,
    )
    # Average the two quadratures of the mathematically symmetric integral V.
    # Report their discrepancy first; this does not certify integration accuracy.
    asymmetry = float(symmetrize_single_layer(v))
    curl = np.cross(panels.geometric_normals[:, None], gradient)
    w = np.zeros((nv, nv))
    for axis in range(3):
        b = coo_matrix((curl[:, :, axis].ravel(),
                        (np.repeat(np.arange(nf), 3), faces.ravel())),
                       shape=(nf, nv)).tocsr()
        w += b.T@(v@b)
    # V is no longer needed after W. Its transpose is a Fortran-contiguous
    # view of the same owned symmetric buffer, so LAPACK can consume it.
    factor = cholesky(v.T, lower=True, overwrite_a=True)
    # V = L L.T, so H.T V^-1 H = (L^-1 H).T (L^-1 H).
    # Forming this Gram matrix needs only one triangular solve.
    auxiliary = solve_triangular(factor, h, lower=True)
    w += auxiliary.T@auxiliary
    return w, dict(
        v_quadrature_asymmetry=asymmetry,
        h_constant_norm=float(np.max(np.abs(h@np.ones(nv)))),
        operator_bytes=dict(single_layer_factor=v.nbytes, mixed_trace=h.nbytes,
                            energy=w.nbytes),
    )
