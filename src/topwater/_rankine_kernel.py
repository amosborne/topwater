"""Fused constant-triangle integration; scalar scratch, no pairwise vectors.

No fast-math reassociation: the cancellation guards rely on these operations.
Compilation is cached; the benchmark reports cold and repeated calls separately.
"""

from math import atan2, log1p, pi, sqrt

from numba import njit
import numpy as np


@njit(cache=True, nogil=True, error_model="numpy")
def analytic_block(targets, triangles, normals, outward, tangents, lengths,
                   areas, diameters, signs, self_indices, gradients, atol, rtol,
                   single, double, gs, gk):
    m, n = len(targets), len(triangles)
    unsafe = np.zeros((m, n), dtype=np.bool_)
    bad = 0
    eps, factor = np.finfo(np.float64).eps, 4*pi
    # Reused fixed-size scratch. Only output arrays grow with pair count.
    r, distance = np.empty((3, 3)), np.empty(3)
    logarithm = np.empty(3)
    for row in range(m):
        for col in range(n):
            self_pair = self_indices[row] == col
            vertex = False
            for edge in range(3):
                squared = 0.
                for j in range(3):
                    r[edge, j] = targets[row, j]-triangles[col, edge, j]
                    squared += r[edge, j]**2
                distance[edge] = sqrt(squared)
                vertex |= squared == 0
            if vertex:
                bad |= 1
                continue
            height = 0.
            for j in range(3):
                height += r[0, j]*normals[col, j]
            if self_pair:
                height = 0.
            near_plane = abs(height) <= 64*eps*diameters[col]
            inside = True
            terms, absolute_terms = 0., 0.
            for edge in range(3):
                d, a = 0., 0.
                for j in range(3):
                    d -= r[edge, j]*outward[col, edge, j]
                    a += r[edge, j]*tangents[col, edge, j]
                inside &= d > 0
                length = lengths[col, edge]
                perpendicular2 = d*d+height*height
                if near_plane and perpendicular2 <= (64*eps*length)**2 and 0 <= a <= length:
                    bad |= 2
                end = length-a
                denominator = (perpendicular2/(distance[edge]+abs(a))
                               + perpendicular2/(distance[(edge+1) % 3]+abs(end))
                               + 2*(max(-a, 0.)+max(a-length, 0.)))
                logarithm[edge] = log1p(2*length/denominator)
                term = d*logarithm[edge]
                terms += term
                absolute_terms += abs(term)
            if near_plane and inside:
                if not self_pair:
                    bad |= 4
                if gradients:
                    bad |= 8
            p01, p12, p20 = 0., 0., 0.
            for j in range(3):
                p01 += r[0, j]*r[1, j]
                p12 += r[1, j]*r[2, j]
                p20 += r[2, j]*r[0, j]
            l0, l1, l2 = distance[0], distance[1], distance[2]
            determinant = 2*areas[col]*height
            base = l0*l1*l2
            denominator = base+p01*l2+p12*l0+p20*l1
            angle = 0. if self_pair else 2*atan2(determinant, denominator)
            squared_angle = denominator**2+determinant**2
            denominator_error = 32*eps*(base+abs(p01*l2)+abs(p12*l0)+abs(p20*l1))
            angle_error = 0. if self_pair else 2*abs(determinant)*denominator_error/squared_angle
            s = (terms-height*angle)/factor
            k = angle*signs[col]/factor
            single[row, col], double[row, col] = s, k
            s_error = (32*eps*(absolute_terms+abs(height*angle))+abs(height)*angle_error)/factor
            unresolved = (not np.isfinite(s) or s <= 0 or not np.isfinite(k)
                          or s_error > atol*diameters[col]+rtol*abs(s)
                          or angle_error/factor > atol+rtol*abs(k))
            if gradients:
                for j in range(3):
                    edge_sum, edge_abs = 0., 0.
                    for edge in range(3):
                        term = logarithm[edge]*outward[col, edge, j]
                        edge_sum += term
                        edge_abs += abs(term)
                    g_s = (-edge_sum-angle*normals[col, j])/factor
                    grad_d = ((l1*l2+p12)*r[0, j]/l0
                              + (l0*l2+p20)*r[1, j]/l1
                              + (l0*l1+p01)*r[2, j]/l2
                              + l2*(r[0, j]+r[1, j])
                              + l0*(r[1, j]+r[2, j])
                              + l1*(r[2, j]+r[0, j]))
                    grad_t = 2*areas[col]*normals[col, j]
                    g_k = 2*(denominator*grad_t-determinant*grad_d)/squared_angle*signs[col]/factor
                    gs[row, col, j], gk[row, col, j] = g_s, g_k
                    gs_error = (32*eps*(edge_abs+abs(angle*normals[col, j]))
                                + angle_error*abs(normals[col, j]))/factor
                    gk_error = (128*eps*(abs(denominator*grad_t)+abs(determinant*grad_d))
                                / squared_angle/factor)
                    gk_error += (4*denominator_error*abs(grad_t)/squared_angle/factor
                                 + 4*abs(g_k)*abs(denominator)*denominator_error/squared_angle)
                    unresolved |= (not np.isfinite(g_s) or not np.isfinite(g_k)
                                   or gs_error > atol+rtol*abs(g_s)
                                   or gk_error*diameters[col] > atol+rtol*abs(g_k)*diameters[col])
            unsafe[row, col] = unresolved
    return unsafe, bad
