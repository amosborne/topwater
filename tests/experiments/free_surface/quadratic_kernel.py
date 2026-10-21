"""Flat-panel P2/P1 Galerkin source moments and fused target integration.

Analytic moments are checked independently; positive far-source quadrature
avoids polynomial cancellation. No certified target-quadrature control or
curved-panel integration is implemented here. Normals are geometric outward
normals of the closed bounded fluid domain; arbitrary normal_sign is unsupported.
"""

from math import atan2, log1p, pi, sqrt
import numpy as np
from numba import njit


@njit(cache=True)
def source(
    x,
    triangle,
    normal,
    outward,
    tangent,
    length,
    area,
    center,
    gradient,
    self_pair,
    far_points,
    far_weights,
    far_shape,
    s,
    k,
    r,
    distance,
):
    height = 0.0
    for axis in range(3):
        height += (x[axis] - triangle[0, axis]) * normal[axis]
    if self_pair:
        height = 0.0
    # At large separation, positive source quadrature avoids cancellation in
    # polynomial moments. Degree-14 quadrature integrates basis times Taylor
    # polynomials; this is not a complete target quadrature error controller.
    radius2 = 0.0
    for axis in range(3):
        radius2 += (x[axis] - center[axis]) ** 2
    diameter = max(length[0], length[1], length[2])
    if radius2 > (8 * diameter) ** 2:
        s[:] = 0.0
        k[:] = 0.0
        for q in range(len(far_weights)):
            d2 = 0.0
            for axis in range(3):
                d2 += (x[axis] - far_points[q, axis]) ** 2
            value = area * far_weights[q] / (4 * pi * sqrt(d2))
            for j in range(3):
                s[j] += value * far_shape[q, j]
            value *= height / d2
            for j in range(3):
                k[j] += value * far_shape[q, j] * (2 * far_shape[q, j] - 1)
            k[3] += value * 4 * far_shape[q, 0] * far_shape[q, 1]
            k[4] += value * 4 * far_shape[q, 1] * far_shape[q, 2]
            k[5] += value * 4 * far_shape[q, 2] * far_shape[q, 0]
        return
    for edge in range(3):
        d2 = 0.0
        for axis in range(3):
            r[edge, axis] = x[axis] - triangle[edge, axis]
            d2 += r[edge, axis] ** 2
        distance[edge] = sqrt(d2)
    angle = 0.0
    if height != 0:
        denominator = distance[0] * distance[1] * distance[2]
        for edge in range(3):
            dot = 0.0
            for axis in range(3):
                dot += r[edge, axis] * r[(edge + 1) % 3, axis]
            denominator += dot * distance[(edge + 2) % 3]
        angle = 2 * atan2(2 * area * height, denominator)
    terms = 0.0
    bx = 0.0
    by = 0.0
    ex = 0.0
    ey = 0.0
    exx = 0.0
    exy = 0.0
    eyy = 0.0
    gx2 = 0.0
    gy2 = 0.0
    gxy = 0.0
    px = 0.0
    py = 0.0
    for axis in range(3):
        gx = gradient[1, axis]
        gy = gradient[2, axis]
        gx2 += gx * gx
        gy2 += gy * gy
        gxy += gx * gy
        offset = x[axis] - center[axis] - height * normal[axis]
        px += offset * gx
        py += offset * gy
    for edge in range(3):
        d = 0.0
        along = 0.0
        gm = 0.0
        gn = 0.0
        tx = 0.0
        ty = 0.0
        for axis in range(3):
            d -= r[edge, axis] * outward[edge, axis]
            along += r[edge, axis] * tangent[edge, axis]
            gm += gradient[1, axis] * outward[edge, axis]
            gn += gradient[2, axis] * outward[edge, axis]
            tx += gradient[1, axis] * tangent[edge, axis]
            ty += gradient[2, axis] * tangent[edge, axis]
        end = length[edge] - along
        perp = d * d + height * height
        denominator = (
            perp / (distance[edge] + abs(along))
            + perp / (distance[(edge + 1) % 3] + abs(end))
            + 2 * (max(-along, 0.0) + max(along - length[edge], 0.0))
        )
        i0 = log1p(2 * length[edge] / denominator)
        i1 = distance[(edge + 1) % 3] - distance[edge]
        e10 = d * gm * i0 + tx * i1
        e01 = d * gn * i0 + ty * i1
        terms += d * i0
        bx += d * e10
        by += d * e01
        ex -= gm * i0
        ey -= gn * i0
        exx -= gm * e10
        exy -= gm * e01
        eyy -= gn * e01
    s0 = terms - height * angle
    sx = (bx - height * height * ex) / 2
    sy = (by - height * height * ey) / 2
    s[1] = (sx + (px + 1 / 3) * s0) / (4 * pi)
    s[2] = (sy + (py + 1 / 3) * s0) / (4 * pi)
    s[0] = s0 / (4 * pi) - s[1] - s[2]
    if height == 0:
        k[:] = 0.0
        return
    kx = height * ex
    ky = height * ey
    kxx = height * (gx2 * s0 + exx)
    kxy = height * (gxy * s0 + exy)
    kyy = height * (gy2 * s0 + eyy)
    a, b = px + 1 / 3, py + 1 / 3
    l1 = kx + a * angle
    l2 = ky + b * angle
    l0 = angle - l1 - l2
    l11 = kxx + 2 * a * kx + a * a * angle
    l22 = kyy + 2 * b * ky + b * b * angle
    l12 = kxy + b * kx + a * ky + a * b * angle
    l00 = angle - 2 * l1 - 2 * l2 + l11 + 2 * l12 + l22
    k[0] = (2 * l00 - l0) / (4 * pi)
    k[1] = (2 * l11 - l1) / (4 * pi)
    k[2] = (2 * l22 - l2) / (4 * pi)
    k[3] = 4 * (l1 - l11 - l12) / (4 * pi)
    k[4] = 4 * l12 / (4 * pi)
    k[5] = 4 * (l2 - l22 - l12) / (4 * pi)


@njit(cache=True)
def assemble_pairs(
    triangles,
    normals,
    outward,
    tangents,
    lengths,
    areas,
    centers,
    gradient,
    dofs,
    points,
    weights,
    bary,
    far_points,
    far_weights,
    far_shape,
    nv,
):
    nf = len(triangles)
    v = np.zeros((3 * nf, 3 * nf))
    h = np.zeros((3 * nf, nv))
    s = np.empty(3)
    k = np.empty(6)
    r = np.empty((3, 3))
    distance = np.empty(3)
    for target in range(nf):
        for other in range(nf):
            for q in range(len(weights)):
                source(
                    points[target, q],
                    triangles[other],
                    normals[other],
                    outward[other],
                    tangents[other],
                    lengths[other],
                    areas[other],
                    centers[other],
                    gradient[other],
                    target == other,
                    far_points[other],
                    far_weights,
                    far_shape,
                    s,
                    k,
                    r,
                    distance,
                )
                weight = areas[target] * weights[q]
                for a in range(3):
                    w = weight * bary[q, a]
                    for b in range(3):
                        v[3 * target + a, 3 * other + b] += w * s[b]
                    for j in range(6):
                        h[3 * target + a, dofs[other, j]] += w * k[j]
        # Exact half mixed mass: integral(lambda_a * P2_basis_j)/2.
        for a in range(3):
            for j in range(3):
                h[3 * target + a, dofs[target, j]] += areas[target] * (
                    1 / 60 if a == j else -1 / 120
                )
            for j in range(3):
                first, second = j, (j + 1) % 3
                h[3 * target + a, dofs[target, 3 + j]] += areas[target] * (
                    1 / 15 if a == first or a == second else 1 / 30
                )
    return v, h
