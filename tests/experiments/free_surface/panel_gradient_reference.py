"""Slow extended-precision reference for rejected bulk panel gradients.

Shared precision recovery for experimental GPU gradients and CPU verification.
The public production evaluator keeps its failure guard. Recompute flagged
entries at two decimal precisions; never accept the original
uncertified result or a fixed quadrature rule near a singularity.
The algebra/logarithms use Decimal; atan2 and pi remain binary64. Agreement of
the two precisions checks cancellation, not a certified interval error bound.
"""
from decimal import Decimal, localcontext
from math import atan2, pi

import numpy as np

from topwater.numerics import LaplaceControls
from topwater.rankine import PanelInfluence, QuadratureError, panel_influence
from topwater._rankine_kernel import analytic_block


def extended(target, triangle, sign, precision):
    def dot(a, b):
        return sum(x*y for x, y in zip(a, b))

    def cross(a, b):
        return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

    def sub(a, b):
        return [x-y for x, y in zip(a, b)]

    with localcontext() as context:
        context.prec = precision
        x = [Decimal(float(v)) for v in target]
        vertices = [[Decimal(float(v)) for v in row] for row in triangle]
        edges = [sub(vertices[(j+1) % 3], vertices[j]) for j in range(3)]
        raw = cross(edges[0], [-v for v in edges[2]])
        area2 = dot(raw, raw).sqrt()
        normal = [v/area2 for v in raw]
        r = [sub(x, v) for v in vertices]
        length = [dot(v, v).sqrt() for v in r]
        height = dot(r[0], normal)
        outward, logs, terms = [], [], Decimal(0)
        for j, edge in enumerate(edges):
            size = dot(edge, edge).sqrt()
            out = cross([v/size for v in edge], normal)
            total = length[j]+length[(j+1) % 3]
            logarithm = ((total+size)/(total-size)).ln()
            terms -= dot(r[j], out)*logarithm
            outward.append(out)
            logs.append(logarithm)
        p01, p12, p20 = dot(r[0], r[1]), dot(r[1], r[2]), dot(r[2], r[0])
        l0, l1, l2 = length
        determinant = area2*height
        denominator = l0*l1*l2+p01*l2+p12*l0+p20*l1
        scale = max(abs(determinant), abs(denominator))
        angle = Decimal(2*atan2(float(determinant/scale), float(denominator/scale)))
        factor = Decimal(4)*Decimal(pi)
        single = (terms-height*angle)/factor
        double = Decimal(int(sign))*angle/factor
        gs, gk = [], []
        for j in range(3):
            gs.append((-sum(logs[e]*outward[e][j] for e in range(3))-angle*normal[j])/factor)
            gd = ((l1*l2+p12)*r[0][j]/l0+(l0*l2+p20)*r[1][j]/l1
                  +(l0*l1+p01)*r[2][j]/l2+l2*(r[0][j]+r[1][j])
                  +l0*(r[1][j]+r[2][j])+l1*(r[2][j]+r[0][j]))
            gt = area2*normal[j]
            gk.append(2*Decimal(int(sign))*(denominator*gt-determinant*gd)
                      /((denominator**2+determinant**2)*factor))
        values = np.array([float(v) for v in [single, double, *gs, *gk]])
        rounding = 64*np.finfo(float).eps
        errors = [rounding*(abs(float(single))+float(abs(height*angle)/factor)),
                  rounding*abs(float(double))]
        errors.extend(rounding*float((sum(abs(logs[e]*outward[e][j]) for e in range(3))
                       +abs(angle*normal[j]))/factor) for j in range(3))
        errors.extend(rounding*abs(float(value)) for value in gk)
        return values, np.array(errors)


def checked_influence(targets, panels):
    try:
        return panel_influence(targets, panels, gradients=True)
    except QuadratureError:
        pass
    controls = LaplaceControls()
    shape = (len(targets), len(panels))
    single, double = np.empty(shape), np.empty(shape)
    gs, gk = np.empty((*shape, 3)), np.empty((*shape, 3))
    unsafe, bad = analytic_block(targets, panels.triangles, panels.geometric_normals,
        panels.edge_outward, panels.tangents, panels.edge_lengths, panels.areas,
        panels.diameters, panels.normal_sign, np.full(len(targets), -1, dtype=np.int64),
        True, controls.quadrature_atol, controls.quadrature_rtol, single, double, gs, gk)
    if bad:
        raise ValueError('reference requires strict bulk targets, away from edges/traces')
    for row, col in zip(*np.nonzero(unsafe)):
        a, _ = extended(targets[row], panels.triangles[col], panels.normal_sign[col], 64)
        b, _ = extended(targets[row], panels.triangles[col], panels.normal_sign[col], 96)
        scale = np.array([1/panels.diameters[col], 1., 1., 1., 1.,
                          panels.diameters[col], panels.diameters[col], panels.diameters[col]])
        tolerance = controls.quadrature_atol+controls.quadrature_rtol*np.abs(b*scale)
        if not np.all(np.isfinite(b)) or np.any(np.abs((a-b)*scale) > tolerance/10):
            raise QuadratureError('extended reference failed precision agreement')
        single[row, col], double[row, col] = b[:2]
        gs[row, col], gk[row, col] = b[2:5], b[5:]
    return PanelInfluence(single, double, gs, gk, int(unsafe.sum()), 0)
