"""Bounded common-z=0 scalar S reference, including edges and vertices.

No gradient, double layer, jump, owner-guard relaxation, or certified interval.
"""
from decimal import Decimal, localcontext
from math import asinh, hypot, log, log1p, pi, sqrt
from time import perf_counter
import numpy as np
from numba import njit
from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels, QuadratureError


@njit(cache=True)
def _asinh_ratio(a, b):
    if a <= b:
        return asinh(a/b)
    ratio = b/a
    return log(a)-log(b)+log1p(sqrt(1+ratio*ratio))


@njit(cache=True)
def _evaluate(points, triangles, tangent, outward, lengths, diameter, atol, rtol):
    m, n = len(points), len(triangles)
    result, error = np.empty((m, n)), np.empty((m, n))
    unsafe = np.zeros((m, n), dtype=np.uint8)
    eps = np.finfo(np.float64).eps
    zero_products = 0
    for i in range(m):
        for j in range(n):
            total, absolute, projection_error = 0., 0., 0.
            for edge in range(3):
                rx = points[i, 0]-triangles[j, edge, 0]
                ry = points[i, 1]-triangles[j, edge, 1]
                d = -rx*outward[j, edge, 0]-ry*outward[j, edge, 1]
                if d == 0:
                    # d*log(...) has limit zero, also at edge endpoints.
                    zero_products += 1
                    # A rounded-zero projection is still subject to the
                    # original-float Decimal recovery when error is material.
                    projection_error += abs(rx)+abs(ry)
                    continue
                a = rx*tangent[j, edge, 0]+ry*tangent[j, edge, 1]
                length, b = lengths[j, edge], abs(d)
                if 0 <= a <= length:
                    edge_log = _asinh_ratio(a, b)+_asinh_ratio(length-a, b)
                else:
                    r0, r1 = hypot(a, b), hypot(length-a, b)
                    denominator = (b*(b/r0)/(1+abs(a)/r0)
                        +b*(b/r1)/(1+abs(length-a)/r1)
                        +2*(max(-a, 0.)+max(a-length, 0.)))
                    if length <= denominator/2:
                        edge_log = log1p(2*(length/denominator))
                    else:
                        edge_log = log(length)-log(denominator)+log(2+denominator/length)
                term = d*edge_log
                total += term
                absolute += abs(term)
                projection_error += (abs(rx)+abs(ry))*(1+abs(edge_log))
            value = total/(4*pi)
            estimate = 64*eps*(absolute+projection_error)/(4*pi)
            result[i, j], error[i, j] = value, estimate
            if (not np.isfinite(value) or value <= 0 or not np.isfinite(estimate)
                    or estimate > atol*diameter[j]+rtol*abs(value)):
                unsafe[i, j] = 1
    return result, error, unsafe, zero_products


def _decimal_pair(point, triangle, precision):
    """Recompute original-float geometry; no cached binary64 normals."""
    with localcontext() as context:
        context.prec = precision
        vertices = [[Decimal.from_float(float(c)) for c in row[:2]] for row in triangle]
        p = [Decimal.from_float(float(c)) for c in point]
        e1 = [vertices[1][k]-vertices[0][k] for k in range(2)]
        e2 = [vertices[2][k]-vertices[0][k] for k in range(2)]
        area2 = e1[0]*e2[1]-e1[1]*e2[0]
        sign = Decimal(1) if area2 > 0 else Decimal(-1)
        total = Decimal(0)
        for e in range(3):
            a0, a1 = vertices[e], vertices[(e+1)%3]
            edge = [a1[k]-a0[k] for k in range(2)]
            r = [p[k]-a0[k] for k in range(2)]
            length = sum(c*c for c in edge).sqrt()
            a = sum(r[k]*edge[k] for k in range(2))/length
            d = sign*(edge[0]*r[1]-edge[1]*r[0])/length
            if d == 0:
                continue
            b = abs(d)
            if 0 <= a <= length:
                u, v = a/b, (length-a)/b
                edge_log = (u+(1+u*u).sqrt()).ln()+(v+(1+v*v).sqrt()).ln()
            else:
                r0, r1 = (a*a+b*b).sqrt(), ((length-a)**2+b*b).sqrt()
                denominator = b*b/(r0+abs(a))+b*b/(r1+abs(length-a))+2*(max(-a, 0)+max(a-length, 0))
                edge_log = (1+2*length/denominator).ln()
            total += d*edge_log
        return float(total/(4*Decimal.from_float(pi)))


def plane_single_layer(points_xy, panels, controls=LaplaceControls(), *, maximum_fallbacks=256):
    """Return positive-Rankine S and diagnostics for arbitrary plane targets.

    Error requests apply to S/panel_diameter. Decimal precision agreement is
    a bounded reference recovery, not a certified integration-error interval.
    Storage bound covers owned targets, scalar outputs/errors/flags and bounded
    fallback indices; excludes caller-owned panels and Numba compile workspace.
    """
    started = perf_counter()
    raw = np.asarray(points_xy)
    if (not isinstance(panels, TrianglePanels) or not isinstance(controls, LaplaceControls)
            or raw.ndim != 2 or raw.shape[1] != 2 or np.iscomplexobj(raw)
            or not np.isfinite(raw).all() or np.any(panels.triangles[:, :, 2] != 0)
            or isinstance(maximum_fallbacks, bool) or not isinstance(maximum_fallbacks, int)
            or maximum_fallbacks < 0):
        raise ValueError('finite real plane targets, z=0 panels and bounded fallback count required')
    m, n = len(raw), len(panels)
    bound = 17*m*n+16*m+8*min(maximum_fallbacks, m*n)
    if bound > controls.max_operator_bytes:
        raise MemoryError('plane scalar operator/scratch exceeds byte budget')
    points = np.array(raw, dtype=float, order='C', copy=True)
    values, errors, unsafe, zero_products = _evaluate(points, panels.triangles, panels.tangents,
        panels.edge_outward, panels.edge_lengths, panels.diameters,
        controls.quadrature_atol, controls.quadrature_rtol)
    count = int(unsafe.sum())
    if count > maximum_fallbacks:
        raise QuadratureError(f'plane scalar fallback budget exceeded: {count}>{maximum_fallbacks}')
    agreement = 0.
    for index in np.flatnonzero(unsafe):
        i, j = divmod(int(index), n)
        a = _decimal_pair(points[i], panels.triangles[j], 64)
        b = _decimal_pair(points[i], panels.triangles[j], 96)
        estimate = abs(a-b)+64*np.finfo(float).eps*abs(b)
        tolerance = controls.quadrature_atol*panels.diameters[j]+controls.quadrature_rtol*abs(b)
        if not np.isfinite([a, b, estimate]).all() or b <= 0 or estimate > tolerance:
            raise QuadratureError('plane scalar Decimal precision/error request unresolved')
        values[i, j], errors[i, j] = b, estimate
        agreement = max(agreement, abs(a-b))
    return values, dict(seconds=perf_counter()-started, source_pairs=m*n,
        fallback_pairs=count, maximum_fallbacks=maximum_fallbacks,
        exact_zero_distance_products=int(zero_products),
        requested_rtol=controls.quadrature_rtol,
        requested_scaled_atol_min=float(controls.quadrature_atol*panels.diameters.min()),
        requested_scaled_atol_max=float(controls.quadrature_atol*panels.diameters.max()),
        maximum_precision_difference=agreement, maximum_error_estimate=float(errors.max(initial=0)),
        output_bytes=values.nbytes, owned_array_bytes_bound=bound,
        certified_error_bound=False, scope='Common-plane scalar only; exact d=0 limit; Decimal agreement is not an interval certificate')
