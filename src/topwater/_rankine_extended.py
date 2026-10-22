"""Rare scalar analytic fallback; never used by the compiled pairwise fast path."""

from decimal import Decimal, localcontext
from math import atan2, pi


def triangle_values(target, triangle, normal, precision):
    """Recompute geometry and logarithms from float inputs at decimal precision.

    Only atan2 and the returned values round to binary64. Re-evaluation at a
    higher precision checks arithmetic cancellation, not spatial quadrature.
    The second result estimates final binary64 rounding, including the angle's
    contribution to the single-layer integral.
    """
    def dot(a, b):
        return sum(x*y for x, y in zip(a, b))

    def subtract(a, b):
        return [x-y for x, y in zip(a, b)]

    def cross(a, b):
        return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
                a[0]*b[1]-a[1]*b[0]]

    with localcontext() as context:
        context.prec = precision
        point = [Decimal(float(x)) for x in target]
        vertices = [[Decimal(float(x)) for x in vertex] for vertex in triangle]
        edges = [subtract(vertices[(i+1) % 3], vertices[i]) for i in range(3)]
        geometric = cross(edges[0], [-x for x in edges[2]])
        twice_area = dot(geometric, geometric).sqrt()
        unit = [x/twice_area for x in geometric]
        sign = 1 if dot(geometric, [Decimal(float(x)) for x in normal]) > 0 else -1
        displacement = [subtract(point, vertex) for vertex in vertices]
        distances = [dot(delta, delta).sqrt() for delta in displacement]
        height = dot(displacement[0], unit)
        terms = Decimal(0)
        for i in range(3):
            length = dot(edges[i], edges[i]).sqrt()
            tangent = [x/length for x in edges[i]]
            outward = cross(tangent, unit)
            offset = -dot(displacement[i], outward)
            total = distances[i]+distances[(i+1) % 3]
            terms += offset*((total+length)/(total-length)).ln()
        determinant = twice_area*height
        denominator = (distances[0]*distances[1]*distances[2]
            + dot(displacement[0], displacement[1])*distances[2]
            + dot(displacement[1], displacement[2])*distances[0]
            + dot(displacement[2], displacement[0])*distances[1])
        scale = max(abs(determinant), abs(denominator))
        angle = 2*atan2(float(determinant/scale), float(denominator/scale))
        factor = Decimal(4)*Decimal(pi)
        single = float((terms-height*Decimal(angle))/factor)
        double = sign*angle/(4*pi)
        # Conservative allowance for atan2 arguments, atan2, pi and float output.
        rounding = 64*2.**-52
        error = (rounding*(abs(single)+float(abs(height*Decimal(angle))/factor)),
                 rounding*abs(double))
        return (single, double), error
