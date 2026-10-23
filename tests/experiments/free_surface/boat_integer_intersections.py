"""Exact integer homogeneous arithmetic for the frozen rational plane owner.

There is no floating predicate or alternate registration/sign cache. Rational
coordinates are constructed only once per intersection for the existing API.
"""
from fractions import Fraction
from math import gcd, lcm


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def _point_homogeneous(point):
    values = tuple(Fraction(value) for value in point)
    denominator = lcm(*(value.denominator for value in values))
    return tuple(value.numerator*(denominator//value.denominator) for value in values), denominator


class IntegerIntersections:
    """One frozen plane arrangement; storage follows the existing exact cap."""
    def __init__(self, planes):
        self.frozen_planes = tuple(tuple(plane) for plane in planes)
        converted = []
        for plane in self.frozen_planes:
            values = tuple(Fraction(value) for value in plane)
            denominator = lcm(*(value.denominator for value in values))
            integers = tuple(value.numerator*(denominator//value.denominator) for value in values)
            divisor = gcd(*integers)
            if not divisor or not any(integers[:3]):
                raise ValueError('nonzero exact plane normal required')
            # Positive scale preserves the original halfspace orientation.
            converted.append(tuple(value//divisor for value in integers))
        self.planes = tuple(converted)
        self.homogeneous = {}

    def check_planes(self, planes):
        if tuple(tuple(plane) for plane in planes) != self.frozen_planes:
            raise ValueError('integer intersection owner requires unchanged frozen planes')

    def intersection(self, triple):
        a, b, c = (self.planes[index] for index in triple)
        crosses = _cross(b[:3], c[:3]), _cross(c[:3], a[:3]), _cross(a[:3], b[:3])
        determinant = _dot(a[:3], crosses[0])
        if not determinant:
            self.homogeneous[triple] = None
            return None
        numerator = tuple(sum(plane[3]*vector[k] for plane, vector in zip((a, b, c), crosses))
                          for k in range(3))
        divisor = gcd(*numerator, determinant)
        if determinant < 0:
            divisor = -divisor
        numerator = tuple(value//divisor for value in numerator)
        determinant //= divisor
        self.homogeneous[triple] = numerator, determinant
        return tuple(Fraction(value, determinant) for value in numerator)

    def sign(self, triple, plane):
        numerator, denominator = self.homogeneous[triple]
        row = self.planes[plane]
        value = _dot(row[:3], numerator)-row[3]*denominator
        return (value > 0)-(value < 0)

    def sign_point(self, point, plane):
        """Preserve the original predicate for a caller-supplied exact point."""
        numerator, denominator = _point_homogeneous(point)
        row = self.planes[plane]
        value = _dot(row[:3], numerator)-row[3]*denominator
        return (value > 0)-(value < 0)

    def signs(self, points):
        """All plane signs; normalize each rational witness exactly once."""
        result = []
        for point in points:
            numerator, denominator = _point_homogeneous(point)
            row = []
            for plane in self.planes:
                value = _dot(plane[:3], numerator)-plane[3]*denominator
                row.append((value > 0)-(value < 0))
            result.append(row)
        return result
