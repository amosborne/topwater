"""Outward-interval rejection of exterior plane triples before exact construction.

Only rejection is decided here. Every survivor still uses the existing rational
intersection, predicates and owner registration. Arithmetic is CPU float64 with
outward nextafter rounding; nonfinite/uncertain results fall through.
"""
from fractions import Fraction

import numpy as np


def _down(value):
    return np.nextafter(value, -np.inf)


def _up(value):
    return np.nextafter(value, np.inf)


def _add(a, b):
    return _down(a[0]+b[0]), _up(a[1]+b[1])


def _sub(a, b):
    return _down(a[0]-b[1]), _up(a[1]-b[0])


def _mul(a, b):
    products = np.stack((a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]))
    return _down(np.min(products, axis=0)), _up(np.max(products, axis=0))


def _component(value, index):
    return value[0][..., index], value[1][..., index]


def _cross(a, b):
    components = [_sub(_mul(_component(a, i), _component(b, j)),
                       _mul(_component(a, j), _component(b, i)))
                  for i, j in ((1, 2), (2, 0), (0, 1))]
    return tuple(np.stack([entry[side] for entry in components], axis=-1) for side in (0, 1))


def _dot(a, b):
    products = [_mul(_component(a, i), _component(b, i)) for i in range(3)]
    return _add(_add(products[0], products[1]), products[2])


def _coefficient_bounds(planes):
    lower, upper = np.empty((len(planes), 4)), np.empty((len(planes), 4))
    for i, plane in enumerate(planes):
        for j, coefficient in enumerate(plane):
            exact = Fraction(coefficient)
            try:
                value = float(exact)
            except OverflowError:
                lower[i, j], upper[i, j] = -np.inf, np.inf
                continue
            if not np.isfinite(value):
                lower[i, j], upper[i, j] = -np.inf, np.inf
                continue
            represented = Fraction.from_float(value)
            lower[i, j] = _down(value) if represented > exact else value
            upper[i, j] = _up(value) if represented < exact else value
    return lower, upper


class IntersectionRejectionFilter:
    """A cache for one already-registered, frozen canonical plane arrangement."""
    def __init__(self, planes, stats, *, max_entries=50000):
        if isinstance(max_entries, bool) or not isinstance(max_entries, int) or max_entries < 0:
            raise ValueError('nonnegative integer interval-cache cap required')
        self.bounds = _coefficient_bounds(planes)
        self.plane_count = len(planes)
        self.homogeneous = {}
        self.stats = stats
        self.max_entries = max_entries
        for name in ('cache_peak_entries', 'cache_peak_numeric_bytes', 'cache_clears', 'oversize_batches'):
            stats.setdefault(name, 0)
        stats['cache_max_entries'] = max_entries
        stats['coefficient_numeric_bytes'] = sum(value.nbytes for value in self.bounds)
        stats['cache_storage_scope'] = ('D/N float64 endpoints only; excludes Python overhead, '
                                        'coefficient bounds and temporary arrays')
        self._cache_stats()

    def _cache_stats(self):
        entries = len(self.homogeneous)
        # Each entry stores two determinant and six numerator endpoints.
        self.stats['cache_entries'] = entries
        self.stats['cache_numeric_bytes'] = entries*8*np.dtype(np.float64).itemsize
        self.stats['cache_peak_entries'] = max(self.stats['cache_peak_entries'], entries)
        self.stats['cache_peak_numeric_bytes'] = max(self.stats['cache_peak_numeric_bytes'],
                                                    self.stats['cache_numeric_bytes'])

    def _construct_intervals(self, triples):
        ids = np.asarray(triples, dtype=int)
        lower, upper = self.bounds
        planes = lower[ids], upper[ids]
        a = planes[0][:, 0, :3], planes[1][:, 0, :3]
        b = planes[0][:, 1, :3], planes[1][:, 1, :3]
        c = planes[0][:, 2, :3], planes[1][:, 2, :3]
        bc, ca, ab = _cross(b, c), _cross(c, a), _cross(a, b)
        determinant = _dot(a, bc)
        terms = []
        for row, cross in enumerate((bc, ca, ab)):
            right = planes[0][:, row, 3, None], planes[1][:, row, 3, None]
            terms.append(_mul(right, cross))
        numerator = _add(_add(terms[0], terms[1]), terms[2])
        for i, triple in enumerate(triples):
            self.homogeneous[triple] = (determinant[0][i], determinant[1][i],
                                       numerator[0][i].copy(), numerator[1][i].copy())
        self.stats['interval_constructions'] += len(triples)
        self._cache_stats()

    def reject(self, triples, constraints):
        """Return true only for a certified strict halfspace violation."""
        constraints = tuple(constraints)
        if any(sign not in (-1, 1) for _, sign in constraints):
            raise ValueError('halfspace constraint signs must be -1 or +1')
        self.stats['calls'] += 1
        self.stats['candidates'] += len(triples)
        if not triples:
            return np.zeros(0, dtype=bool)
        if len(triples) > self.max_entries or not constraints:
            if len(triples) > self.max_entries:
                self.stats['oversize_batches'] += 1
            self.stats['exact_fallback'] += len(triples)
            return np.zeros(len(triples), dtype=bool)
        missing = list(dict.fromkeys(triple for triple in triples if triple not in self.homogeneous))
        if len(self.homogeneous)+len(missing) > self.max_entries:
            self.homogeneous.clear()
            self.stats['cache_clears'] += 1
            self._cache_stats()
            missing = list(dict.fromkeys(triples))
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            if missing:
                self._construct_intervals(missing)
            cached = [self.homogeneous[triple] for triple in triples]
            determinant = (np.array([entry[0] for entry in cached])[:, None],
                           np.array([entry[1] for entry in cached])[:, None])
            numerator = (np.stack([entry[2] for entry in cached])[:, None],
                         np.stack([entry[3] for entry in cached])[:, None])
            plane_ids, signs = np.asarray(constraints, dtype=int).T
            lower, upper = self.bounds
            normal = lower[None, plane_ids, :3], upper[None, plane_ids, :3]
            offset = lower[None, plane_ids, 3], upper[None, plane_ids, 3]
            residual = _sub(_dot(normal, numerator), _mul(offset, determinant))
            # sign(n.x-d) = sign(n.N-d*D)*sign(D); zero D never rejects.
            positive_d, negative_d = determinant[0] > 0, determinant[1] < 0
            positive_r, negative_r = residual[0] > 0, residual[1] < 0
            outside = np.where(signs[None] > 0,
                (positive_d & positive_r) | (negative_d & negative_r),
                (positive_d & negative_r) | (negative_d & positive_r))
            # Invalid interval endpoints can never certify a rejection.
            valid = (np.isfinite(determinant[0]) & np.isfinite(determinant[1])
                     & np.isfinite(residual[0]) & np.isfinite(residual[1]))
            rejected = np.any(outside & valid, axis=1)
        self.stats['rejected'] += int(rejected.sum())
        self.stats['exact_fallback'] += int((~rejected).sum())
        return rejected
