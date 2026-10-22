"""Analytic 3D coordinate/conormal controls before any canoe PML solve.

Frequency-domain only. These controls do not establish a discrete wave error,
Capytaine agreement, or a causal nonlinear time-domain boundary.
"""
import numpy as np

from .pml_strip import stretching


def coordinates(points, onset, width, strength, omega):
    """Signed outward stretch in x/y, with an exactly physical central box."""
    points = np.asarray(points)
    onset, width = np.broadcast_arrays(np.asarray(onset, float), np.asarray(width, float))
    if onset.shape != (2,) or width.shape != (2,) or not np.all(np.isfinite(onset)) or np.any(onset < 0):
        raise ValueError('two finite nonnegative onsets and two positive widths required')
    if points.shape[-1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError('finite three-dimensional points required')
    transformed = np.array(points, dtype=complex)
    derivative = np.ones(points.shape[:-1]+(2,), complex)
    for axis in (0, 1):
        distance = np.maximum(abs(points[..., axis])-onset[axis], 0.)
        value, slope, _ = stretching(distance, width[axis], strength, omega)
        # Subtract distance to preserve the input coordinate exactly in the
        # quiet region. Negative-side imaginary coordinates must be negative.
        transformed[..., axis] += np.sign(points[..., axis])*(value-distance)
        derivative[..., axis] = slope
    return transformed, derivative


def coefficients(stretch):
    """Diagonal volume tensor and flat free-surface measure, including corners."""
    sx, sy = np.moveaxis(np.asarray(stretch), -1, 0)
    return np.stack((sy/sx, sx/sy, sx*sy), axis=-1), sx*sy
