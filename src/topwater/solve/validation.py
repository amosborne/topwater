"""Finite solver arrays and ordered frequency/time grids."""
import numpy as np


def _finite(value, shape, name, *, complex_value=False):
    result = np.asarray(value, dtype=complex if complex_value else float)
    if result.ndim != len(shape) or any(n is not None and n != m for n, m in zip(shape, result.shape)):
        raise ValueError(f'{name} must have shape {shape}')
    if not np.isfinite(result).all():
        raise ValueError(f'{name} must be finite')
    return result

def _grid(values, name, *, zero=False):
    values = _finite(values, (None,), name)
    if len(values) < 2 or np.any(np.diff(values) <= 0) or values[0] < 0 or (not zero and values[0] == 0):
        raise ValueError(f'{name} must increase strictly from a {"nonnegative" if zero else "positive"} value')
    return values
