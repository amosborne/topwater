"""Shared validation for physical inputs; no unit conversion is implicit."""

import numpy as np
from numpy.typing import ArrayLike, NDArray


def array(value: ArrayLike, shape: tuple[int | None, ...], name: str) -> NDArray:
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numbers")
    result = np.array(raw, dtype=np.float64, copy=True, order="C")
    if result.ndim != len(shape) or any(
        expected is not None and actual != expected
        for actual, expected in zip(result.shape, shape)
    ):
        raise ValueError(f"{name} must have shape {shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    result.setflags(write=False)
    return result


def quaternion(value: ArrayLike) -> NDArray:
    result = array(value, (4,), "orientation")
    norm = np.linalg.norm(result)
    if not np.isclose(norm, 1.0, rtol=0, atol=1e-10):
        raise ValueError("orientation must be a unit quaternion in (w, x, y, z) order")
    # Remove roundoff only; a materially non-unit input is an error.
    result = result / norm
    result.setflags(write=False)
    return result


def indices(value: ArrayLike, size: int, name: str) -> NDArray:
    """Own a unique one-dimensional index list, including an empty list."""
    raw = np.asarray(value)
    if raw.ndim != 1 or (raw.size and raw.dtype.kind not in 'iu'):
        raise ValueError(f"{name} must be a one-dimensional integer index list")
    if np.any(raw < 0) or np.any(raw >= size):
        raise ValueError(f"{name} contains out-of-range indices")
    result = np.array(raw, dtype=np.int64, copy=True)
    if len(np.unique(result)) != len(result):
        raise ValueError(f"{name} contains duplicate indices")
    result.setflags(write=False)
    return result
