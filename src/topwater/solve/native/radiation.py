"""Full-history radiation with causal, stage-local time quadrature."""
import numpy as np
from numba import njit
from scipy.special import spherical_jn
from topwater.solve.settings import DEFAULTS


def radiation_kernel(omega, damping, *, duration_s, dt_s=DEFAULTS.time_step_s):
    omega, damping = np.asarray(omega), np.asarray(damping)
    if (omega.ndim != 1 or len(omega) < 2 or omega[0] <= 0 or np.any(np.diff(omega) <= 0)
            or damping.shape != (len(omega), 6, 6)
            or not np.isfinite(omega).all() or not np.isfinite(damping).all()):
        raise ValueError('Radiation requires finite increasing frequencies and (frequency,6,6) damping')
    if dt_s <= 0 or duration_s <= 0 or abs(round(duration_s/dt_s)*dt_s-duration_s) > 1e-10:
        raise ValueError('Radiation kernel duration must be a positive integer number of timesteps')
    time = np.arange(round(duration_s/dt_s)+1)*dt_s
    # Integrate the existing piecewise-linear B/omega interpolant exactly.
    # Centered interval moments expressed as spherical Bessel functions avoid
    # cancellation at zero lag (and at very short lags).
    midpoint = (omega[:-1]+omega[1:])/2
    halfwidth = np.diff(omega)/2
    z = time[:, None]*halfwidth
    phase = time[:, None]*midpoint
    j0, j1, j2 = (spherical_jn(order, z) for order in range(3))
    cosine, sine = np.cos(phase), np.sin(phase)
    even = halfwidth*(2*midpoint*j0*cosine-2*halfwidth*j1*sine)
    odd = halfwidth*(halfwidth*(2*j0/3-4*j2/3)*cosine-2*midpoint*j1*sine)
    values = (damping/omega[:, None, None]).reshape(len(omega), 36)
    from threadpoolctl import threadpool_limits
    with threadpool_limits(1):
        kernel = ((even-odd) @ values[:-1]+(even+odd) @ values[1:])/np.pi
    kernel = kernel.reshape(len(time), 6, 6)
    return time, kernel


def stage_weights(kernel, dt_s, fraction):
    """Linear operator on accepted history plus the provisional stage velocity.

Integrate on the union of shifted history knots, zero lag, and the memory
cutoff. Interpolate at a partial oldest interval; never extend past the cutoff.
"""
    n = len(kernel)-1
    shift = fraction*dt_s
    end = n*dt_s
    lag = np.unique(np.r_[0., (fraction+np.arange(n+1))*dt_s, end])
    lag = lag[(lag >= 0) & (lag <= end)]
    quadrature = np.r_[np.diff(lag), 0.]/2 + np.r_[0., np.diff(lag)]/2
    sampled = np.empty((len(lag), 6, 6))
    for i in range(6):
        for j in range(6):
            sampled[:, i, j] = np.interp(lag, np.arange(n+1)*dt_s, kernel[:, i, j])
    weights = np.zeros((n+2, 6, 6))
    current = np.zeros((6, 6))
    for age, value in zip(lag, quadrature[:, None, None]*sampled):
        if fraction > 0 and age == 0:
            current += value
        else:
            index = max(0., (age-shift)/dt_s)
            left = int(np.floor(index))
            alpha = index-left
            weights[left] += (1-alpha)*value
            weights[left+1] += alpha*value
    return weights, current


@njit(cache=True)
def history_force(history, step, weights):
    result = np.zeros(6)
    for lag in range(min(step+1, len(weights))):
        for i in range(6):
            for j in range(6):
                result[i] += weights[lag, i, j]*history[step-lag, j]
    return result
