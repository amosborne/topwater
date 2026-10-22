"""Array-only Cummins TIME reference; never a Topwater fluid/time-boundary model.

F_rad(t)=-A_inf*a(t)-integral_0^t K(t-tau)*V(tau)dtau,
K(t)=(2/pi)*integral_0^infinity B(omega)*cos(omega*t)domega.
See https://wec-sim.github.io/WEC-Sim/main/theory/theory.html#convolution-integral-formulation
and https://capytaine.org/stable/user_manual/problem_setup.html .

Only the SUPPLIED frequency band is integrated. No missing low/high tail,
window, fitted mode, mean removal, matrix symmetrization or PSD projection is
inserted. Three existing omega2/4/6 samples do not establish a time reference.
The force history needs independent band/spacing/time-grid/memory and Capytaine
panel controls. A_inf is REQUIRED from the SAME depth4 hull, CG, density and
six-mode ordering; the older infinite-depth helper is not that input. Export
provenance must bind those quantities and all36 raw A/B/A_inf entries before
comparison. This module does not run Capytaine or claim those controls passed.
"""
import numpy as np
from scipy.integrate import trapezoid


def _real(values, shape, name):
    if np.iscomplexobj(values):
        raise ValueError(name+' requires real physical coefficients/history')
    values = np.asarray(values, dtype=float)
    if values.shape != shape or not np.all(np.isfinite(values)):
        raise ValueError(name+' has incompatible shape or nonfinite values')
    return values


def _times(times):
    if np.iscomplexobj(times):
        raise ValueError('real physical reference times required')
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0
            or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0)):
        raise ValueError('finite increasing reference times starting at zero required')
    return times


def radiation_kernel(omega, damping, times):
    """Finite-band cosine quadrature of all36 dimensional damping entries.

    Omega is in rad/s, not Hz. A positive first sample leaves [0,omega_min)
    omitted explicitly; a supplied omega0 value is never inferred from a
    finite-depth Capytaine solve (that exact endpoint is unsupported there).
    Times can be nonuniform. Returns K[time,force_row,velocity_column].
    """
    if np.iscomplexobj(omega):
        raise ValueError('real angular-frequency samples required')
    omega = np.asarray(omega, dtype=float)
    times = _times(times)
    if (omega.ndim != 1 or len(omega) < 2 or not np.all(np.isfinite(omega))
            or omega[0] < 0 or np.any(np.diff(omega) <= 0)):
        raise ValueError('finite increasing nonnegative angular-frequency samples required')
    damping = _real(damping, (len(omega), 6, 6), 'radiation damping')
    result = np.empty((len(times), 6, 6))
    # Small time batches keep the reference quadrature tensor bounded. This is
    # direct integration of supplied data, not a fitted radiation state model.
    for start in range(0, len(times), 16):
        cosine = np.cos(times[start:start+16, None]*omega)
        result[start:start+16] = 2/np.pi*trapezoid(
            cosine[:, :, None, None]*damping[None], omega, axis=1)
    return result


def radiation_force(times, velocity, acceleration, added_mass_infinite, kernel, *,
                    added_mass_infinite_depth):
    """Zero-memory causal full-coupling force on a uniform reference time grid.

    The required depth assertion prevents silently using the existing deep-
    water A_inf helper for this depth4 case; it does not replace source/mesh/CG
    provenance. Velocity/acceleration are WORLD axes about the SAME fixed CG.
    Supplied history starts at t0 with zero prehistory radiation memory; it is
    not a steady periodic convolution. No hydrostatic or incident-wave load is
    included: compare radiation alone or add independently known restoring.
    """
    if added_mass_infinite_depth != 4.:
        raise ValueError('this canoe reference requires explicit depth4 A_inf data')
    times = _times(times)
    step = times[1]-times[0]
    if np.max(abs(np.diff(times)-step)) > 64*np.finfo(float).eps*max(1., times[-1]):
        raise ValueError('uniform reference times required; refine/interpolate explicitly outside this adapter')
    velocity = _real(velocity, (len(times), 6), 'world velocity')
    acceleration = _real(acceleration, (len(times), 6), 'world acceleration')
    added = _real(added_mass_infinite, (6, 6), 'infinite-frequency added mass')
    kernel = _real(kernel, (len(times), 6, 6), 'causal radiation kernel')
    memory = np.zeros((len(times), 6))
    for row in range(6):
        for column in range(6):
            memory[:, row] += np.convolve(kernel[:, row, column], velocity[:, column])[:len(times)]
    # Composite causal trapezoid: both moving endpoints have half weight.
    # At t0 these corrections cancel the sole convolution term exactly.
    memory -= .5*np.einsum('tij,j->ti', kernel, velocity[0])
    memory -= .5*(velocity@kernel[0].T)
    memory[0] = 0.  # Exactly zero integration interval, including nonzero initial V.
    return -acceleration@added.T-step*memory


def reconstruct_coefficients(times, kernel, added_mass_infinite, omega):
    """Finite-memory A/B back-check, not a fit or a Topwater frequency solve.

    For exp(-i*omega*t), integral K(t)*exp(i*omega*t)dt
    = B(omega)+i*omega*(A_inf-A(omega)). Nonzero residuals expose truncation,
    grid, A_inf or reference inconsistencies; no correction is applied.
    """
    times = _times(times)
    kernel = _real(kernel, (len(times), 6, 6), 'radiation kernel')
    added = _real(added_mass_infinite, (6, 6), 'infinite-frequency added mass')
    if np.iscomplexobj(omega):
        raise ValueError('real angular frequencies required for A/B back-check')
    omega = np.asarray(omega, dtype=float)
    if omega.ndim != 1 or not len(omega) or not np.all(np.isfinite(omega)) or np.any(omega <= 0):
        raise ValueError('positive finite angular frequencies required for A/B back-check')
    response = np.array([trapezoid(np.exp(1j*w*times)[:, None, None]*kernel, times, axis=0)
                         for w in omega])
    return added-response.imag/omega[:, None, None], response.real
