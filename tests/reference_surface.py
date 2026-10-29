"""Historical causal surface operators retained only as numerical test controls."""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import trapezoid
from scipy.fft import rfft, irfft, next_fast_len
from topwater.solve.validation import _finite, _grid


@dataclass
class SurfaceKernel:
    lag_s: np.ndarray
    values: np.ndarray
    direct_velocity: np.ndarray
    accepted: bool = False
    direct_acceleration: np.ndarray | None = None

def cosine_kernel(omega, velocity_transfer, lag_s, *, direct_velocity, direct_acceleration=None):
    """Candidate regular causal kernel, never an automatic acceptance gate.

    An explicitly identified instantaneous velocity term is subtracted before
    cosine inversion.  Other pole/derivative terms must be removed analytically
    by the caller; this function must not be used to hide them with finite-band
    integration.  Supplied frequency endpoints/tails remain caller-reviewed.
    The returned candidate cannot publish a total field until accepted.
    """
    omega = _grid(omega, 'omega', zero=True)
    transfer = _finite(velocity_transfer, (len(omega), None, 6), 'velocity_transfer', complex_value=True)
    direct = _finite(direct_velocity, transfer.shape[1:], 'direct_velocity')
    acceleration = (np.zeros_like(direct) if direct_acceleration is None else
                    _finite(direct_acceleration, transfer.shape[1:], 'direct_acceleration'))
    lag = _grid(lag_s, 'lag_s', zero=True)
    if lag[0] != 0:
        raise ValueError('lag_s must start at zero')
    regular = (transfer-direct[None]+1j*omega[:, None, None]*acceleration[None]).real
    delta = np.diff(omega)
    weights = np.r_[delta[0]/2, (delta[:-1]+delta[1:])/2, delta[-1]/2]
    weighted = (weights[:, None, None]*regular).reshape(len(omega), -1)
    values = np.empty((len(lag), *transfer.shape[1:]))
    for start in range(0, len(lag), 64):
        cosine = np.cos(lag[start:start+64, None]*omega)
        values[start:start+64] = (2/np.pi*(cosine @ weighted)).reshape(-1, *transfer.shape[1:])
    return SurfaceKernel(lag, values, direct, direct_acceleration=acceleration)

def recover_transfer(kernel, omega):
    """Full complex response of a finite candidate kernel at held-out omega."""
    omega = _finite(omega, (None,), 'omega')
    output = np.empty((len(omega), *kernel.values.shape[1:]), dtype=complex)
    for i, w in enumerate(omega):
        output[i] = kernel.direct_velocity + trapezoid(
            kernel.values*np.exp(1j*w*kernel.lag_s)[:, None, None], x=kernel.lag_s, axis=0)
        if kernel.direct_acceleration is not None:
            output[i] -= 1j*w*kernel.direct_acceleration
    return output

def convolve_radiation(kernel, time_s, velocity, acceleration=None):
    """Zero-prehistory radiation from actual fixed-mean-CG modal velocities.

    These are the six native radiation-input velocities, not Euler angle rates.
    Lag coverage must extend through the whole requested trajectory, regardless
    of the different memory cutoff used by WEC-Sim's body-force operator.
    """
    if not kernel.accepted:
        raise ValueError('Transient surface kernel has not passed endpoint/tail and full-complex-transfer review')
    times = _grid(time_s, 'time_s', zero=True)
    v = _finite(velocity, (len(times), 6), 'velocity')
    dt = float(times[1]-times[0])
    if times[0] != 0 or not np.allclose(np.diff(times), dt, rtol=0, atol=1e-12):
        raise ValueError('Radiation history requires a uniform grid starting at zero')
    if np.any(v[0] != 0):
        raise ValueError('This reconstruction requires exact zero initial modal velocity')
    if len(kernel.lag_s) < len(times) or not np.allclose(kernel.lag_s[:len(times)], times, rtol=0, atol=1e-12):
        raise ValueError('Surface kernel must cover the entire requested history at its saved cadence')
    result = v @ kernel.direct_velocity.T
    if kernel.direct_acceleration is not None and np.any(kernel.direct_acceleration != 0):
        if acceleration is None:
            raise ValueError('Actual modal accelerations are required for the identified instantaneous acceleration term')
        a = _finite(acceleration, v.shape, 'acceleration')
        result += a @ kernel.direct_acceleration.T
    h = kernel.values[:len(times)]
    n = next_fast_len(2*len(times)-1)
    product = np.sum(rfft(h, n=n, axis=0)*rfft(v, n=n, axis=0)[:, None, :], axis=2)
    full = irfft(product, n=n, axis=0)[:len(times)]
    # Endpoint half weights implement integral_0^t at each t. Both endpoints
    # coincide at t=0, where the integral is exactly zero.
    result += dt*(full-.5*(v @ h[0].T)-.5*np.sum(h*v[0][None, None, :], axis=2))
    return result


from topwater.solve.types import MODES
from topwater.solve.progress import _StageProgress, _summarize_capytaine_warnings
from topwater.solve.surface import exterior_field_points, _single_threaded

@_summarize_capytaine_warnings
@_single_threaded
def radiation_velocity_transfer(hydro, xy, omega=None, *, checkpoint=None, shared_green=False, allow_boundary=False):
    """Return H_eta,v=i H_eta,q/omega, shape (omega, point, six modes)."""
    import capytaine as cpt
    checkpoint = (lambda: None) if checkpoint is None else checkpoint
    xy = exterior_field_points(hydro, xy, allow_boundary=allow_boundary)
    frequencies = hydro.radiation_omega if omega is None else _grid(omega, 'omega')
    values = np.empty((len(frequencies), len(xy), 6), dtype=complex)
    if shared_green:
        # DefaultMatrixEngine.build_S_matrix merges (copies) this immutable
        # source mesh on every call, discarding all cached geometric properties.
        # Use its exact public Green-function evaluator on one merged mesh.
        # No quadrature, field points or source coefficients change.
        source_mesh = hydro.body.mesh_including_lid.merged()
        field_points = np.c_[xy, np.zeros(len(xy))]
    progress = _StageProgress('radiation transfers', len(frequencies), 'frequencies')
    for iw, w in enumerate(frequencies):
        results = []
        for j, mode in enumerate(MODES):
            checkpoint()
            key = (float(w), mode)
            if key not in hydro.radiation_results:
                hydro.radiation_results[key] = hydro.solver.solve(cpt.RadiationProblem(
                    body=hydro.body, rho=hydro.density, g=hydro.gravity,
                    water_depth=hydro.depth, omega=float(w), radiating_dof=mode), keep_details=True)
            results.append(hydro.radiation_results[key])
            if not shared_green:
                eta_q = hydro.solver.compute_free_surface_elevation(xy, hydro.radiation_results[key])
                values[iw, :, j] = 1j*eta_q/float(w)
        if shared_green:
            # This is precisely the public field evaluator's S @ sources,
            # batching the six source columns. No new BEM equation or solver.
            representative = results[0]
            S, _ = hydro.solver.engine.green_function.evaluate(field_points,
                source_mesh, free_surface=representative.free_surface,
                water_depth=representative.water_depth, wavenumber=representative.encounter_wavenumber)
            if not np.isfinite(S).all():
                raise ValueError('Nonfinite Capytaine single-layer field matrix')
            sources = np.column_stack([r.sources for r in results])
            # i/omega * (i*omega/g) * potential = -potential/g.
            values[iw] = -(S @ sources)/hydro.gravity
        progress.update(iw+1)
    progress.finish()
    return values
