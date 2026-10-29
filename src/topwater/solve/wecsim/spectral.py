"""Physical FullDir waves and ordinary coefficient inputs for WEC-Sim.

This module prepares waves and coefficients separately. It never multiplies
wave coefficients by hydrodynamic coefficients or constructs applied forces.
"""
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator


@dataclass(frozen=True)
class SpectralInputs:
    wave: dict
    hydro: dict
    coefficients: np.ndarray
    original_omega_rad_s: np.ndarray
    native_omega_rad_s: np.ndarray
    omega_roundtrip_error_rad_s: np.ndarray
    directions_rad: np.ndarray
    domega_rad_s: np.ndarray
    dtheta_rad: np.ndarray


def _finite(value, name, *, real=False):
    result = np.asarray(value)
    if (not np.issubdtype(result.dtype, np.number) or not result.size
            or not np.isfinite(result).all() or (real and not np.isrealobj(result))):
        raise ValueError(f'{name} must contain finite numeric values')
    return result


def _grid(value, name, *, first=None):
    result = _finite(value, name, real=True)
    if result.ndim != 1 or len(result) < 2 or np.any(np.diff(result) <= 0):
        raise ValueError(f'{name} must be a strictly increasing vector with at least two nodes')
    if first is not None and result[0] != first:
        raise ValueError(f'{name} must start at {first}')
    return result.astype(float, copy=False)


def _native_frequency_widths(omega):
    # waveClass ImportedFullDir endpoint rule, not the packet's original
    # integration weights (which are ALREADY part of its complex C).
    return np.r_[omega[1]-omega[0], (omega[2:]-omega[:-2])/2, omega[-1]-omega[-2]]


def _native_direction_widths(degrees):
    gaps = np.r_[degrees[1]-degrees[0], (degrees[2:]-degrees[:-2])/2,
                 degrees[-1]-degrees[-2]]
    return np.deg2rad(gaps)


def _coherent_waves(train):
    original = _grid(train.omega, 'packet frequencies', first=0.)
    theta = np.asarray(train.theta)
    amplitudes = np.asarray(train.carriers.amplitudes)
    C = np.asarray(train.C)
    if (any(not np.issubdtype(value.dtype, np.number) for value in (theta, amplitudes, C))
            or not np.isrealobj(theta) or not np.isrealobj(amplitudes)):
        raise ValueError('Packet directions/amplitudes must be real numeric and coefficients numeric')
    if theta.size or amplitudes.size or C.size:
        theta = _finite(theta, 'physical carrier directions', real=True)
        amplitudes = _finite(amplitudes, 'carrier amplitudes', real=True)
        C = _finite(C, 'physical packet coefficients')
    if (theta.ndim != 1 or amplitudes.shape != theta.shape or np.any(amplitudes < 0)
            or C.shape != (len(theta), 2, len(original))):
        raise ValueError('Packet coefficients must have carrier, signed branch, frequency axes')
    if np.any(C[:, :, 0] != 0):
        raise ValueError('Only an exactly zero physical omega=0 row may be omitted')
    if np.any(C[amplitudes == 0] != 0):
        raise ValueError('Inactive physical carriers contain nonzero coefficients')
    active = np.flatnonzero(amplitudes)
    if len(active):
        directions = np.mod(theta[active, None]+np.array([0., np.pi]), 2*np.pi).reshape(-1)
        physical_rows = C[active, :, 1:].reshape(-1, len(original)-1)
        nonzero = np.any(physical_rows != 0, axis=1)
        directions, physical_rows = directions[nonzero], physical_rows[nonzero]
        unique, inverse = np.unique(directions, return_inverse=True)
        rows = np.zeros((len(unique), len(original)-1), dtype=complex)
        np.add.at(rows, inverse, physical_rows)
        if not len(unique):
            unique = np.array([0., np.pi])
            rows = np.zeros((2, len(original)-1), complex)
        elif len(unique) == 1:
            # FullDir's directional quadrature needs two headings. The added
            # column is exactly zero and creates no counter-propagating wave.
            unique = np.r_[unique, np.mod(unique[0]+np.pi, 2*np.pi)]
            rows = np.vstack([rows, np.zeros((1, len(original)-1), complex)])
    else:
        # A calm sea needs a valid two-column FullDir representation. Both
        # columns are exactly zero, so their headings add no physical waves.
        unique = np.array([0., np.pi])
        rows = np.zeros((2, len(original)-1), dtype=complex)
    wrapped = np.where(unique >= np.pi, unique-2*np.pi, unique)
    degrees = np.rad2deg(wrapped)
    order = np.argsort(degrees)
    degrees, unique, rows = degrees[order], unique[order], rows[order]
    if len(degrees) < 2 or np.any(np.diff(degrees) <= 0):
        raise ValueError('Distinct physical headings collapse in the native degree representation')
    positive = original[1:]
    if len(positive) < 2 or positive[0] <= 0:
        raise ValueError('FullDir requires at least two positive physical packet frequencies')
    frequency_hz = positive/(2*np.pi)
    # Match waveClass's EXACT left-to-right imported-frequency operation.
    native = (frequency_hz*2)*np.pi
    if np.any(np.diff(native) <= 0) or not np.isfinite(native).all():
        raise ValueError('Native Hz frequency conversion loses distinct packet nodes')
    return positive, native, frequency_hz, unique, degrees, rows.T


def _spectral_hydro(hydro, native_omega, directions_rad, directions_deg):
    source_w = _grid(hydro.excitation_omega, 'native coefficient frequencies', first=.1)
    source_theta = _grid(hydro.excitation_headings, 'native coefficient headings', first=0.)
    H = _finite(hydro.excitation, 'native dimensional excitation')
    K = _finite(hydro.K, 'physical stiffness', real=True)
    density = _finite(hydro.density, 'density', real=True)
    gravity = _finite(hydro.gravity, 'gravity', real=True)
    if (source_theta[-1] != 2*np.pi or H.shape != (len(source_w), len(source_theta), 6)
            or K.shape != (6, 6) or density.size != 1 or gravity.size != 1
            or density.item() <= 0 or gravity.item() <= 0):
        raise ValueError('Incomplete six-mode periodic excitation operator or physical constants')
    if not np.array_equal(H[:, 0], H[:, -1]):
        raise ValueError('Native excitation must retain its exact periodic seam')
    if native_omega[-1] > source_w[-1]:
        raise ValueError('The complete physical packet exceeds genuine coefficient support')
    # Include a genuine upper coverage boundary, not an extra wave/quadrature
    # node. This prevents a second Hz comparison from dropping the last packet row.
    table_w = np.unique(np.r_[0., native_omega, source_w[-1]])
    table_theta = np.unique(np.r_[-180., directions_deg, 180.])
    by_degree = dict(zip(directions_deg, directions_rad))
    query_theta = np.array([np.pi if abs(value) == 180. else by_degree[value]
                            for value in table_theta])
    query = np.empty((len(table_w), len(table_theta), 2))
    query[:, :, 0] = np.maximum(table_w, source_w[0])[:, None]
    query[:, :, 1] = query_theta[None, :]
    sampled = RegularGridInterpolator((source_w, source_theta), H, bounds_error=True)(query)
    low = table_w < source_w[0]
    limit = K[:, 2]
    sampled[low] = limit+(table_w[low, None, None]/source_w[0])**2*(sampled[low]-limit)
    factor = float(density.item()*gravity.item())
    real = sampled.real.transpose(2, 1, 0).copy()/factor
    imaginary = -sampled.imag.transpose(2, 1, 0).copy()/factor
    if not np.isfinite(real).all() or not np.isfinite(imaginary).all():
        raise ValueError('Derived excitation coefficients are nonfinite')
    # These are existing public bodyClass(hydroStruct) fields, not a WEC extension.
    return dict(simulation_parameters=dict(w=table_w[:, None], direction=table_theta[:, None]),
                hydro_coeffs=dict(excitation=dict(re=real, im=imaginary),
                                  mean_drift=np.zeros_like(real)))


def prepare_spectral_inputs(setup, hydro):
    """Prepare exact coherent waves and the old coefficient operator's ordinary input table."""
    original, native, frequency_hz, directions, degrees, C = _coherent_waves(setup.train)
    dw, dtheta = _native_frequency_widths(native), _native_direction_widths(degrees)
    if np.any(dw <= 0) or np.any(dtheta <= 0):
        raise ValueError('Native spectrum weights must be positive without floors')
    squared = abs(C)**2
    if np.any((C != 0) & (squared == 0)):
        raise ValueError('Nonzero physical coefficients underflow in the spectrum representation')
    energy = np.sum(squared, axis=1)
    spectrum = np.pi*energy/dw
    spread = np.zeros_like(squared)
    np.divide(squared, energy[:, None]*dtheta[None, :], out=spread, where=energy[:, None] > 0)
    if np.any((squared > 0) & (spread == 0)):
        raise ValueError('Nonzero physical coefficients underflow in directional spectrum weights')
    phase = np.where(squared > 0, -np.angle(C), 0.)
    if not all(np.isfinite(value).all() for value in (C, spectrum, spread, phase)):
        raise ValueError('Physical spectrum representation overflowed or became nonfinite')
    wave = dict(frequencies=frequency_hz[:, None], directions=degrees[:, None],
                spectrum=spectrum[:, None], spread=spread, phase=phase)
    standard = _spectral_hydro(hydro, native, directions, degrees)
    return SpectralInputs(wave, standard, C.copy(), original.copy(), native.copy(),
                          native-original, directions.copy(), dw, dtheta)
