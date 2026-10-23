"""CPU reference HOS Dirichlet-to-velocity map on a periodic deep-water graph.

This is an elliptic component for time-domain evolution, not a frequency-domain
wave/body solver. No hull, radiation boundary, or near/far coupling is included.
Arrays use (y, x); lengths and elevations are metres, potential is m²/s.
"""
from dataclasses import dataclass
from math import factorial

import numpy as np
from scipy.signal import resample


def _resize(values, shape):
    # Fourier resampling handles the even-grid Nyquist split explicitly.
    return resample(resample(values, shape[0], axis=0), shape[1], axis=1)


@dataclass(frozen=True)
class GraphField:
    velocity: np.ndarray
    normal_flux: np.ndarray
    flat_terms: tuple[np.ndarray, ...]


def graph_field(eta, psi, lengths, *, order=3, padding=2.):
    """Return Cartesian velocity and order-consistent upward graph flux.

    HOS counts eta and psi as first order: f_1=psi and
    f_m=−sum(eta**j/j! * |D|**j f_(m−j), j=1..m−1).
    Each Cartesian derivative of the harmonic extension is Taylor-expanded at
    z=eta through the same total order. Products are formed on the padded grid;
    increasing padding is a required aliasing control, not an accuracy promise.
    flat_terms retain each f_m on that grid for independent continuation checks.
    The horizontal cross product in the flux retains orders through M−1:
    multiplying its order-M term by the slope would add an unmatched order
    M+1 term. The returned Cartesian velocity still retains all M orders.
    The zero Fourier mode has zero velocity (potential gauge only).
    """
    eta, psi = np.asarray(eta, dtype=float), np.asarray(psi, dtype=float)
    lengths = np.asarray(lengths, dtype=float)
    if (eta.ndim != 2 or min(eta.shape) < 3 or psi.shape != eta.shape
            or not np.isfinite(eta).all() or not np.isfinite(psi).all()):
        raise ValueError("eta and psi must be finite matching two-dimensional grids")
    if lengths.shape != (2,) or not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError("lengths must contain positive finite (Lx, Ly)")
    if isinstance(order, bool) or not isinstance(order, int) or not 1 <= order <= 8:
        raise ValueError("order must be an integer from 1 to 8")
    if not np.isfinite(padding) or padding < 1:
        raise ValueError("padding must be finite and at least one")
    shape = tuple(int(np.ceil(padding*n)) for n in eta.shape)
    height, trace = _resize(eta, shape), _resize(psi, shape)
    kx = 2*np.pi*np.fft.fftfreq(shape[1], d=lengths[0]/shape[1])[None, :]
    ky = 2*np.pi*np.fft.fftfreq(shape[0], d=lengths[1]/shape[0])[:, None]
    wave = np.hypot(kx, ky)

    def derivative(spectrum, power, horizontal=1):
        return np.fft.ifft2(spectrum*wave**power*horizontal).real

    powers = [height**j/factorial(j) for j in range(order)]
    terms, spectra = [trace], [np.fft.fft2(trace)]
    for m in range(2, order+1):
        term = -sum(powers[j]*derivative(spectra[m-j-1], j)
                    for j in range(1, m))
        terms.append(term)
        spectra.append(np.fft.fft2(term))
    velocity = np.zeros((*shape, 3))
    highest_horizontal = np.zeros((*shape, 2))
    for m, spectrum in enumerate(spectra, start=1):
        for j in range(order-m+1):
            horizontal_x = powers[j]*derivative(spectrum, j, 1j*kx)
            horizontal_y = powers[j]*derivative(spectrum, j, 1j*ky)
            velocity[..., 0] += horizontal_x
            velocity[..., 1] += horizontal_y
            if m+j == order:
                highest_horizontal[..., 0] += horizontal_x
                highest_horizontal[..., 1] += horizontal_y
            velocity[..., 2] += powers[j]*derivative(spectrum, j+1)
    surface = np.fft.fft2(height)
    flux = (velocity[..., 2]
            -derivative(surface, 0, 1j*kx)*(velocity[..., 0]-highest_horizontal[..., 0])
            -derivative(surface, 0, 1j*ky)*(velocity[..., 1]-highest_horizontal[..., 1]))
    return GraphField(
        np.stack([_resize(velocity[..., j], eta.shape) for j in range(3)], axis=-1),
        _resize(flux, eta.shape), tuple(terms),
    )
