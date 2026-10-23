"""Fixed-graph GPU HOS Dirichlet map; odd grids, float64/complex128 only.

Coefficients are padded FFTs of the harmonic flat-plane terms. Their sum can
be continued to submerged points by the caller; this module owns no body solve.
"""
from math import factorial
import numpy as np
from .cuda_operator import cp


_body_maps = cp.ElementwiseKernel(
    'raw float64 points, raw float64 normals, raw float64 origin, '
    'raw float64 kx, raw float64 ky, raw float64 wave, int64 modes',
    'complex128 potential, complex128 normal',
    '''
    const long long target = i / modes;
    const long long mode = i % modes;
    const double phase = (points[3*target]-origin[0])*kx[mode]
                       + (points[3*target+1]-origin[1])*ky[mode];
    const double amplitude = exp(points[3*target+2]*wave[mode])/modes;
    potential = complex<double>(amplitude*cos(phase), amplitude*sin(phase));
    normal = potential * complex<double>(normals[3*target+2]*wave[mode],
        normals[3*target]*kx[mode] + normals[3*target+1]*ky[mode]);
    ''', 'graph_body_harmonic_maps')


def body_maps(points, normals, size, spacing, origin):
    """Build flat harmonic potential/normal maps directly in GPU storage.

    Points use z relative to the expansion plane; normals are fluid outward.
    The independent NumPy formula remains in the radiation verification harness.
    """
    points, normals, origin = (np.ascontiguousarray(a, dtype=np.float64)
                               for a in (points, normals, origin))
    if (points.ndim != 2 or points.shape[1:] != (3,) or normals.shape != points.shape
            or origin.shape != (2,) or not all(np.isfinite(a).all() for a in (points, normals, origin))
            or np.any(points[:, 2] > 0)):
        raise ValueError('finite submerged points, matching normals and two-dimensional origin required')
    if isinstance(size, bool) or not isinstance(size, int) or size < 3 or size % 2 == 0 or not np.isfinite(spacing) or spacing <= 0:
        raise ValueError('odd integer size>=3 and positive finite spacing required')
    kx, ky = cp.meshgrid(2*np.pi*cp.fft.fftfreq(size, d=spacing),
                         2*np.pi*cp.fft.fftfreq(size, d=spacing))
    wave = cp.hypot(kx, ky)
    potential, normal = _body_maps(cp.asarray(points), cp.asarray(normals), cp.asarray(origin),
        kx.ravel(), ky.ravel(), wave.ravel(), np.int64(size*size), size=len(points)*size*size)
    return potential.reshape(len(points), size*size), normal.reshape(len(points), size*size)


def _resize(values, shape):
    """Fourier interpolate/project between odd grids, preserving amplitudes."""
    old = values.shape
    spectrum = cp.fft.fftshift(cp.fft.fft2(values))
    target = cp.zeros(shape, dtype=cp.complex128)
    width = tuple(min(a, b) for a, b in zip(old, shape, strict=True))
    source_slice = tuple(slice((a-w)//2, (a+w)//2) for a, w in zip(old, width, strict=True))
    target_slice = tuple(slice((a-w)//2, (a+w)//2) for a, w in zip(shape, width, strict=True))
    target[target_slice] = spectrum[source_slice]
    return cp.fft.ifft2(cp.fft.ifftshift(target)).real*(np.prod(shape)/np.prod(old))


class GraphDirichlet:
    """One owned current graph; expansion order counts eta and psi as order one."""
    def __init__(self, eta, lengths, *, order=4, padding=2.):
        eta = cp.array(eta, copy=True)
        lengths = np.asarray(lengths, dtype=float)
        if (eta.ndim != 2 or min(eta.shape) < 3 or any(n % 2 == 0 for n in eta.shape)
                or eta.dtype != cp.float64 or not bool(cp.isfinite(eta).all())):
            raise ValueError('eta must be a finite float64 odd two-dimensional grid')
        if lengths.shape != (2,) or not np.isfinite(lengths).all() or np.any(lengths <= 0):
            raise ValueError('positive finite (Lx, Ly) required')
        if isinstance(order, bool) or not isinstance(order, int) or not 1 <= order <= 8:
            raise ValueError('integer order 1..8 required')
        if not np.isfinite(padding) or padding < 1:
            raise ValueError('finite padding >= 1 required')
        self.shape, self.order = eta.shape, order
        # Smallest odd dimension >= padding*n; padding2 gives 2*n+1.
        self.padded_shape = tuple(2*int(np.ceil((padding*n-1)/2))+1 for n in eta.shape)
        height = _resize(eta, self.padded_shape)
        self._powers = tuple(height**j/factorial(j) for j in range(order))
        self._kx = 2*np.pi*cp.fft.fftfreq(self.padded_shape[1], d=lengths[0]/self.padded_shape[1])[None, :]
        self._ky = 2*np.pi*cp.fft.fftfreq(self.padded_shape[0], d=lengths[1]/self.padded_shape[0])[:, None]
        self._wave = cp.hypot(self._kx, self._ky)
        spectrum = cp.fft.fft2(height)
        self._slope = tuple(self._derivative(spectrum, 0, 1j*k) for k in (self._kx, self._ky))

    def _derivative(self, spectrum, power, horizontal=1):
        return cp.fft.ifft2(spectrum*self._wave**power*horizontal).real

    def coefficients(self, trace):
        """Return padded complex128 FFT(f_m), linear in trace at this graph."""
        trace = cp.asarray(trace)
        if trace.shape != self.shape or trace.dtype != cp.float64 or not bool(cp.isfinite(trace).all()):
            raise ValueError('finite float64 trace matching the current graph required')
        spectra = [cp.fft.fft2(_resize(trace, self.padded_shape))]
        for m in range(2, self.order+1):
            term = -sum(self._powers[j]*self._derivative(spectra[m-j-1], j)
                        for j in range(1, m))
            spectra.append(cp.fft.fft2(term))
        return tuple(spectra)

    def surface(self, spectra):
        """Return (Cartesian velocity, order-consistent upward graph flux).

        Velocity retains M orders. Flux excludes slope times the order-M
        horizontal term: that product would be unmatched order M+1. No mean
        projection is applied. All returned arrays are on the original grid.
        """
        if (len(spectra) != self.order or any(s.shape != self.padded_shape
                or s.dtype != cp.complex128 for s in spectra)):
            raise ValueError('matching padded complex128 coefficient spectra required')
        velocity = cp.zeros((*self.padded_shape, 3), dtype=cp.float64)
        highest = cp.zeros((*self.padded_shape, 2), dtype=cp.float64)
        for m, spectrum in enumerate(spectra, start=1):
            for j in range(self.order-m+1):
                for axis, k in enumerate((self._kx, self._ky)):
                    term = self._powers[j]*self._derivative(spectrum, j, 1j*k)
                    velocity[..., axis] += term
                    if m+j == self.order:
                        highest[..., axis] += term
                velocity[..., 2] += self._powers[j]*self._derivative(spectrum, j+1)
        flux = velocity[..., 2]-sum(slope*(velocity[..., axis]-highest[..., axis])
                                   for axis, slope in enumerate(self._slope))
        return (cp.stack([_resize(velocity[..., j], self.shape) for j in range(3)], axis=-1),
                _resize(flux, self.shape))
