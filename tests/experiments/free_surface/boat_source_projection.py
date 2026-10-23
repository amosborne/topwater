"""Project oversampled analytic source traces onto the owned odd solver grid."""
import numpy as np
from .cuda_graph import cp


def project_trace(values, sample_size, size):
    """Retain Fourier modes; values have one row per sample, optional columns."""
    values = cp.asarray(values)
    if (values.ndim not in (1, 2) or values.shape[0] != sample_size**2
            or size < 3 or size % 2 != 1 or sample_size < size or sample_size % 2 != 1):
        raise ValueError('odd source/solver grids and one source value row per sample required')
    columns = values.shape[1:]
    spectrum = cp.fft.fft2(values.reshape(sample_size, sample_size, *columns), axes=(0, 1))
    indices = cp.asarray(np.r_[np.arange((size+1)//2),
                               np.arange(sample_size-(size-1)//2, sample_size)])
    retained = spectrum[indices[:, None], indices[None, :]]
    projected = cp.fft.ifft2(retained, axes=(0, 1)).real*(size/sample_size)**2
    return projected.reshape(size**2, *columns)
