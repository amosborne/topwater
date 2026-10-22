"""Fixed-reference surface test space; no guarantee of current physical bandwidth.

Cosines satisfy the reference rectangle's sidewall normal condition. The hull
hole and the ALE map prevent treating them as exact physical tank eigenmodes.
"""
import numpy as np


def rectangle_cosines(reference_xy, beam, *, enrichment=1.):
    """Dimensionless cosine columns and integer (m,n), ordered by radial k.

    The reference cutoff is enrichment*pi/beam. Nonconstant factors sqrt(2)
    balance rectangle column norms; the perforated-domain mass is not diagonal.
    R is fixed on master labels throughout the physical trajectory.
    """
    xy = np.asarray(reference_xy, dtype=float)
    if (xy.ndim != 2 or xy.shape[1] != 2 or not np.all(np.isfinite(xy))
            or not np.isfinite(beam) or beam <= 0
            or not np.isfinite(enrichment) or enrichment < 1):
        raise ValueError('finite reference XY, positive beam and enrichment at least one required')
    origin = xy.min(axis=0)
    length = xy.max(axis=0)-origin
    if np.any(length <= 0):
        raise ValueError('reference rectangle must have positive extent in both axes')
    cutoff = enrichment*np.pi/beam
    maximum = np.ceil(cutoff*length/np.pi).astype(int)
    modes = [(m, n) for m in range(maximum[0]+1) for n in range(maximum[1]+1)
             if (m*np.pi/length[0])**2+(n*np.pi/length[1])**2 <= cutoff**2*(1+32*np.finfo(float).eps)]
    modes.sort(key=lambda pair: ((pair[0]/length[0])**2+(pair[1]/length[1])**2, *pair))
    modes = np.asarray(modes, dtype=int)
    basis = np.prod(np.cos((xy[:, None]-origin)*np.pi*modes[None]/length)
                    *np.where(modes[None] == 0, 1., np.sqrt(2.)), axis=2)
    return basis, modes
