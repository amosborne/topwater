"""Independent geometry integrals used by hydrostatics regressions."""
import numpy as np


def surface_integrals(vertices, faces):
    """Independent divergence-theorem reference integrals of planar triangles.

    Uses z*n_z for volume (the input validator uses tetrahedral determinants).
    First moments use exact triangle polynomial averages.
    """
    triangles = vertices[faces]
    p, q, r = triangles.transpose(1, 0, 2)
    normal_area = np.cross(q - p, r - p) / 2
    volume = np.sum(normal_area[:, 2] * triangles[:, :, 2].mean(axis=1))
    average_squared = (p*p + q*q + r*r + p*q + p*r + q*r) / 6
    first_moments = np.sum(normal_area * average_squared, axis=0) / 2
    return volume, first_moments / volume
