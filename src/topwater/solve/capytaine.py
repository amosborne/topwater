"""Version-specific Capytaine potential data shared by numeric backends."""
from importlib.metadata import version
from types import SimpleNamespace
import numpy as np


def potential_inputs(green_function, source_mesh):
    """Read version-specific Capytaine objects at this adapter boundary only."""
    if (version('capytaine') != '3.0.0'
            or green_function.gf_singularities != 'low_freq'
            or green_function.tabulation_grid_shape != 'scaled_nemoh3'
            or green_function.floating_point_precision != 'float64'):
        raise ValueError('Potential fields require Capytaine 3.0.0 float64 low_freq scaled_nemoh3 data')
    mesh = source_mesh.merged()
    quadrature, weights = mesh.quadrature_points
    if np.any(quadrature[:, :, 2] > 0):
        raise ValueError('Potential fields require submerged quadrature points')
    return SimpleNamespace(pi=np.float64(green_function.fortran_core.constants.pi),
        gamma=np.float64(green_function.fortran_core.constants.euler_gamma),
        nodes=mesh.vertices[mesh.faces], centers=mesh.faces_centers,
        normals=mesh.faces_normals, areas=mesh.faces_areas, radii=mesh.faces_radiuses,
        quadrature=quadrature, weights=weights,
        rr=green_function.tabulated_r_range, zz=green_function.tabulated_z_range,
        table=green_function.tabulated_integrals[:, :, [0, 2]])
