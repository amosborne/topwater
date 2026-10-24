"""Instantaneous fluid inertia at a flat zero-acceleration-potential surface."""

from dataclasses import dataclass
from time import perf_counter

import numpy as np

from ._validation import array
from .bem import assemble_dirichlet_halfspace
from .boundary import BoundarySystem


@dataclass(frozen=True, eq=False)
class AddedMassResult:
    """World-axis six-motion inertia about reference_point, without symmetrizing.

    Columns/rows are surge, sway, heave, roll, pitch, yaw. Translational entries
    have units kg, mixed entries kg*m, rotational entries kg*m^2. The auxiliary
    acceleration_potential has one column per unit generalized acceleration.
    Fluid acceleration reaction is -added_mass @ acceleration.
    """

    added_mass: np.ndarray
    acceleration_potential: np.ndarray
    reference_point: np.ndarray
    diagnostics: dict


def instantaneous_added_mass(panels, reference_point, *, water_level=0., density=1000., controls=None):
    """Solve six auxiliary problems on cached physical wetted panels at once.

    Supply fluid-outward normals: TrianglePanels(wetted_triangles, -1) when
    triangles retain body-outward winding. reference_point is the world load
    reference, normally the CG.
    The infinite-depth zero-potential plane models calm startup and the
    infinite-frequency linear limit. It does not model finite-frequency waves or an evolving
    surface. This inertia is part of pressure/acceleration coupling, never an
    additional load placed on top of the same total-pressure integral.
    """
    start = perf_counter()
    reference = array(reference_point, (3,), 'reference_point')
    rho = array(density, (), 'density').item()
    if rho <= 0:
        raise ValueError('density must be positive')
    operators = assemble_dirichlet_halfspace(panels, water_level=water_level, controls=controls)
    length = float(np.ptp(panels.triangles.reshape(-1, 3), axis=0).max())
    system = BoundarySystem(operators, dirichlet=np.zeros(len(panels), bool),
        areas=panels.areas, length_scale=length, bounded=False)
    modes = np.column_stack((panels.normals,
                             np.cross(panels.centroids-reference, panels.normals)))
    solution = system.solve(normal_velocity=modes)
    mass = rho*(modes.T@(panels.areas[:, None]*solution.potential))
    mass.setflags(write=False)
    scale = np.array([1., 1., 1., length, length, length])
    normalized = mass/scale[:, None]/scale[None, :]/(rho*length**3)
    diagnostics = dict(assembly=operators.diagnostics, factorization=system.diagnostics,
        solve=solution.diagnostics, density=rho, normalization_length=length,
        normalized_reciprocity_defect=float(np.max(np.abs(normalized-normalized.T))),
        normalized_energy_eigenvalues=np.linalg.eigvalsh((normalized+normalized.T)/2).tolist())
    diagnostics['total_seconds'] = perf_counter()-start
    return AddedMassResult(mass, solution.potential, reference, diagnostics)
