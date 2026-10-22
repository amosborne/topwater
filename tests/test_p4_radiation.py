"""Physical P4 surface integration controls, no tank factorization required."""
import numpy as np

from .experiments.free_surface import p4_tetra
from .experiments.free_surface.canoe_radiation import sponge
from .experiments.free_surface.p4_radiation import physical_sponge, surface_matrices


def test_physical_sponge_matches_existing_ramp_without_nodal_interpolation():
    points = np.array([[2*np.pi, np.pi, 0.], [0., 0., 0.],
                       [2.7, 1.8, 0.], [-4.1, .8, 0.]])
    actual, widths = physical_sponge(points, np.array([2*np.pi, np.pi]))
    expected, expected_widths = sponge(points, 4., .5)
    assert np.array_equal(actual, expected)
    assert np.array_equal(widths, expected_widths)
    assert actual[1] == 0


def test_complete_p4_surface_mass_and_positive_damping():
    a, b = 2*np.pi, np.pi
    points = np.array([[-a, -b, 0.], [a, -b, 0.], [a, b, 0.],
                       [-a, b, 0.], [0., 0., -4.]])
    mesh = p4_tetra.make_mesh(points, np.array([[0, 1, 2, 4], [0, 2, 3, 4]]))
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    free, mass, weighted, _ = surface_matrices(mesh, faces, 12)
    assert len(free) == 25
    area = 4*a*b
    assert abs(mass.sum()-area) < 1e-11
    # Independent exact rectangle integral of (x²)² verifies complete mass,
    # not a lumped/diagonal approximation at the interpolation nodes.
    polynomial = mesh.points[free, 0]**2
    assert abs(polynomial@(mass@polynomial)/(area*a**4/5)-1) < 2e-14
    assert np.linalg.eigvalsh(mass.toarray()).min() > 0
    assert np.linalg.eigvalsh(weighted.toarray()).min() > -1e-12
    assert np.linalg.eigvalsh((4*mass-weighted).toarray()).min() > -1e-12
