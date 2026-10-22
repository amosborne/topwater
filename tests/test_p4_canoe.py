"""Independent weak-form and rigid-load checks without optional meshing tools."""
import numpy as np

from .experiments.free_surface import p4_tetra
from .experiments.free_surface.p4_canoe import harmonic_quartic, loads
from .experiments.free_surface.pocket import boundary


def test_harmonic_quartic_and_uniform_flow_boundary_identity():
    points = np.array([[0., 0., 0.], [1., 0., 0.],
                       [0., 1., 0.], [0., 0., 1.]])
    tets = np.array([[0, 1, 2, 3]])
    mesh = p4_tetra.make_mesh(points, tets)
    matrix = p4_tetra.stiffness(mesh)
    center = np.array([.2, -.3, .1])
    mode_load, patch = loads(mesh, boundary(points, tets), center)
    polynomial, _ = harmonic_quartic(mesh.points, center)
    assert np.linalg.norm(matrix@polynomial-patch)/np.linalg.norm(patch) < 2e-13
    assert np.max(abs(matrix@mesh.points-mode_load[:, :3])) < 2e-14

    # Independent divergence-theorem load for phi=z on a closed volume:
    # integral(phi*n)=volume*ez, integral(phi*(r x n))=(CB-CG) x force.
    force = np.array([0., 0., 1/6])
    expected = np.r_[force, np.cross(points.mean(axis=0)-center, force)]
    assert np.max(abs(mode_load.T@mesh.points[:, 2]-expected)) < 2e-14
