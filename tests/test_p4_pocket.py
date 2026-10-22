"""Independent polynomial and physical-field checks of the contained P4 probe."""
import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from tests.laplace_reference import triangle_samples
from .experiments.free_surface import p4_tetra
from .experiments.free_surface.p4_pocket import solve_case


def test_exact_tetrahedral_polynomial_energy():
    points = np.array([[0., 0., 0.], [1., 0., 0.],
                       [0., 1., 0.], [0., 0., 1.]])
    mesh = p4_tetra.make_mesh(points, np.array([[0, 1, 2, 3]]))
    matrix = p4_tetra.stiffness(mesh)
    potential = np.sum(mesh.points**2, axis=1)
    # Integral over the unit simplex of 4*(x²+y²+z²) is 1/5.
    assert potential@(matrix@potential) == pytest.approx(.2, abs=2e-14)
    affine = mesh.points@np.array([.7, -.4, .2])
    assert affine@(matrix@affine) == pytest.approx(.69/6, abs=2e-14)
    assert np.max(abs(matrix@np.ones(len(mesh.points)))) < 1e-13
    values, _ = p4_tetra.basis(p4_tetra.indices(4)/4)
    assert np.max(abs(values-np.eye(35))) < 2e-14


def test_interface_restriction_matches_the_whole_triangle():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
                       [0., 0., 1.], [0., 0., -1.]])
    mesh = p4_tetra.make_mesh(points, np.array([[0, 1, 2, 3], [0, 2, 1, 4]]))
    interface = np.array([[0, 1, 2]])
    restriction, retained = p4_tetra.restrict_trace(mesh, interface)
    # Unrelated retained volume/other-face coefficients must not leak into the
    # constrained interface trace. Compare interiors, not merely shared nodes.
    coefficient = np.random.default_rng(981).normal(size=len(retained))
    full = restriction@coefficient
    barycentric, _ = triangle_samples(np.eye(3)[None], 9)
    barycentric = barycentric[0]
    values, _ = p4_tetra.basis(barycentric)
    face_nodes = p4_tetra.face_nodes(mesh, interface)[0]
    assert np.max(abs(values@full[face_nodes]-barycentric@full[interface[0]])) < 2e-14
    reversed_nodes = p4_tetra.face_nodes(mesh, interface[:, ::-1])[0]
    assert np.array_equal(np.sort(face_nodes), np.sort(reversed_nodes))
    assert np.max(abs(restriction@np.ones(len(retained))-1)) < 1e-15


@pytest.fixture(scope='module')
def physical_refinement():
    with threadpool_limits(1):
        return {(n, coupled): solve_case(n, coupled)
                for n in (6, 12) for coupled in (False, True)}


def test_physical_harmonic_field_refines_without_claiming_qualification(physical_refinement):
    for coupled in (False, True):
        coarse, fine = (physical_refinement[n, coupled] for n in (6, 12))
        for name in ('boundary_potential_relative_rms',
                     'whole_free_surface_flux_relative_rms',
                     'projected_weak_free_surface_flux_relative_error',
                     'energy_relative_error'):
            assert fine[name] < .5*coarse[name], (coupled, name, coarse[name], fine[name])
        assert fine['relative_unknown_residual'] < 1e-11
        assert fine['qualified'] is False
    # Physical refinement and interface algebra are separate evidence.
    fine_mixed = physical_refinement[12, True]
    assert fine_mixed['shared_interface_potential_max_abs'] < 1e-12
    assert fine_mixed['relative_shared_weak_flux_residual'] < 1e-11
    assert physical_refinement[12, False]['whole_free_surface_flux_relative_rms'] < .001
