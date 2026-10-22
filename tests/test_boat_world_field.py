"""Physical affine harmonic velocity on an inclined, partially exposed graph."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

from .experiments.free_surface.boat_world_field import surface_rates


@pytest.mark.parametrize('clipped', [False, True])
def test_affine_harmonic_rates_on_full_or_partly_exposed_parent(clipped):
    parent_xy = np.array([[0., 0.], [2., 0.], [0., 2.]])
    if clipped:
        # Remove a finite corner from the parent; integrate only exposed water.
        xy = np.array([[.5, 0.], [2., 0.], [0., 2.], [0., .5]])
        faces = np.array([[0, 1, 2], [0, 2, 3]])
    else:
        xy, faces = parent_xy, np.array([[0, 1, 2]])
    slope, height = np.array([.2, -.13]), .3
    horizontal, vertical = np.array([.7, -.4]), .23
    points = np.column_stack((xy, xy@slope+height))
    mapping = np.column_stack((1-xy.sum(axis=1)/2, xy/2))
    parent_height = parent_xy@slope+height
    psi = parent_xy@(horizontal+vertical*slope)+vertical*height+1.9
    exact_flux = vertical-horizontal@slope
    flux = np.zeros(len(points))
    for face in faces:
        triangle = xy[face]
        area = np.linalg.det(triangle[1:]-triangle[:1])/2
        flux[face] += area/3*exact_flux
    eta_rate, psi_rate, diagnostics = surface_rates(points, faces, np.arange(len(points)),
        mapping, psi, flux, gravity=9.4)
    exact_psi_rate = (-9.4*parent_height-.5*(horizontal@horizontal)
                      -vertical*(horizontal@slope)+.5*vertical**2)
    assert_allclose(eta_rate, exact_flux, atol=3e-15)
    assert_allclose(psi_rate, exact_psi_rate, atol=3e-14)
    assert diagnostics['unresolved_kinematic_mass_relative'] < 2e-14


def test_unresolved_local_flux_is_not_hidden_from_bernoulli():
    # One flat parent with an interior cut-node degree of freedom. The local
    # flux cannot fit the parent affine state. A zero-mean local motion still
    # contributes positive vertical kinetic head to Bernoulli.
    points = np.array([[0., 0., 0.], [2., 0., 0.], [0., 2., 0.], [.5, .5, 0.]])
    faces = np.array([[0, 1, 3], [1, 2, 3], [2, 0, 3]])
    mapping = np.column_stack((1-points[:, :2].sum(axis=1)/2, points[:, :2]/2))
    # Construct a flux orthogonal to every parent trace, without assuming a
    # nodal mass formula or reproducing the rate assembly.
    _, _, vh = np.linalg.svd(mapping.T, full_matrices=True)
    flux = vh[-1]
    eta_rate, psi_rate, diagnostics = surface_rates(points, faces, np.arange(4),
        mapping, np.zeros(3), flux)
    assert_allclose(eta_rate, 0., atol=1e-14)
    assert diagnostics['unresolved_kinematic_mass_relative'] > .999999
    # Integral of Bernoulli over the parent is positive. Its nodal projection
    # may have negative entries; the equal parent basis integrals give the mean.
    assert psi_rate.sum() > 0
