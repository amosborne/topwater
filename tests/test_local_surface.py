"""Analytic volume moments independently verify manufactured boundary loads."""
import numpy as np
from numpy.testing import assert_allclose

from topwater.constants import WATER_DENSITY
from tests.experiments.free_surface.boat_local_surface import exact_data, field


def test_tetrahedron_boundary_work_and_six_pressure_moments():
    points = np.vstack((np.zeros(3), np.eye(3)))
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    beam = 1.7
    load, moments, energy = exact_data(points, faces, [], np.arange(4), np.zeros(3), beam)
    volume, first = 1/6, np.full(3, 1/24)
    second = (np.ones((3, 3))+np.eye(3))/120
    g = np.array([.2, .3, .4])/beam
    hessians = np.array([[[0, 1, 0], [1, 0, 0], [0, 0, 0]],
                        [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
                        [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
                        [[2, 0, 0], [0, -2, 0], [0, 0, 0]],
                        [[2, 0, 0], [0, 0, 0], [0, 0, -2]]])/beam**2
    expected = np.empty((6, 6))
    expected[:3, 0] = volume*g
    expected[3:, 0] = np.cross(first, g)
    energies = [volume*(g@g)]
    for column, h in enumerate(hessians, start=1):
        expected[:3, column] = h@first
        expected[3:, column] = sum(second[i, j]*np.cross(np.eye(3)[i], h[:, j])
                                    for i in range(3) for j in range(3))
        energies.append(np.trace(h.T@h@second))
    assert_allclose(moments, -WATER_DENSITY*expected, rtol=2e-13, atol=1e-12)
    assert_allclose(energy, energies, rtol=2e-13, atol=1e-14)
    assert_allclose(load.sum(axis=0), 0, atol=1e-15)


def test_harmonic_fields_have_independent_potential_gradients():
    points = np.array([[.2, -.7, .3], [.8, .1, -.6]])
    center, beam = np.array([.1, .2, -.3]), 2.4
    _, actual = field(points, center, beam)
    for axis in range(3):
        delta = np.eye(3)[axis]*1e-5
        derivative = (field(points+delta, center, beam)[0]-field(points-delta, center, beam)[0])/2e-5
        assert_allclose(actual[..., axis], derivative, rtol=2e-10, atol=2e-12)
