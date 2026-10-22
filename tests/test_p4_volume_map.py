"""Coordinate-map work identities before coupling to physical canoe dynamics."""
from itertools import permutations

import numpy as np

from .experiments.free_surface.p4_tetra import make_mesh
from .experiments.free_surface.p4_volume_map import make_map
from .experiments.free_surface.p4_geometry import integrals


def cube_mesh():
    points = np.array(list(np.ndindex(2, 2, 2)), dtype=float)
    numbering = np.arange(8).reshape(2, 2, 2)
    paths = np.zeros((6, 4, 3), dtype=int)
    paths[:, 1:] = np.cumsum(np.eye(3, dtype=int)[list(permutations(range(3)))], axis=1)
    tets = np.array([numbering[tuple(path.T)] for path in paths])
    # Nonlinear cell calculus uses positive reference orientations.
    negative = np.linalg.det(points[tets[:, 1:]]-points[tets[:, :1]]) < 0
    tets[negative, 0], tets[negative, 1] = tets[negative, 1].copy(), tets[negative, 0].copy()
    return make_mesh(points, tets)


def test_boundary_extension_translation_and_exact_transpose_work():
    mesh = cube_mesh()
    mapping = make_map(mesh)
    assert len(mapping.interior) > 0
    shift = np.array([.03, -.02, .01])
    translated = mapping.lift(np.tile(shift, (len(mapping.boundary), 1)))
    np.testing.assert_allclose(translated, np.tile(shift, (len(mesh.points), 1)), atol=2e-16)
    rng = np.random.default_rng(147)
    velocity = rng.normal(size=(len(mapping.boundary), 3))
    force = rng.normal(size=mesh.points.shape)
    np.testing.assert_allclose(np.sum(force*mapping.lift(velocity)),
                               np.sum(mapping.pullback(force)*velocity), atol=2e-13, rtol=2e-14)


def test_actual_curved_fluid_energy_gradient_pulls_back_to_moving_boundary():
    mesh = cube_mesh()
    mapping = make_map(mesh)
    points = mesh.points[mapping.boundary]
    displacement = np.zeros_like(points)
    displacement[:, 2] = .04*points[:, 2]*np.sin(np.pi*points[:, 0])*np.sin(np.pi*points[:, 1])
    coordinates = mesh.points+mapping.lift(displacement)
    phi = np.sin(mesh.points@np.array([.7, -.5, .4]))
    gradient = np.zeros_like(mesh.points)
    energy = 0.
    for nodes in mesh.elements:
        result = integrals(coordinates[nodes], phi[nodes])
        energy += result['kinetic_energy']
        np.add.at(gradient, nodes, result['kinetic_gradient'])
    # A genuine top-boundary bubble moves volume nodes while all side and
    # bottom coordinates stay fixed. Nodal phi remains fixed in every solve.
    direction = np.zeros_like(points)
    direction[:, 2] = points[:, 2]*points[:, 0]*(1-points[:, 0])*points[:, 1]*(1-points[:, 1])
    rate = mapping.lift(direction)
    assert np.linalg.norm(rate[mapping.interior]) > 0
    expected = np.sum(mapping.pullback(gradient)*direction)
    errors = []
    for step in (2e-3, 1e-3):
        values = []
        for sign in (-1, 1):
            current = coordinates+sign*step*rate
            values.append(sum(integrals(current[nodes], phi[nodes])['kinetic_energy']
                              for nodes in mesh.elements))
        errors.append(abs((values[1]-values[0])/(2*step)-expected))
    assert errors[-1] < 2e-8*max(1., abs(expected))
    assert errors[-1] < .4*errors[0]+1e-9
