"""Flat-surface image physics on the same physical canoe used throughout."""

import json

import numpy as np
from numpy.testing import assert_allclose
import pytest

from topwater import SimulationInputs
from topwater.bem import LaplaceOperators, assemble_dirichlet_halfspace, assemble_laplace
from topwater.boundary import BoundarySystem
from topwater.hydrostatics import equilibrate_heave
from topwater.inertia import instantaneous_added_mass
from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels, panel_influence
from .reference_hull import input_data
from .laplace_reference import point_source, triangle_samples


@pytest.fixture(scope='module')
def canoe():
    data = SimulationInputs.from_json(json.dumps(input_data(3)))
    rest, hydro = equilibrate_heave(data.body, data.initial_state)
    return TrianglePanels(hydro.wetted_triangles, -1), rest.position


def test_target_reflection_matches_source_reflection_and_keeps_smooth_diagonal(canoe, tmp_path):
    panels, _ = canoe
    actual = assemble_dirichlet_halfspace(panels)
    base = assemble_laplace(panels)
    triangles = panels.triangles.copy()
    triangles[:, :, 2] *= -1
    # Reflection reverses winding. +1 on reflected triangles is the reflection
    # of the original fluid normal (-1); no source-normal sign is applied twice.
    image = panel_influence(panels.centroids, TrianglePanels(triangles, 1),
                            controls=LaplaceControls(target_block_size=len(panels)))
    assert_allclose(actual.single_layer, base.single_layer-image.single_layer, atol=1e-13, rtol=1e-11)
    assert_allclose(actual.double_layer, base.double_layer-image.double_layer, atol=1e-13, rtol=1e-11)
    assert np.any(np.abs(np.diag(actual.double_layer)) > 1e-4)
    assert not actual.single_layer.flags.writeable
    path = tmp_path/'halfspace.npz'
    actual.save(path)
    restored = LaplaceOperators.load(path)
    assert restored.diagnostics['green_function'] == 'dirichlet_halfspace'
    assert restored.diagnostics['water_level'] == 0
    assert_allclose(restored.double_layer, actual.double_layer, atol=0, rtol=0)


def test_half_jump_and_operator_norm_with_legitimate_image_self_terms(canoe):
    panels, _ = canoe
    op = assemble_dirichlet_halfspace(panels)
    matrix = op.double_layer+.5*np.eye(len(panels))
    q = np.column_stack((panels.normals[:, 0], panels.normals[:, 2]))
    for mask in (np.zeros(len(panels), bool), panels.centroids[:, 0] > 0):
        system = BoundarySystem(op, dirichlet=mask, areas=panels.areas,
                                length_scale=3.5, bounded=False)
        expected = np.linalg.solve(matrix, op.single_layer@q)
        actual = system.solve(potential=expected, normal_velocity=q)
        assert_allclose(actual.potential, expected, atol=1e-12, rtol=1e-11)
        assert_allclose(actual.normal_velocity, q, atol=1e-12, rtol=1e-11)
        assert_allclose(system._h_norm, np.linalg.norm(matrix, 1), rtol=1e-14)
        assert max(actual.diagnostics['original_bie_residual']) < 1e-12


def test_halfspace_requires_panels_below_plane(canoe):
    panels, _ = canoe
    with pytest.raises(ValueError, match='below'):
        assemble_dirichlet_halfspace(panels, water_level=-1)
    with pytest.raises(ValueError, match='water_level'):
        assemble_dirichlet_halfspace(panels, water_level=np.nan)


def test_independent_image_pair_quadrature_near_waterline(canoe):
    panels, _ = canoe
    op = assemble_dirichlet_halfspace(panels)
    row = int(np.argmax(panels.centroids[:, 2]))
    target = panels.centroids[row]
    columns = np.argsort(np.linalg.norm(panels.centroids-target, axis=1))[-4:]
    points, weights = triangle_samples(panels.triangles[columns], 24)
    delta = points-target
    mirrored = target.copy()
    mirrored[2] *= -1
    image_delta = points-mirrored
    r, image_r = np.linalg.norm(delta, axis=2), np.linalg.norm(image_delta, axis=2)
    h_source, h_target = -points[:, :, 2], -target[2]
    # Rationalized difference of inverse distances avoids subtracting nearly
    # equal kernels. Differentiate it directly with respect to source position.
    coefficient = h_target/(np.pi*r*image_r*(r+image_r))
    green = coefficient*h_source
    log_derivative = (delta/r[:, :, None]**2+image_delta/image_r[:, :, None]**2
        +(delta/r[:, :, None]+image_delta/image_r[:, :, None])/(r+image_r)[:, :, None])
    normal = panels.normals[columns]
    double = coefficient*(-normal[:, None, 2]
        -h_source*np.einsum('fqj,fj->fq', log_derivative, normal))
    single_reference = panels.areas[columns]*(green@weights)
    double_reference = panels.areas[columns]*(double@weights)
    assert_allclose(op.single_layer[row, columns], single_reference, atol=1e-13, rtol=1e-9)
    assert_allclose(op.double_layer[row, columns], double_reference, atol=1e-13, rtol=1e-9)
    free_space = panels.areas[columns]*((1/(4*np.pi*r))@weights)
    assert np.max(single_reference/free_space) < .01


def test_source_plus_mirror_bulk_convergence_and_zero_surface_potential():
    source = np.array([0., 0., -.06])
    mirror = source*np.array([1., 1., -1.])
    bulk = np.array([[2., .1, -.1], [0., .8, -.1], [-2., -.1, -.3]])
    surface = np.array([[0., .8, 0.], [2., .1, 0.], [-2., -.1, 0.]])
    targets = np.vstack((bulk, surface))
    exact_phi, exact_velocity = point_source(bulk, source)
    image_phi, image_velocity = point_source(bulk, mirror)
    exact_phi -= image_phi
    exact_velocity -= image_velocity
    errors = []
    for refinement in (2, 3, 4):
        data = SimulationInputs.from_json(json.dumps(input_data(refinement)))
        _, wet = equilibrate_heave(data.body, data.initial_state)
        panels = TrianglePanels(wet.wetted_triangles, -1)
        op = assemble_dirichlet_halfspace(panels)
        _, velocity = point_source(panels.centroids, source)
        _, image_velocity = point_source(panels.centroids, mirror)
        q = np.einsum('frj,fj->fr', velocity-image_velocity, panels.normals)
        solution = BoundarySystem(op, dirichlet=np.zeros(len(panels), bool),
            areas=panels.areas, length_scale=3.5, bounded=False).solve(normal_velocity=q)
        real = panel_influence(targets, panels, gradients=True)
        reflected = targets.copy()
        reflected[:, 2] *= -1
        image = panel_influence(reflected, panels, gradients=True)
        potential = ((real.single_layer-image.single_layer)@q
                     -(real.double_layer-image.double_layer)@solution.potential)
        reflection = np.array([1., 1., -1.])
        gradient = (np.einsum('mnj,nr->mrj', real.gradient_single_layer
                    -image.gradient_single_layer*reflection, q)
                    -np.einsum('mnj,nr->mrj', real.gradient_double_layer
                    -image.gradient_double_layer*reflection, solution.potential))
        assert_allclose(potential[len(bulk):], 0, atol=0, rtol=0)
        assert_allclose(gradient[len(bulk):, :, :2], 0, atol=0, rtol=0)
        errors.append([np.linalg.norm(potential[:len(bulk)]-exact_phi)/np.linalg.norm(exact_phi),
                       np.linalg.norm(gradient[:len(bulk)]-exact_velocity)/np.linalg.norm(exact_velocity)])
    assert np.all(np.diff(errors, axis=0) < 0)


def test_inertia_yaw_translation_density_and_moment_reference_transform(canoe):
    panels, reference = canoe
    baseline = instantaneous_added_mass(panels, reference)
    angle = .37
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1.]])
    transform = np.zeros((6, 6))
    transform[:3, :3] = transform[3:, 3:] = rotation
    shift = np.array([2., -3., .6])
    moved = TrianglePanels(panels.triangles@rotation.T+shift, -1)
    result = instantaneous_added_mass(moved, reference@rotation.T+shift, water_level=shift[2], density=1025.)
    assert_allclose(result.added_mass, 1.025*transform@baseline.added_mass@transform.T,
                    atol=1e-10, rtol=1e-10)
    # Moving the rotation reference changes unit-mode velocities, hence A by a
    # congruence transform. It must not change any physical pressure solution.
    displacement = np.array([.2, -.1, .05])
    x, y, z = displacement
    cross = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    change = np.eye(6)
    change[:3, 3:] = cross
    relocated = instantaneous_added_mass(panels, reference+displacement)
    assert_allclose(relocated.added_mass, change.T@baseline.added_mass@change, atol=1e-11, rtol=1e-10)
    size = 1.7
    resized = instantaneous_added_mass(TrianglePanels(size*panels.triangles, -1), size*reference)
    powers = np.array([0, 0, 0, 1, 1, 1])
    assert_allclose(resized.added_mass, baseline.added_mass*size**(3+powers[:, None]+powers[None]),
                    atol=1e-10, rtol=1e-10)
    assert not baseline.added_mass.flags.writeable
    assert not baseline.acceleration_potential.flags.writeable
    assert min(baseline.diagnostics['normalized_energy_eigenvalues']) > 0
