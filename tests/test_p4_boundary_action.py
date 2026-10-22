"""Physical work, curved geometry partials and objectivity, without a PDE solve."""
import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from tests.laplace_reference import triangle_samples
from .experiments.free_surface.p4_tetra import indices, basis, make_mesh, face_nodes
from .experiments.free_surface.pocket import boundary
from .experiments.free_surface.p4_boundary_action import integrals


def test_plane_normal_cross_signs_and_exact_physical_action():
    barycentric = indices(3)/4
    r, s = barycentric[:, 1], barycentric[:, 2]
    a, b = .3, -.2
    coordinates = np.column_stack((r, s, .4+a*r+b*s))
    potential = .7+.2*r-.1*s
    velocity = np.tile([.4, -.3, .7], (15, 1))
    result = integrals(coordinates, potential, velocity)
    integral_phi = .7/2+(.2-.1)/6
    normal = np.array([-a, -b, 1.])
    assert result['action'] == pytest.approx(integral_phi*normal@velocity[0], abs=2e-14)
    np.testing.assert_allclose(result['velocity_gradient'].sum(axis=0), normal*integral_phi, atol=2e-14)
    # Tilting z by +epsilon*x rotates the area vector toward -x, and
    # tilting by +epsilon*y rotates it toward -y. These fix both cross signs.
    assert result['position_gradient'][:, 2]@r == pytest.approx(-.4*integral_phi, abs=2e-14)
    assert result['position_gradient'][:, 2]@s == pytest.approx(.3*integral_phi, abs=2e-14)


def test_curved_independent_geometry_velocity_partials_and_rigid_objectivity():
    barycentric = indices(3)/4
    r, s = barycentric[:, 1], barycentric[:, 2]
    bubble = 27*np.prod(barycentric, axis=1)
    coordinates = np.column_stack((r+.04*bubble, s-.025*bubble*(r-s),
                                   .2+.3*r*r-.2*r*s+.1*s**3))
    potential = np.sin(np.arange(15)*.7)+.3
    velocity = np.column_stack((.2+r*s, -.1+r*r-s, .4+.2*s*s+r*bubble))
    result = integrals(coordinates, potential, velocity)
    rng = np.random.default_rng(593)
    direction = rng.normal(size=(15, 3))*.1
    expected = np.sum(result['position_gradient']*direction)
    errors = []
    for step in (2e-3, 1e-3):
        actual = (integrals(coordinates+step*direction, potential, velocity)['action']
                  -integrals(coordinates-step*direction, potential, velocity)['action'])/(2*step)
        errors.append(abs(actual-expected))
    # F is quadratic in X at fixed U: central differences are exact except
    # roundoff. Do not demand fictitious O(h²) truncation for this partial.
    assert max(errors) < 2e-12
    expected_velocity = np.sum(result['velocity_gradient']*direction)
    for step in (2e-3, 1e-3):
        actual = (integrals(coordinates, potential, velocity+step*direction)['action']
                  -integrals(coordinates, potential, velocity-step*direction)['action'])/(2*step)
        assert abs(actual-expected_velocity) < 2e-12
    rotation = Rotation.from_rotvec([.3, -.4, .2]).as_matrix()
    transformed = integrals(coordinates@rotation.T+[.1, .2, -.3], potential, velocity@rotation.T)
    assert transformed['action'] == pytest.approx(result['action'], abs=2e-14)
    for key in ('position_gradient', 'velocity_gradient'):
        np.testing.assert_allclose(transformed[key], result[key]@rotation.T, atol=2e-14, rtol=2e-12)
    np.testing.assert_allclose(result['position_gradient'].sum(axis=0), 0, atol=2e-14)
    refined = integrals(coordinates, potential, velocity, order=12)
    for key in ('action', 'position_gradient', 'velocity_gradient'):
        np.testing.assert_allclose(refined[key], result[key], atol=2e-14, rtol=2e-12)


def test_closed_curved_cell_constant_potential_is_physical_volume_work():
    mesh = make_mesh(np.array([[0., 0., 0.], [1., 0., 0.],
                               [0., 1., 0.], [0., 0., 1.]]), np.array([[0, 1, 2, 3]]))
    coordinates = mesh.points.copy()
    coordinates[:, 2] += .2*(coordinates[:, 0]**2-coordinates[:, 1]**2+coordinates[:, 0])
    faces = face_nodes(mesh, boundary(mesh.points, mesh.tets))
    # Isotropic dilation has U=X and dV/dt=3V; constant phi is not a
    # zero-action gauge on an isolated flux-carrying patch or expanding cell.
    action = sum(integrals(coordinates[nodes], np.ones(15), coordinates[nodes])['action'] for nodes in faces)
    assert action == pytest.approx(.5, abs=3e-14)

    def volume(points):
        barycentric, weights = triangle_samples(np.eye(3)[None], 10)
        values, derivative = basis(barycentric[0])
        derivative = derivative[:, :, 1:]-derivative[:, :, :1]
        patches = points[faces]
        height = np.einsum('qa,fa->fq', values, patches[:, :, 2])
        tangent = np.einsum('qaj,fai->fqji', derivative, patches)
        nz = np.cross(tangent[:, :, 0], tangent[:, :, 1])[:, :, 2]
        return .5*np.einsum('fq,fq,q->', height, nz, weights)

    assert volume(coordinates) == pytest.approx(1/6, abs=2e-14)
    errors = []
    for step in (2e-3, 1e-3):
        actual = (volume((1+step)*coordinates)-volume((1-step)*coordinates))/(2*step)
        errors.append(abs(actual-action))
    assert errors[-1] < 2e-7
    assert errors[-1] < .3*errors[0]


def test_zero_area_samples_fail_without_rejecting_vertical_or_downward_faces():
    coordinates = np.column_stack((indices(3)[:, 1:]/4, np.zeros(15)))
    velocity = np.ones((15, 3))
    vertical = coordinates[:, [0, 2, 1]]
    downward = coordinates.copy(); downward[:, 0] *= -1
    assert np.isfinite(integrals(vertical, np.ones(15), velocity)['action'])
    assert np.isfinite(integrals(downward, np.ones(15), velocity)['action'])
    with pytest.raises(ValueError, match='boundary Jacobian'):
        integrals(np.zeros((15, 3)), np.ones(15), velocity)
