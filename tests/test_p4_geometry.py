"""Curved physical cell identities before nonlinear canoe integration."""
import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from .experiments.free_surface import p4_tetra
from .experiments.free_surface.p4_geometry import integrals, InvalidGeometryError
from .experiments.free_surface.pocket import boundary


def curved_cell(amount):
    vertices = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    mesh = p4_tetra.make_mesh(vertices, np.array([[0, 1, 2, 3]]))
    coordinates = mesh.points.copy()
    coordinates[:, 2] *= 1+amount*coordinates[:, 0]
    return mesh, coordinates


def boundary_moments(mesh, coordinates):
    # Independent divergence theorem: V=int z*nz, int_volume z=.5*int z²*nz.
    x, w = leggauss(9)
    x, w = (x+1)/2, w/2
    u, v = np.meshgrid(x, x, indexing='ij')
    wu, wv = np.meshgrid(w, w, indexing='ij')
    r, s = u.ravel(), ((1-u)*v).ravel()
    values, derivative = p4_tetra.basis(np.column_stack((1-r-s, r, s)))
    gradient = derivative[:, :, 1:]-derivative[:, :, :1]
    surface = coordinates[p4_tetra.face_nodes(mesh, boundary(mesh.points, mesh.tets))]
    height = np.einsum('qa,fa->fq', values, surface[:, :, 2])
    tangent = np.einsum('qaj,fai->fqji', gradient, surface)
    nz = np.cross(tangent[:, :, 0], tangent[:, :, 1])[:, :, 2]
    weights = (wu*wv*(1-u)).ravel()
    return np.sum(height*nz*weights), .5*np.sum(height**2*nz*weights)


def test_curved_cell_exact_physical_flow_volume_and_gravity():
    amount = .7
    mesh, coordinates = curved_cell(amount)
    nodes = mesh.elements[0]
    expected_volume = 1/6+amount/24
    expected_first_moment = 1/24+amount/60+amount**2/360
    for axis in range(3):
        result = integrals(coordinates[nodes], coordinates[nodes, axis])
        assert result['volume'] == pytest.approx(expected_volume, abs=2e-14)
        assert result['gravity'] == pytest.approx(9810*expected_first_moment, abs=2e-10)
        # phi is exactly the physical coordinate, despite the curved map.
        assert result['kinetic_energy'] == pytest.approx(500*expected_volume, abs=2e-10)
    nonlinear = integrals(coordinates[nodes], np.prod(coordinates[nodes, :2], axis=1))
    assert nonlinear['kinetic_energy'] == pytest.approx(500*(1/30+amount/90), abs=2e-11)
    np.testing.assert_allclose(boundary_moments(mesh, coordinates),
                               [expected_volume, expected_first_moment], rtol=2e-13, atol=2e-14)


def test_curved_harmonic_flow_has_exact_physical_energy_and_shape_work():
    mesh, coordinates = curved_cell(0.)
    reference = coordinates.copy()
    alpha = .2
    shape = reference[:, 0]**2-reference[:, 1]**2+reference[:, 0]
    coordinates[:, 2] += alpha*shape
    nodes = mesh.elements[0]
    # X=(r,s,t+alpha*(r²-s²+r)), detJ=1. Nodal phi=t is held
    # fixed, giving the nonuniform harmonic physical velocity
    # (-alpha*(2x+1), 2alpha*y, 1). All closed forms are continuum integrals.
    result = integrals(coordinates[nodes], reference[nodes, 2])
    assert result['volume'] == pytest.approx(1/6, abs=2e-14)
    assert result['gravity'] == pytest.approx(9810*(1+alpha)/24, abs=2e-10)
    assert result['kinetic_energy'] == pytest.approx(500*(1/6+7*alpha**2/15), abs=2e-11)
    np.testing.assert_allclose(boundary_moments(mesh, coordinates),
                               [1/6, (1+alpha)/24], rtol=2e-13, atol=2e-14)
    assert np.dot(result['volume_gradient'][:, 2], shape[nodes]) == pytest.approx(0., abs=2e-14)
    assert np.dot(result['gravity_gradient'][:, 2], shape[nodes]) == pytest.approx(9810/24, abs=2e-10)
    assert np.dot(result['kinetic_gradient'][:, 2], shape[nodes]) == pytest.approx(7000*alpha/15, abs=2e-10)


def test_curved_geometry_fixed_potential_shape_work_and_quadrature_refinement():
    mesh, coordinates = curved_cell(.7)
    nodes = mesh.elements[0]
    x = coordinates[nodes]
    reference = p4_tetra.indices(4)/4
    # Curved geometry, generic fixed P4 potential, and a genuine high-order
    # shape direction: finite differences change geometry but never phi.
    potential = np.sin(reference@np.array([.2, 1.1, -.7, .9]))
    direction = np.column_stack((reference[:, 1]*reference[:, 2],
        reference[:, 2]*reference[:, 3], np.prod(reference, axis=1)))
    base = integrals(x, potential, order=10)
    for name in ('kinetic', 'volume', 'gravity'):
        scalar = 'kinetic_energy' if name == 'kinetic' else name
        expected = np.sum(base[name+'_gradient']*direction)
        errors = []
        for step in (1e-3, 5e-4):
            actual = (integrals(x+step*direction, potential, order=10)[scalar]
                      -integrals(x-step*direction, potential, order=10)[scalar])/(2*step)
            errors.append(abs(actual-expected))
        assert errors[-1] < 1e-7*max(1., abs(expected)), (name, errors, expected)
        assert errors[-1] < .4*errors[0]+1e-9, (name, errors)
    # Rigid translation cannot alter kinetic energy or volume; raising the
    # complete cell increases physical gravity by rho*g*V.
    np.testing.assert_allclose(base['kinetic_gradient'].sum(axis=0), 0., atol=2e-11)
    np.testing.assert_allclose(base['volume_gradient'].sum(axis=0), 0., atol=2e-14)
    np.testing.assert_allclose(base['gravity_gradient'].sum(axis=0),
                               [0., 0., 9810*base['volume']], atol=2e-10)
    reference_integrals = integrals(x, potential, order=12)
    energy_errors, gradient_errors = [], []
    for order in (4, 6, 8):
        result = integrals(x, potential, order=order)
        energy_errors.append(abs(result['kinetic_energy']/reference_integrals['kinetic_energy']-1))
        gradient_errors.append(np.linalg.norm(result['kinetic_gradient']-reference_integrals['kinetic_gradient'])
                               /np.linalg.norm(reference_integrals['kinetic_gradient']))
    assert energy_errors[2] < energy_errors[1] < energy_errors[0]
    assert gradient_errors[2] < gradient_errors[1] < gradient_errors[0]
    assert energy_errors[2] < 1e-10
    assert gradient_errors[2] < 1e-8
    shifted = integrals(x, potential+.37, order=10)
    assert shifted['kinetic_energy'] == pytest.approx(base['kinetic_energy'], abs=2e-11)
    np.testing.assert_allclose(shifted['kinetic_gradient'], base['kinetic_gradient'], atol=2e-10, rtol=2e-12)


def test_curved_cell_sampled_inversion_is_rejected():
    mesh, coordinates = curved_cell(-3.)
    with pytest.raises(InvalidGeometryError, match='Jacobian'):
        integrals(coordinates[mesh.elements[0]], np.zeros(35))
