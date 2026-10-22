"""Actual whole-canoe geometry proofs; optional mesher, no fluid operator."""
from importlib.util import find_spec
import json

import numpy as np
import pytest
from scipy.sparse import eye
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

if find_spec('meshpy') is None:
    pytest.skip('optional MeshPy preprocessing required', allow_module_level=True)

from tests.laplace_reference import triangle_samples
from .experiments.free_surface import p4_tetra
from .experiments.free_surface.p4_whole_geometry import P4WholeCanoeGeometry


def moved(orientation, omega, time):
    return (Rotation.from_rotvec(time*np.asarray(omega))*Rotation.from_quat(
        orientation, scalar_first=True)).as_quat(scalar_first=True)


@pytest.fixture(scope='module')
def case():
    with threadpool_limits(1):
        owner = P4WholeCanoeGeometry.make()
    print(json.dumps(dict(stage='geometry graph ready', **owner.preflight)), flush=True)
    return owner


def finite_state(owner):
    rest = owner.domain.rest
    position = rest.position+np.array([.003, -.002, .005])
    orientation = moved(rest.orientation, [.008, -.006, .004], 1.)
    xy = owner.mesh.points[owner.free_nodes, :2]
    eta = 8e-5*np.cos(1.3*xy[:, 0])+3e-5*np.sin(2*xy[:, 1])
    high = owner.chart.high_contact_nodes
    columns = owner.chart_eta_columns[np.searchsorted(owner.chart.free_nodes, high)]
    eta[columns] += 2e-6*np.sin(np.arange(len(columns)))
    return position, orientation, eta


def test_reference_topology_positive_orientation_and_single_boundary_ownership(case):
    owner = case
    state = owner.evaluate(owner.domain.rest.position, owner.domain.rest.orientation,
                           np.zeros(len(owner.free_nodes)))
    np.testing.assert_allclose(state.points, owner.mesh.points, atol=3e-15, rtol=0)
    determinant = np.linalg.det(owner.mesh.points[owner.mesh.tets[:, 1:]]
                                -owner.mesh.points[owner.mesh.tets[:, :1]])
    assert np.all(determinant > 0)
    assert len(np.unique(owner.chart_to_global)) == len(owner.chart_to_global)
    interface_only = np.setdiff1d(np.arange(len(owner.chart.reference_points)), owner.override_chart_nodes)
    assert len(interface_only) > 0
    assert not np.any(np.isin(owner.chart_to_global[interface_only], owner.volume_map.boundary))
    height_identity = state.boundary_height_jacobians[2][owner.free_boundary_rows]-eye(len(owner.free_nodes))
    assert height_identity.nnz == 0 or np.max(abs(height_identity.data)) == 0
    assert np.max(abs(state.boundary_body_jacobian[owner.free_boundary_rows, 2])) == 0
    # A new waterline height physically changes its own node, not an interpolant
    # from the original P1 controls, and has one owner shared by both surfaces.
    high = owner.chart.high_contact_nodes[len(owner.chart.high_contact_nodes)//2]
    column = owner.chart_eta_columns[np.searchsorted(owner.chart.free_nodes, high)]
    eta = np.zeros(len(owner.free_nodes)); eta[column] = 2e-6
    changed = owner.evaluate(owner.domain.rest.position, owner.domain.rest.orientation, eta)
    assert changed.points[owner.chart_to_global[high], 2] == eta[column]


def test_all_six_world_twists_and_global_heights_match_geometry_fd(case):
    owner = case
    position, orientation, eta = finite_state(owner)
    state = owner.evaluate(position, orientation, eta)
    np.testing.assert_allclose(state.points[owner.free_nodes, 2], eta, atol=1e-18, rtol=0)
    height_rate = 3e-4*np.cos(np.arange(len(eta))*.9)
    directions = [(np.eye(6)[i], np.zeros(len(eta))) for i in range(6)]
    directions.append((np.array([.07, -.04, .03, .02, -.03, .015]), height_rate))
    with threadpool_limits(1):
        for twist, rate in directions:
            expected = state.velocity(rate, twist)
            errors = []
            for step in (1e-5, 5e-6):
                plus = owner.evaluate(position+step*twist[:3], moved(orientation, twist[3:], step), eta+step*rate)
                minus = owner.evaluate(position-step*twist[:3], moved(orientation, twist[3:], -step), eta-step*rate)
                errors.append(np.max(abs((plus.points-minus.points)/(2*step)-expected)))
            assert errors[-1] < 2e-8, (twist, errors)
            assert errors[-1] < .8*errors[0]+2e-10, (twist, errors)
        force = np.random.default_rng(147).normal(size=state.points.shape)
        twist, rate = directions[-1]
        body_force, height_force = state.pullback(force)
        np.testing.assert_allclose(np.sum(force*state.velocity(rate, twist)),
                                   body_force@twist+height_force@rate, atol=2e-11, rtol=2e-13)


def test_boundary_transport_matches_full_map_for_actual_pressure_work(case):
    from .experiments.free_surface.canonical_canoe import body_map
    from .experiments.free_surface.p4_variational_canoe import (
        _configuration_state, _pullback, _boundary_configuration_pullback)

    owner = case;position, orientation, eta = finite_state(owner)
    state = owner.evaluate(position, orientation, eta)
    boundary = owner.volume_map.boundary
    twist = np.array([.07, -.04, .03, .02, -.03, .015])
    height_rate = 3e-4*np.cos(np.arange(len(eta))*.9)
    with threadpool_limits(1):
        direct_velocity = state.boundary_velocity(height_rate, twist)
        np.testing.assert_allclose(direct_velocity, state.velocity(height_rate, twist)[boundary],
                                   atol=3e-15, rtol=3e-14)
        # Independent physical pressure work p=rho*g*z on the complete real
        # boundary, including hull and every free-surface height master.
        force = owner.boundary_moments(state)['gravity_gradient']
        direct = state.boundary_pullback(force[boundary])
        full = state.pullback(force)
        for actual, expected in zip(direct, full):
            np.testing.assert_allclose(actual, expected, atol=3e-10, rtol=3e-14)
        work = np.sum(force[boundary]*direct_velocity)
        assert work == pytest.approx(direct[0]@twist+direct[1]@height_rate, abs=3e-10, rel=3e-14)
        rotation = (Rotation.from_quat(orientation, scalar_first=True)
                    *Rotation.from_quat(owner.domain.rest.orientation, scalar_first=True).inv()).as_rotvec()
        q = np.r_[position, rotation, eta]
        qdot = np.r_[np.linalg.solve(body_map(rotation), twist), height_rate]
        for step in (1e-5, 5e-6):
            plus, minus = q+step*qdot, q-step*qdot
            old = (_pullback(_configuration_state(owner, plus), plus, force)
                   -_pullback(_configuration_state(owner, minus), minus, force))/(2*step)
            new = (_boundary_configuration_pullback(owner, plus, force[boundary])
                   -_boundary_configuration_pullback(owner, minus, force[boundary]))/(2*step)
            np.testing.assert_allclose(new, old, atol=3e-8, rtol=3e-13)


def moment_difference(owner, plus, minus):
    """Independent physical quadrature, subtracting integrands before totals."""
    barycentric, weights = triangle_samples(np.eye(3)[None], 10)
    values, derivative = p4_tetra.basis(barycentric[0])
    derivative = derivative[:, :, 1:]-derivative[:, :, :1]
    nodes = p4_tetra.face_nodes(owner.mesh, owner.faces)
    terms = []
    for points in (plus.points, minus.points):
        coordinates = points[nodes]
        tangent = np.einsum('qaj,fai->fqji', derivative, coordinates)
        normal_z = np.cross(tangent[:, :, 0], tangent[:, :, 1])[:, :, 2]
        z = values@coordinates[:, :, 2, None]
        z = z[:, :, 0]
        terms.append((z*normal_z, 4905*z*z*normal_z))
    return tuple(.5*np.einsum('q,fq->', weights, a-b) for a, b in zip(terms[0], terms[1]))


def test_current_closed_physical_volume_gravity_and_shape_work(case):
    owner = case
    position, orientation, eta = finite_state(owner)
    state = owner.evaluate(position, orientation, eta)
    with threadpool_limits(1):
        surface = owner.boundary_moments(state)
        volume = owner.volume_moments(state)
    assert volume['minimum_sampled_jacobian'] > 0
    assert abs(volume['volume']/surface['volume']-1) < 2e-12
    assert abs(volume['gravity']/surface['gravity']-1) < 2e-12
    twist = np.array([.07, -.04, .03, .02, -.03, .015])
    rate = 3e-4*np.cos(np.arange(len(eta))*.9)
    expected = []
    with threadpool_limits(1):
        for key in ('volume_gradient', 'gravity_gradient'):
            body_force, height_force = state.pullback(surface[key])
            expected.append(body_force@twist+height_force@rate)
        errors = []
        for step in (2e-4, 1e-4):
            plus = owner.evaluate(position+step*twist[:3], moved(orientation, twist[3:], step), eta+step*rate)
            minus = owner.evaluate(position-step*twist[:3], moved(orientation, twist[3:], -step), eta-step*rate)
            actual = np.asarray(moment_difference(owner, plus, minus))/(2*step)
            errors.append(abs(actual-expected))
    assert errors[-1][0] < 2e-8, errors
    assert errors[-1][1] < 2e-5, errors
    assert np.all(errors[-1] < .5*errors[0]+np.array([1e-10, 1e-7])), errors
