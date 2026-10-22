"""Actual canonical hull's P4 boundary chart; no fluid operator or dynamics."""
import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from topwater.frames import rotation_matrix
from .experiments.free_surface.moving_canoe import make_domain
from .experiments.free_surface.p4_canoe_contact import P4CanoeContact
from .experiments.free_surface.canoe_contact import ContactTopologyError
from .experiments.free_surface.p4_tetra import basis
from tests.laplace_reference import triangle_samples


@pytest.fixture(scope='module')
def case():
    domain = make_domain(12, False, 2)
    return domain, P4CanoeContact.make(domain.contact)


def moved(orientation, omega, time):
    return (Rotation.from_rotvec(time*np.asarray(omega))*Rotation.from_quat(
        orientation, scalar_first=True)).as_quat(scalar_first=True)


def test_reference_is_exact_topological_lift_and_all_surface_heights_are_independent(case):
    domain, chart = case
    initial = chart.evaluate(domain.rest.position, domain.rest.orientation,
                             np.zeros(len(chart.free_nodes)))
    np.testing.assert_allclose(initial.points, chart.reference_points, atol=4e-16, rtol=0)
    assert np.max(abs(np.asarray(chart.lift.sum(axis=1)).ravel()-1)) == 0
    assert len(chart.free_nodes) > len(domain.contact.free_nodes)
    assert len(chart.high_contact_nodes) > 0
    identity = initial.height_jacobians[2][chart.free_nodes].toarray()
    np.testing.assert_array_equal(identity, np.eye(len(chart.free_nodes)))
    np.testing.assert_array_equal(initial.body_jacobian[chart.free_nodes, 2], 0)
    # An independently varied new waterline node physically moves, with no
    # changes to old eta controls. It is not a silently interpolated height.
    eta = np.zeros(len(chart.free_nodes))
    column = np.searchsorted(chart.free_nodes, chart.high_contact_nodes[len(chart.high_contact_nodes)//2])
    eta[column] = 2e-6
    changed = chart.evaluate(domain.rest.position, domain.rest.orientation, eta)
    assert changed.points[chart.free_nodes[column], 2] == 2e-6
    assert np.linalg.norm(changed.points-initial.points) > 2e-6
    np.testing.assert_allclose(changed.points[:len(domain.contact.reference_points)],
                               initial.points[:len(domain.contact.reference_points)], atol=0, rtol=0)


def test_finite_pose_nonflat_heights_and_all_six_world_twists_match_geometry_fd(case):
    domain, chart = case
    position = domain.rest.position+np.array([.003, -.002, .005])
    orientation = moved(domain.rest.orientation, [.008, -.006, .004], 1.)
    xy = chart.reference_points[chart.free_nodes, :2]
    eta = 8e-5*np.cos(1.3*xy[:, 0])+3e-5*np.sin(2*xy[:, 1])
    own = np.searchsorted(chart.free_nodes, chart.high_contact_nodes)
    eta[own] += 2e-6*np.sin(np.arange(len(own)))
    geometry = chart.evaluate(position, orientation, eta)
    np.testing.assert_array_equal(geometry.points[chart.free_nodes, 2], eta)
    eta_rate = 3e-4*np.cos(np.arange(len(eta))*.9)
    directions = [(np.eye(6)[i], np.zeros(len(eta))) for i in range(6)]
    directions.append((np.array([.07, -.04, .03, .02, -.03, .015]), eta_rate))
    for twist, height_rate in directions:
        expected = geometry.velocity(height_rate, twist)
        errors = []
        for step in (1e-5, 5e-6):
            plus = chart.evaluate(position+step*twist[:3], moved(orientation, twist[3:], step), eta+step*height_rate)
            minus = chart.evaluate(position-step*twist[:3], moved(orientation, twist[3:], -step), eta-step*height_rate)
            difference = (plus.points-minus.points)/(2*step)
            errors.append(np.max(abs(difference-expected)))
        assert errors[-1] < 8e-9, (twist, errors)
        assert errors[-1] < .7*errors[0]+8e-11, (twist, errors)
    # Independent off-grid physical facet test, not the chart's guard points.
    barycentric, _ = triangle_samples(np.eye(3)[None], 6)
    values, _ = basis(barycentric[0])
    sample = np.einsum('qa,fai->fqi', values, geometry.points[chart.faces[chart.linear.body_faces]])
    world = (chart.linear.physical_vertices-chart.linear.cg)@rotation_matrix(orientation).T+position
    parent = world[chart.linear.physical_faces[chart.linear.body_source_faces]]
    edge = parent[:, 1:]-parent[:, :1]
    normal = np.cross(edge[:, 0], edge[:, 1])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    assert np.max(abs(np.einsum('fqi,fi->fq', sample-parent[:, :1], normal))) < 3e-15
    uv = np.einsum('fqi,fji->fqj', sample-parent[:, :1], edge)
    uv = np.einsum('fqj,fjk->fqk', uv, np.linalg.inv(edge@edge.transpose(0, 2, 1)))
    assert np.min(uv) > -1e-12
    assert np.max(uv.sum(axis=2)) < 1+1e-12
    twist = directions[-1][0]
    velocity = geometry.velocity(eta_rate, twist)
    sampled_velocity = np.einsum('qa,fai->fqi', values,
        velocity[chart.faces[chart.linear.body_faces]])
    rigid = twist[:3]+np.cross(twist[3:], sample-position)
    assert np.max(abs(np.einsum('fqi,fi->fq', sampled_velocity-rigid, normal))) < 2e-14


def test_curved_waterline_trace_is_shared_by_hull_and_free_surface(case):
    domain, chart = case
    node = chart.high_contact_nodes[len(chart.high_contact_nodes)//2]
    eta = np.zeros(len(chart.free_nodes)); eta[np.searchsorted(chart.free_nodes, node)] = 2e-6
    geometry = chart.evaluate(domain.rest.position, domain.rest.orientation, eta)
    support = chart.lift.indices[chart.lift.indptr[node]:chart.lift.indptr[node+1]]
    traces = []
    for family in (chart.linear.body_faces, chart.linear.free_faces):
        face_id = next(i for i in family if np.all(np.isin(support, chart.linear.faces[i])))
        face = chart.linear.faces[face_id]
        t = np.linspace(0, 1, 19)
        bary = np.zeros((len(t), 3))
        bary[:, np.flatnonzero(face == support[0])[0]] = 1-t
        bary[:, np.flatnonzero(face == support[1])[0]] = t
        values, _ = basis(bary)
        traces.append(values@geometry.points[chart.faces[face_id]])
    np.testing.assert_allclose(traces[0], traces[1], atol=2e-15, rtol=0)
    assert np.max(abs(traces[0][:, 2])) > 1e-6
    assert traces[0][0, 2] == traces[0][-1, 2] == 0


def test_an_excessive_independent_contact_mode_stops_at_the_sampled_chart_limit(case):
    domain, chart = case
    eta = np.zeros(len(chart.free_nodes))
    node = chart.high_contact_nodes[len(chart.high_contact_nodes)//2]
    eta[np.searchsorted(chart.free_nodes, node)] = .05
    with pytest.raises(ContactTopologyError, match='sampled P4 wet patch'):
        chart.evaluate(domain.rest.position, domain.rest.orientation, eta)
