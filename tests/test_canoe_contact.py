"""Actual canoe waterline geometry and its six-axis ALE velocity."""
import json

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from topwater import HullMesh, SimulationInputs
from topwater.frames import rotation_matrix
from dataclasses import replace
from topwater.hydrostatics import equilibrate_heave, hydrostatics
from tests.reference_hull import input_data
from tests.experiments.free_surface.canoe_contact import CanoeContact, ContactTopologyError
from tests.experiments.free_surface.pocket import grid, boundary


@pytest.fixture(scope='module', params=[2, 3], ids=['wet-subdivision-2', 'wet-subdivision-3'])
def canoe(request):
    data = SimulationInputs.from_json(json.dumps(input_data(2)))
    rest, wet = equilibrate_heave(data.body, data.initial_state)
    points, tets, near, _ = grid(12, 0., 0)
    near &= np.all(points[tets, 2] >= -1., axis=1)
    faces = boundary(points, tets[near])
    top = np.all(points[faces, 2] == 0, axis=1)
    contact = CanoeContact.make(data.body, rest, points[faces[~top]],
                               points[np.unique(faces[top])], subdivisions=request.param)
    return data, rest, wet, contact


def test_rest_matches_physical_canoe_and_is_closed(canoe):
    _, rest, wet, contact = canoe
    geometry = contact.evaluate(rest.position, rest.orientation,
                               np.zeros(len(contact.free_nodes)))
    np.testing.assert_allclose(geometry.points, contact.reference_points, atol=5e-16)
    HullMesh(geometry.points, contact.faces)
    body = geometry.points[contact.faces[contact.body_faces]]
    area = np.linalg.norm(np.cross(body[:, 1]-body[:, 0], body[:, 2]-body[:, 0]), axis=1).sum()/2
    expected = wet.wetted_triangles
    expected_area = np.linalg.norm(np.cross(expected[:, 1]-expected[:, 0],
                                          expected[:, 2]-expected[:, 0]), axis=1).sum()/2
    np.testing.assert_allclose(area, expected_area, rtol=2e-14)
    volume = -np.einsum('ij,ij->i', body[:, 0],
                        np.cross(body[:, 1], body[:, 2])).sum()/6
    np.testing.assert_allclose(volume, wet.displaced_volume, rtol=2e-14)
    assert set(contact.contact_nodes) == set(np.intersect1d(
        contact.free_nodes, np.unique(contact.faces[contact.body_faces])))


def test_moved_canoe_nonflat_surface_and_analytic_ale_derivatives(canoe):
    _, rest, _, contact = canoe
    position = rest.position+np.array([.003, -.002, .0002])
    rotation = Rotation.from_rotvec([.001, -.0007, .0013])
    orientation = rotation.as_quat(scalar_first=True)
    xy = contact.reference_points[contact.free_nodes, :2]
    eta = .0003*np.sin(.6*xy[:, 0]+.3*xy[:, 1]+.2)
    fine_contact = np.setdiff1d(contact.contact_nodes, contact.coarse_nodes)
    eta[np.searchsorted(contact.free_nodes, fine_contact)] += .0002*np.sin(np.arange(len(fine_contact)))
    coarse_eta = eta[np.searchsorted(contact.free_nodes, contact.coarse_nodes)]
    difference = eta[np.searchsorted(contact.free_nodes, contact.contact_nodes)]-contact.contact_blend@coarse_eta
    assert np.max(np.abs(difference)) > .0001
    geometry = contact.evaluate(position, orientation, eta)
    HullMesh(geometry.points, contact.faces)
    world = (contact.physical_vertices-contact.cg)@rotation_matrix(orientation).T+position
    a, b = world[contact.contact_edges[:, 0]], world[contact.contact_edges[:, 1]]
    p = geometry.points[contact.coarse_nodes]
    np.testing.assert_allclose(np.cross(p-a, b-a), 0., atol=1e-16)
    np.testing.assert_array_equal(geometry.points[contact.free_nodes, 2], eta)
    h = 1e-6
    for axis in range(6):
        direction = np.eye(6)[axis]
        plus_q = (Rotation.from_rotvec(h*direction[3:])*rotation).as_quat(scalar_first=True)
        minus_q = (Rotation.from_rotvec(-h*direction[3:])*rotation).as_quat(scalar_first=True)
        plus = contact.evaluate(position+h*direction[:3], plus_q, eta)
        minus = contact.evaluate(position-h*direction[:3], minus_q, eta)
        np.testing.assert_allclose((plus.points-minus.points)/(2*h),
            geometry.body_jacobian[:, :, axis], rtol=3e-7, atol=3e-10)
    eta_rate = np.cos(xy[:, 0])
    actual = (contact.evaluate(position, orientation, eta+h*eta_rate).points
              -contact.evaluate(position, orientation, eta-h*eta_rate).points)/(2*h)
    expected = np.column_stack([matrix@eta_rate for matrix in geometry.height_jacobians])
    np.testing.assert_allclose(actual, expected, rtol=3e-7, atol=3e-10)
    np.testing.assert_array_equal(geometry.height_jacobians[2][contact.free_nodes].toarray(),
                                  np.eye(len(contact.free_nodes)))
    np.testing.assert_array_equal(geometry.body_jacobian[contact.free_nodes, 2], 0.)
    # Wet-node redistribution is actual ALE motion in the physical facet.
    twist = np.array([.2, -.1, .3, .1, -.2, .3])
    velocity = geometry.velocity(contact.free_nodes, eta_rate, twist)
    nodes = contact.faces[contact.body_faces]
    parent = world[contact.physical_faces[contact.body_source_faces]]
    normal = np.cross(parent[:, 1]-parent[:, 0], parent[:, 2]-parent[:, 0])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    rigid = twist[:3]+np.cross(twist[3:], geometry.points[nodes]-position)
    residual = np.einsum('fci,fi->fc', velocity[nodes]-rigid, normal)
    np.testing.assert_allclose(residual, 0., atol=2e-13)


@pytest.mark.parametrize('heave,rotvec', [(.01, [0, 0, 0]), (-.01, [0, 0, 0]),
    (0, [np.pi/180, 0, 0]), (0, [-np.pi/180, 0, 0]),
    (0, [0, np.pi/180, 0]), (0, [0, -np.pi/180, 0])])
def test_finite_motion_preserves_physical_wet_hull(canoe, heave, rotvec):
    data, rest, _, contact = canoe
    state = replace(rest, position=rest.position+np.array([0., 0., heave]),
        orientation=Rotation.from_rotvec(rotvec).as_quat(scalar_first=True))
    geometry = contact.evaluate(state.position, state.orientation,
                                np.zeros(len(contact.free_nodes)))
    HullMesh(geometry.points, contact.faces)
    wet = hydrostatics(data.body, state)
    body = geometry.points[contact.faces[contact.body_faces]]
    area = np.linalg.norm(np.cross(body[:, 1]-body[:, 0], body[:, 2]-body[:, 0]), axis=1).sum()/2
    volume = -np.einsum('ij,ij->i', body[:, 0], np.cross(body[:, 1], body[:, 2])).sum()/6
    np.testing.assert_allclose(area, wet.wetted_area, rtol=2e-13)
    np.testing.assert_allclose(volume, wet.displaced_volume, rtol=2e-13)
    world = (contact.physical_vertices-contact.cg)@rotation_matrix(state.orientation).T+state.position
    parent = world[contact.physical_faces[contact.body_source_faces]]
    edges = parent[:, 1:]-parent[:, :1]
    normal = np.cross(edges[:, 0], edges[:, 1])
    signed = np.einsum('fi,fi->f', normal,
        np.cross(body[:, 1]-body[:, 0], body[:, 2]-body[:, 0]))
    assert np.all(signed < 0)
    relative = body-parent[:, :1]
    np.testing.assert_allclose(np.einsum('fci,fi->fc', relative, normal), 0., atol=2e-16)
    uv = np.einsum('fci,fji->fcj', relative, edges)
    uv = np.einsum('fcj,fjk->fck', uv, np.linalg.inv(np.einsum('fji,fki->fjk', edges, edges)))
    assert np.min(uv) > -1e-12
    assert np.max(uv.sum(axis=2)) < 1+1e-12


def test_contact_topology_event_is_explicit(canoe):
    _, rest, _, contact = canoe
    eta = np.zeros(len(contact.free_nodes))
    with pytest.raises(ContactTopologyError, match='physical vertex'):
        contact.evaluate(rest.position+np.array([0., 0., .1]), rest.orientation, eta)
    eta[0] = np.nan
    with pytest.raises(ValueError, match='finite'):
        contact.evaluate(rest.position, rest.orientation, eta)
