"""Actual canoe waterline geometry and its six-axis ALE velocity."""
import json

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from topwater import HullMesh, SimulationInputs
from topwater.frames import rotation_matrix
from topwater.hydrostatics import equilibrate_heave
from tests.reference_hull import input_data
from tests.experiments.free_surface.canoe_contact import CanoeContact, ContactTopologyError
from tests.experiments.free_surface.pocket import grid, boundary


@pytest.fixture(scope='module', params=[2, 3], ids=['material-2', 'material-3'])
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
    geometry = contact.evaluate(position, orientation, eta)
    HullMesh(geometry.points, contact.faces)
    world = (contact.material_vertices-contact.cg)@rotation_matrix(orientation).T+position
    a, b = world[contact.contact_edges[:, 0]], world[contact.contact_edges[:, 1]]
    p = geometry.points[contact.contact_nodes]
    np.testing.assert_allclose(np.cross(p-a, b-a), 0., atol=1e-17)
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
    expected = np.zeros_like(actual)
    expected[contact.free_nodes] = geometry.eta_jacobian[contact.free_nodes]*eta_rate[:, None]
    np.testing.assert_allclose(actual, expected, rtol=3e-7, atol=3e-10)


def test_contact_topology_event_is_explicit(canoe):
    _, rest, _, contact = canoe
    eta = np.zeros(len(contact.free_nodes))
    with pytest.raises(ContactTopologyError, match='material vertex'):
        contact.evaluate(rest.position+np.array([0., 0., .1]), rest.orientation, eta)
    eta[0] = np.nan
    with pytest.raises(ValueError, match='finite'):
        contact.evaluate(rest.position, rest.orientation, eta)
