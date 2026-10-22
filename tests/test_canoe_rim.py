"""Physical rim immersion uses the represented surface, including edge interiors."""
import numpy as np
import pytest

from topwater.frames import rotation_matrix
from topwater.hydrostatics import FloodingError
from tests.test_canoe_contact import canoe
from tests.experiments.free_surface.canoe_contact import RimCoverageError


def test_rim_clearance_on_flat_and_tilted_surfaces(canoe):
    _, rest, _, contact = canoe
    points = contact.reference_points.copy()
    rim = (contact.rim_vertices-contact.cg)@rotation_matrix(rest.orientation).T+rest.position
    flat = contact.rim_clearance(points, rest.position, rest.orientation)
    np.testing.assert_allclose(flat, rim[:, 2].min(), atol=2e-15)
    slope = np.array([.01, -.02])
    points[contact.free_nodes, 2] = points[contact.free_nodes, :2]@slope+.003
    expected = np.min(rim[:, 2]-rim[:, :2]@slope-.003)
    np.testing.assert_allclose(contact.rim_clearance(points, rest.position, rest.orientation),
                               expected, atol=2e-15)
    # Query the same defined graph with a lower physical rim. This isolates the
    # physical event from the separate fixed-contact-chart numerical limit.
    with pytest.raises(FloodingError):
        contact.rim_clearance(points, rest.position-np.array([0., 0., .3]), rest.orientation)
    with pytest.raises(RimCoverageError):
        contact.rim_clearance(points, rest.position+np.array([10., 0., 0.]), rest.orientation)


def test_immersion_between_rim_vertices_is_detected(canoe):
    _, rest, _, contact = canoe
    points = contact.reference_points.copy()
    rim = (contact.rim_vertices-contact.cg)@rotation_matrix(rest.orientation).T+rest.position
    faces = contact.faces[contact.free_faces]
    triangles = points[faces, :2]
    inverse = np.linalg.inv(triangles[:, 1:]-triangles[:, :1])

    def interpolate_nodes(xy):
        """Dense sample evaluator used only to construct an interior peak."""
        values = np.zeros((len(xy), len(points)))
        for i, location in enumerate(xy):
            tail = np.einsum('fi,fij->fj', location-triangles[:, 0], inverse)
            bary = np.column_stack((1-tail.sum(axis=1), tail))
            match = np.flatnonzero(np.all(bary >= -1e-12, axis=1))
            assert len(match)
            face = match[0]
            values[i, faces[face]] = bary[face]
        return values

    endpoint_weights = interpolate_nodes(rim[:, :2])
    fractions = np.linspace(0, 1, 31)
    samples = (rim[:, None]+fractions[None, :, None]*(
        np.roll(rim, -1, axis=0)-rim)[:, None]).reshape(-1, 3)
    sample_weights = interpolate_nodes(samples[:, :2])
    endpoint_max = endpoint_weights.max(axis=0)
    interior_max = sample_weights.max(axis=0)
    node = np.argmax(interior_max-endpoint_max)
    assert interior_max[node] > endpoint_max[node]+.01
    # Raise one existing surface node. The threshold lies between its maximum
    # influence at physical rim vertices and its larger influence inside an edge.
    height = rim[0, 2]/((interior_max[node]+endpoint_max[node])/2)
    points[node, 2] = height
    assert np.min(rim[:, 2]-endpoint_weights@points[:, 2]) > 0
    assert np.min(samples[:, 2]-sample_weights@points[:, 2]) < 0
    with pytest.raises(FloodingError):
        contact.rim_clearance(points, rest.position, rest.orientation)
