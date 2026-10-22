"""Affine physical vertical velocity survives piecewise graph slopes."""
import numpy as np
from numpy.testing import assert_allclose

from tests.experiments.free_surface.boat_contact_trace import compare_traces


def test_auxiliary_affine_velocity_on_piecewise_sloping_graph():
    points = np.array([[0., 0., 0.], [1., 0., .2], [0., 1., -.1], [1., 1., .4]])
    faces = np.array([[0, 1, 2], [1, 3, 2]])
    velocity = np.array([.3, -.2, .4])
    normals = np.cross(points[faces[:, 1]]-points[faces[:, 0]], points[faces[:, 2]]-points[faces[:, 0]])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    result = compare_traces(points, faces, np.arange(4), points@velocity+.7,
        np.ones(4)*.4, normals@velocity, np.array([[0, 1]]), gravity=0.)
    rows = [row for row in result['bins'] if row['quadrature_projected_area_m2']]
    assert_allclose(sum(row['quadrature_projected_area_m2'] for row in rows), 1.)
    for row in rows:
        assert_allclose(row['auxiliary']['maximum_absolute_vertical_m_per_s'], .4, atol=2e-16)
        assert_allclose(row['auxiliary']['vertical_square_integral_m4_per_s2'],
                        row['quadrature_projected_area_m2']*.4**2, atol=2e-16)
