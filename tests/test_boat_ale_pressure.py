"""Exact translating harmonic flow with tangentially sliding hull mesh nodes."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

from .experiments.free_surface.boat_ale_pressure import differentiate
from .experiments.free_surface.body_pressure import evaluate_body_pressure


def test_sliding_chart_recovers_exact_eulerian_pressure():
    # Oblique patch, four triangles and one interior node whose material label
    # slides. Uniform flow matches the translating wall's normal velocity.
    a, b = np.array([1., 0., .4]), np.array([0., 1., -.3])
    points = np.array([np.zeros(3), a, a+b, b, .5*(a+b)])
    faces = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
    base = dict(points=points, faces=faces, body_faces=np.arange(4),
        free_faces=np.array([], dtype=int), free_nodes=np.array([], dtype=int),
        body_source_group=np.zeros(4, dtype=int), free_source_parent=np.array([], dtype=int),
        free_parent_indices=np.empty((0, 3), dtype=int), node_chart_key=np.array('exact-patch'))
    translation = np.array([.3, -.2, .1])
    velocity = np.broadcast_to(translation, points.shape).copy()
    velocity[-1] += 2*a-b
    h, slope = .001, .7
    plus, minus = dict(base, points=points+h*velocity), dict(base, points=points-h*velocity)
    rate, mesh_velocity, diagnostics = differentiate(base, plus, minus,
        plus['points']@translation+h*slope, minus['points']@translation-h*slope,
        h, np.zeros(3), np.r_[translation, np.zeros(3)])
    result = evaluate_body_pressure(points, faces, base['body_faces'], np.zeros(3),
        10., np.eye(3), np.r_[translation, np.zeros(3)], points@translation,
        rate, mesh_velocity, np.zeros((len(points), 6)), density=1000., gravity=0.)
    area_vector = np.cross(a, b)
    pressure = -1000*(slope+.5*translation@translation)
    assert_allclose(result.known_pressure_load[:3], pressure*area_vector, atol=1e-9)
    assert_allclose(result.known_pressure_load[3:],
                    np.cross(.5*(a+b), pressure*area_vector), atol=1e-9)
    assert diagnostics['maximum_body_normal_velocity_error_m_per_s'] < 1e-12
    with pytest.raises(ValueError, match='same chart'):
        differentiate(base, dict(plus, node_chart_key=np.array('different')), minus,
                      np.zeros(5), np.zeros(5), h, np.zeros(3), np.zeros(6))
    with pytest.raises(ValueError, match='same chart'):
        differentiate(base, dict(plus, faces=faces[::-1]), minus,
                      np.zeros(5), np.zeros(5), h, np.zeros(3), np.zeros(6))
