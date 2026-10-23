"""Guarded source values against CPU and unsupported-geometry failure cases."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

pytest.importorskip('cupy')
from topwater.rankine import TrianglePanels, panel_influence, QuadratureError
from topwater.numerics import LaplaceControls
from tests.experiments.free_surface.cuda_source_potential import paired_potential
from tests.experiments.free_surface.cuda_graph import cp


def panels():
    triangle = np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.]])
    return TrianglePanels(np.array([triangle, triangle+[2., .3, -.2]]), -1)


def test_nonself_graph_targets_and_reflection_match_cpu():
    source = panels()
    points = np.array([[.2,.3,.2], [2.2,.6,.1], [-.4,.6,.3], [4.,-2.,.2]])
    reflected = points.copy(); reflected[:, 2] = 1.-points[:, 2]
    expected = panel_influence(points, source).single_layer-panel_influence(reflected, source).single_layer
    result = paired_potential(points, source, .5)
    assert_allclose(cp.asnumpy(result.values), expected, atol=3e-15, rtol=3e-13)
    assert result.diagnostics['gpu_flagged_pairs'] == 0


@pytest.mark.parametrize('point', [[0.,0.,0.], [.5,0.,0.], [1/3,1/3,0.], [0.,0.,1.]])
def test_vertices_edges_unowned_centroids_and_reflected_vertices_still_reject(point):
    with pytest.raises(ValueError):
        paired_potential(np.array([point]), panels(), .5)


def test_exact_image_plane_cancellation_has_zero_value_after_cpu_geometry_check():
    result = paired_potential(np.array([[.2,.3,.5]]), panels(), .5,
        LaplaceControls(quadrature_atol=1e-20, quadrature_rtol=1e-12))
    assert_allclose(cp.asnumpy(result.values), 0., atol=0., rtol=0.)
    assert result.diagnostics['paired_accuracy_flags'] > 0


def test_unresolved_near_cancellation_fails_instead_of_returning_unqualified_digits():
    with pytest.raises(QuadratureError, match='paired source cancellation'):
        paired_potential(np.array([[.2,.3,.5+1e-12]]), panels(), .5,
            LaplaceControls(quadrature_atol=1e-25, quadrature_rtol=1e-12))


def test_operator_memory_preflight_runs_before_device_allocation():
    with pytest.raises(MemoryError):
        paired_potential(np.zeros((5,3)), panels(), .5, LaplaceControls(max_operator_bytes=1))


def test_constant_source_field_is_additive_under_panel_refinement_and_renumbering():
    from tests.experiments.free_surface.boat_spatial_work import subdivide
    source = panels()
    points = np.array([[.2,.3,.2], [2.2,.6,.1], [-.4,.6,.3]])
    coarse = cp.asnumpy(paired_potential(points,source,.5).values)
    reordered = TrianglePanels(source.triangles[::-1], 1)
    assert_allclose(cp.asnumpy(paired_potential(points,reordered,.5).values)[:,::-1],coarse,atol=2e-15,rtol=2e-13)
    refined = TrianglePanels(subdivide(source.triangles), -1)
    assert_allclose(cp.asnumpy(paired_potential(points,refined,.5).values).sum(axis=1),
                    coarse.sum(axis=1),atol=3e-15,rtol=3e-13)


def test_far_accuracy_flags_recover_through_cpu_guards():
    source = panels()
    points = np.array([[1e16,1e16,.2]])
    reflected = points.copy(); reflected[:,2] = 1.-points[:,2]
    expected = panel_influence(points,source).single_layer-panel_influence(reflected,source).single_layer
    result = paired_potential(points,source,.5)
    assert result.diagnostics['individual_accuracy_flags']>0
    assert_allclose(cp.asnumpy(result.values),expected,atol=1e-25,rtol=0.)
