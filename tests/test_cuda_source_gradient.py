"""Bulk image-chain rule, actual scalar derivative and rare accuracy recovery."""
import numpy as np
from types import SimpleNamespace
from numpy.testing import assert_allclose
import pytest

pytest.importorskip('cupy')
from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels, QuadratureError
from tests.experiments.free_surface.panel_gradient_reference import checked_influence
from tests.experiments.free_surface.cuda_source_potential import paired_gradient,paired_velocity,paired_potential
from tests.experiments.free_surface.cuda_graph import cp


def source():
    return TrianglePanels(np.array([[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]],
        [[2.,.2,-.5],[3.1,.3,-.4],[2.2,1.2,.2]]]),-1)


def points():
    return np.array([[.2,.3,.7],[2.3,.6,.45],[-.4,.7,.2],[.3,.2,-.3],[3.,-2.,.4]])


def test_bulk_gradient_matches_cpu_on_both_sides_of_sources_and_image_chain():
    panels=source();p=points();mirror=p.copy();mirror[:,2]=1.6-p[:,2]
    expected=checked_influence(p,panels).gradient_single_layer \
        -checked_influence(mirror,panels).gradient_single_layer*[1.,1.,-1.]
    result=paired_gradient(p,panels,.8)
    assert_allclose(cp.asnumpy(result.values),expected,atol=3e-14,rtol=3e-12)
    assert result.diagnostics['gpu_flagged_pairs']==0
    h=1e-4
    for axis in range(3):
        shift=np.eye(3)[axis]*h
        difference=cp.asnumpy(paired_potential(p+shift,panels,.8).values
                            -paired_potential(p-shift,panels,.8).values)/(2*h)
        assert_allclose(cp.asnumpy(result.values[:,:,axis]),difference,atol=2e-9,rtol=3e-7)


def test_density_contraction_is_memory_bounded_and_preserves_all_components():
    panels=source();p=points();density=np.array([1.3,-.4])
    expected=cp.asnumpy(paired_gradient(p,panels,.8).values)
    expected=np.einsum('tnj,n->tj',expected,density)
    result=paired_velocity(p,panels,.8,density,LaplaceControls(max_operator_bytes=100))
    assert result.diagnostics['batches']==3 and result.diagnostics['max_matrix_bytes']<=100
    assert_allclose(cp.asnumpy(result.values),expected,atol=3e-15,rtol=3e-13)


def test_graph_bulk_velocity_uses_same_guarded_evaluator():
    from tests.experiments.free_surface.boat_graph_fields import source_velocity_bulk
    panels=source();p=points();density=np.array([1.3,-.4])
    context=SimpleNamespace(panels=panels,image_height=.8)
    reflected=p.copy();reflected[:,2]=1.6-p[:,2]
    gradient=checked_influence(p,panels).gradient_single_layer \
        -checked_influence(reflected,panels).gradient_single_layer*[1.,1.,-1.]
    assert_allclose(source_velocity_bulk(context,density,p),
        np.einsum('tnj,n->tj',gradient,density),atol=3e-14,rtol=3e-12)


@pytest.mark.parametrize('point',[[0.,0.,0.],[.5,0.,0.],[1/3,1/3,0.],[0.,0.,1.6]])
def test_bulk_gradient_never_invents_edge_vertex_or_unowned_trace_limits(point):
    with pytest.raises(ValueError):
        paired_gradient(np.array([point]),source(),.8)


def test_exact_image_plane_has_zero_tangential_velocity_and_doubled_vertical():
    panels=source();p=np.array([[.2,.3,.8]])
    result=paired_gradient(p,panels,.8,LaplaceControls(quadrature_atol=1e-25,quadrature_rtol=1e-12))
    expected=checked_influence(p,panels).gradient_single_layer
    assert_allclose(cp.asnumpy(result.values)[:,:,:2],0.,atol=0.,rtol=0.)
    assert_allclose(cp.asnumpy(result.values)[:,:,2],2*expected[:,:,2],atol=3e-14,rtol=3e-12)
    assert result.diagnostics['paired_accuracy_flags']>0


def test_unresolved_close_paired_gradient_does_not_claim_impossible_accuracy():
    with pytest.raises(QuadratureError,match='paired gradient'):
        paired_gradient(np.array([[.2,.3,.8+1e-12]]),source(),.8,
            LaplaceControls(quadrature_atol=1e-25,quadrature_rtol=1e-16))


def test_far_recovery_obeys_area_weighted_point_source_limit():
    panels=TrianglePanels(source().triangles[:1],-1)
    p=np.array([[1e8,1e8,.2]])
    result=paired_gradient(p,panels,.8,LaplaceControls(quadrature_atol=1e-35,quadrature_rtol=1e-10))
    delta=p[0]-panels.centroids[0];radius=np.linalg.norm(delta)
    coefficient=panels.areas[0]*(.8-panels.centroids[0,2])/(4*np.pi)
    expected=-6*(.8-p[0,2])*coefficient*delta/radius**5
    expected[2]=-2*coefficient/radius**3
    assert result.diagnostics['gpu_flagged_pairs']>0
    # O((panel diameter/radius)^2) finite-panel correction is negligible here.
    assert_allclose(cp.asnumpy(result.values)[0,0],expected,atol=1e-42,rtol=1e-7)


def test_far_vertical_cancellation_does_not_hide_constituent_roundoff():
    panels=TrianglePanels(source().triangles[:1],-1)
    # Dipole vertical velocity vanishes near z/r=1/sqrt(3). Its two
    # constituent terms remain finite, so relative output error is misleading.
    p=panels.centroids.copy()
    p[0,0]+=1e8
    p[0,2]=.8+1e8/np.sqrt(2)
    with pytest.raises(QuadratureError,match='paired gradient'):
        paired_gradient(p,panels,.8,
            LaplaceControls(quadrature_atol=1e-40,quadrature_rtol=1e-12))
