"""Small GPU-optional consistency falsifiers, not a boat accuracy certificate."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

pytest.importorskip('cupy')
from .experiments.free_surface import cuda_operator as gpu
from .experiments.free_surface.cuda_kernel_linearization import RawKernelAnchor
from .experiments.free_surface.cuda_shape_gradient import vertex_gradient
from topwater.rankine import TrianglePanels


@pytest.fixture
def tetra():
    try:
        if gpu.cp.cuda.runtime.getDeviceCount()==0:
            pytest.skip('CUDA device unavailable')
    except gpu.cp.cuda.runtime.CUDARuntimeError:
        pytest.skip('CUDA runtime unavailable')
    points=np.array([[.1,.2,-.3],[1.3,.1,-.2],[.2,1.1,-.1],[.3,.2,.9]])
    faces=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]],dtype=np.int64)
    return points,faces


def test_anchor_matches_fresh_and_translation(tetra):
    points,faces=tetra
    phi=np.array([.2,-.4,.6,.1]);rho=1025.
    fresh=gpu.prepare_energy(faces,TrianglePanels(points[faces],1),order=3)
    cache=RawKernelAnchor(points,faces,chart_token=('tetra',1),order=3)
    current=cache.evaluate(points,faces,chart_token=('tetra',1))
    assert cache.pair_bytes==4*4*4*19*8
    assert_allclose(current.energy,fresh.energy,rtol=2e-12,atol=2e-13)
    exact_gradient=vertex_gradient(points,faces,phi,assembly=fresh,density=rho,order=3)
    gradient=current.vertex_gradient(phi,density=rho)
    assert_allclose(gradient,exact_gradient,rtol=3e-10,atol=3e-10)
    translated=cache.evaluate(points+[.31,-.27,.18],faces,chart_token=('tetra',1))
    assert_allclose(translated.energy,current.energy,rtol=3e-11,atol=3e-12)
    assert_allclose(translated.vertex_gradient(phi,density=rho),gradient,rtol=3e-10,atol=3e-10)
    assert_allclose(gradient.sum(axis=0),0.,atol=3e-10)
    assert_allclose(current.mixed_trace.get()@np.ones(4),fresh.mixed_trace.get()@np.ones(4),atol=2e-14)


def test_current_gradient_differentiates_same_approximate_energy(tetra):
    points,faces=tetra
    cache=RawKernelAnchor(points,faces,chart_token='fixed-tetra',order=3)
    rng=np.random.default_rng(35)
    current_points=points+.015*rng.normal(size=points.shape)
    phi=rng.normal(size=len(points));direction=rng.normal(size=points.shape);rho=1025.
    current=cache.evaluate(current_points,faces,chart_token='fixed-tetra')
    actual=np.sum(current.vertex_gradient(phi,density=rho)*direction)
    for h in (1e-5,5e-6):
        plus=cache.evaluate(current_points+h*direction,faces,chart_token='fixed-tetra')
        minus=cache.evaluate(current_points-h*direction,faces,chart_token='fixed-tetra')
        expected=.5*rho*(phi@(plus.energy-minus.energy)@phi)/(2*h)
        assert_allclose(actual,expected,rtol=2e-7,atol=2e-6)
    # Changing phi must recompute q, not reuse an anchor/current previous q.
    other=current.vertex_gradient(phi*1.7,density=rho)
    assert_allclose(other,current.vertex_gradient(phi,density=rho)*1.7**2,rtol=2e-12,atol=2e-10)


def test_ownership_and_allocation_rejections(tetra):
    points,faces=tetra
    with pytest.raises(MemoryError):
        RawKernelAnchor(points,faces,chart_token='owned',max_pair_bytes=0)
    with pytest.raises(ValueError,match='chart token'):
        RawKernelAnchor(points,faces,chart_token=None)
    cache=RawKernelAnchor(points,faces,chart_token='owned')
    with pytest.raises(ValueError,match='semantic chart'):
        cache.evaluate(points,faces,chart_token='other')
    with pytest.raises(ValueError,match='semantic chart'):
        cache.evaluate(points,faces[:,::-1],chart_token='owned')
    current=cache.evaluate(points,faces,chart_token='owned')
    with pytest.raises(ValueError,match='exact current'):
        current.auxiliary(points+1e-8,faces,np.ones(4))
    with pytest.raises(ValueError,match='exact current'):
        current.auxiliary(points,faces,np.ones(4),order=6)
