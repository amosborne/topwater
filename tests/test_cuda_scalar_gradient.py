"""Shared sample owner / contracted scalar gate, CUDA optional."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

pytest.importorskip('cupy')
from .experiments.free_surface import cuda_operator as gpu
from .experiments.free_surface.cuda_shape_gradient import vertex_gradient as raw_gradient
from .experiments.free_surface.cuda_scalar_gradient import vertex_gradient as scalar_gradient
from .experiments.free_surface.cuda_reverse_gradient import vertex_gradient as reverse_gradient
from .experiments.free_surface.cuda_pair_gradient_probe import selected_pair_derivatives
from .experiments.free_surface.cuda_directional import energy_direction_action
from .experiments.free_surface.cuda_kernel_linearization import RawKernelAnchor
from topwater.rankine import TrianglePanels


@pytest.mark.parametrize('order',[3,6])
def test_contracted_gradient_preserves_current_energy(order):
    try:
        if gpu.cp.cuda.runtime.getDeviceCount()==0:
            pytest.skip('CUDA device unavailable')
    except gpu.cp.cuda.runtime.CUDARuntimeError:
        pytest.skip('CUDA runtime unavailable')
    points=np.array([[.1,.2,-.3],[1.3,.1,-.2],[.2,1.1,-.1],[.3,.2,.9]])
    faces=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]],dtype=np.int64)
    phi=np.array([.2,-.4,.6,.1]);rho=1025.
    assembly=gpu.prepare_energy(faces,TrianglePanels(points[faces],1),order=order)
    raw=raw_gradient(points,faces,phi,assembly=assembly,order=order,density=rho,batch_width=1)
    contracted=scalar_gradient(points,faces,phi,assembly=assembly,order=order,density=rho)
    assert_allclose(contracted,raw,rtol=3e-10,atol=3e-10)
    reverse=reverse_gradient(points,faces,phi,assembly=assembly,order=order,density=rho)
    assert_allclose(reverse,raw,rtol=3e-10,atol=3e-10)
    direction=np.random.default_rng(53).normal(size=points.shape)
    action,_=energy_direction_action(points,faces,direction,phi,assembly=assembly)
    assert_allclose(np.sum(contracted*direction),.5*rho*phi@action,rtol=3e-10,atol=3e-10)
    assert_allclose(contracted.sum(axis=0),0.,atol=3e-10)
    assert_allclose(np.cross(points,contracted).sum(axis=0),0.,atol=3e-10)
    assert_allclose(np.sum(reverse*direction),.5*rho*phi@action,rtol=3e-10,atol=3e-10)
    forward_pairs,reverse_pairs=selected_pair_derivatives(points,faces,phi,
        np.array([[0,0],[0,1],[1,0],[1,1]]),assembly=assembly,density=rho,order=order)
    assert_allclose(reverse_pairs,forward_pairs,rtol=3e-10,atol=3e-10)
    # The existing one-anchor raw adapter remains API-compatible after refactor.
    anchor=RawKernelAnchor(points,faces,chart_token='same-tetra',order=order)
    owned=anchor.evaluate(points,faces,chart_token='same-tetra')
    assert_allclose(owned.energy,assembly.energy,rtol=3e-12,atol=3e-13)
    with pytest.raises(ValueError,match='matching quadrature'):
        scalar_gradient(points+1e-8,faces,phi,assembly=assembly,order=order)
