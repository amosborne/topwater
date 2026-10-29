"""Actual native CPU/GPU rasterization, independent of video encoding."""
import os
import numpy as np
import pytest
from scipy.ndimage import maximum_filter,minimum_filter

from topwater.render.backend import DrawingBackend


def test_auto_falls_back_to_native_cpu_and_gpu_request_never_emulates(monkeypatch):
    from topwater.render.gpu import backend as gpu
    def unavailable():raise RuntimeError('no physical GPU')
    monkeypatch.setattr(gpu,'select_adapter',unavailable)
    cpu=DrawingBackend((16,16),[-1.,-1.],[1.,1.])
    try:
        assert cpu.info['backend_type']=='Skia' and cpu.info['adapter_type']=='CPU'
    finally:cpu.close()
    with pytest.raises(RuntimeError,match='Hardware WebGPU rendering is unavailable'):
        DrawingBackend((16,16),[-1.,-1.],[1.,1.],device='gpu')


def draw_probe(device,size):
    backend=DrawingBackend(size,[-1.,-1.],[1.,1.],device=device)
    try:
        vertices=backend.vertices(backend.fullscreen,[.2,.4,.8,1])
        # Clockwise concave outer boundary; counterclockwise nested hole.
        outer=np.array([[-.8,-.8],[-.8,.8],[0.,.8],[0.,.4],[.8,.4],[.8,-.8]])
        hole=np.array([[-.4,-.4],[.4,-.4],[.4,0.],[-.4,0.]])
        line=np.array([[[-.9,-.9],[.9,.9]]])
        pixels=backend.draw(vertices,np.arange(6).reshape(2,3),[outer,hole],
            [('line',line,[0.,0.,1.,1.],1.8,False)]).copy()
        return pixels,backend.info
    finally:backend.close()


@pytest.mark.parametrize('size',[(642,642),(64,66)])
def test_native_cpu_winding_fill_hole_concavity_and_readback(size):
    pixels,info=draw_probe('cpu',size)
    assert info['backend_type']=='Skia' and info['adapter_type']=='CPU'
    assert pixels.shape==(size[1],size[0],4) and pixels.dtype==np.uint8
    for x,y,fill in ((-.6,.1,True),(.1,-.2,False),(.5,.65,False),(.6,-.3,True)):
        column=int((x+1)/2*size[0]);row=int((1-y)/2*size[1])
        assert bool(pixels[row,column,0]>200)==fill
    assert np.all(pixels[:,:,3]==255)


def test_native_cpu_remaps_indices_above_uint16_without_losing_triangles():
    backend=DrawingBackend((64,64),[-1.,-1.],[1.,1.],device='cpu')
    try:
        vertices=backend.vertices(np.zeros((70000,2)),[0.,0.,1.,1.])
        vertices[[0,65536,69999],:2]=[[-1,-1],[1,-1],[0,1]]
        pixels=backend.draw(vertices,np.array([[0,65536,69999]]),[],[])
        np.testing.assert_array_equal(pixels[32,32],[0,0,255,255])
    finally:backend.close()


@pytest.mark.skipif(os.environ.get('TOPWATER_TEST_GPU')!='1',reason='explicit GPU opt-in required')
@pytest.mark.parametrize('size',[(642,642),(64,66)])
def test_gpu_and_native_cpu_match_scene_interiors_and_preserve_holes(size):
    cpu,_=draw_probe('cpu',size)
    gpu,info=draw_probe('gpu',size)
    assert info['adapter_type'] in ('DiscreteGPU','IntegratedGPU')
    # Skia analytic edge coverage and GPU four-sample coverage differ at edges.
    # Everywhere outside their one-pixel edge footprint, colors must agree.
    uniform=np.all(maximum_filter(cpu,size=(3,3,1))==minimum_filter(cpu,size=(3,3,1)),axis=2)
    error=abs(gpu.astype(int)-cpu.astype(int))
    assert error[uniform].max()<=1
    assert np.all(gpu[:,:,3]==255)
