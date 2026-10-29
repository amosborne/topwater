"""Opt-in hardware tests: TOPWATER_TEST_GPU=1 uv run --with pyopencl ..."""
import os
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
import pytest


pytestmark = pytest.mark.skipif(os.environ.get('TOPWATER_TEST_GPU') != '1',
                                reason='explicit GPU hardware test opt-in required')


@pytest.fixture(scope='module')
def evaluator():
    pytest.importorskip('pyopencl')
    import capytaine as cpt
    from capytaine.green_functions.delhommeau import Delhommeau
    from topwater.solve.gpu import OpenCLPotential
    # Submerged and free-surface panels exercise near, far, edge and self terms.
    vertices = np.array([[-1.,-1.,-.4],[1.,-1.,-.4],[1.,1.,-.4],[-1.,1.,-.4],
                         [-1.,-1.,0.],[1.,-1.,0.],[1.,1.,0.],[-1.,1.,0.]])
    mesh = cpt.Mesh(vertices=vertices, faces=np.array([[0,1,2,3],[4,5,6,7]]),
                    auto_clean=False, auto_check=False)
    green = Delhommeau(tabulation_cache_dir=str(Path.home()/'.cache/topwater/capytaine-3.0.0'))
    gpu = OpenCLPotential(green, mesh, device_name=os.environ.get('TOPWATER_GPU_DEVICE'))
    return green, mesh, gpu


def test_potential_matches_cpu_at_surface_edges_and_far_field(evaluator):
    green, mesh, gpu = evaluator
    points = np.array([[0.,0.,0.],[1.,1.,0.],[1.,0.,0.],[.123,.456,0.],
                       [20.,20.,0.],[1e3,0.,0.]])
    ks = np.array([1e-8,.025,.5,2.,50.,235.])
    expected = np.array([green.evaluate(points,mesh,wavenumber=float(k))[0] for k in ks])
    actual = gpu.evaluate(points,ks)
    assert_allclose(actual,expected,rtol=1e-10,atol=2e-13)
    # Reuse and invalidation of resident targets must preserve the potential.
    assert_allclose(gpu.evaluate(points,ks),actual,rtol=0,atol=0)
    assert_allclose(gpu.evaluate(points[:2],ks),actual[:,:2],rtol=0,atol=0)


def test_incident_gpu_preserves_signed_branches_phase_and_coherence(evaluator):
    from types import SimpleNamespace
    from topwater.solve.surface import incident_spectrum
    random=np.random.default_rng(21)
    C=random.normal(size=(4,2,5))+1j*random.normal(size=(4,2,5))
    C[:,:,0]=0
    C[3]=0
    train=SimpleNamespace(theta=np.array([0.,.23,.23,2.]),k=np.array([0.,.1,2.,20.,200.]),
                          omega=np.array([0.,1.,2.,3.,4.]),C=C)
    xy=np.array([[-1.,0.],[0.,0.],[.1,.2],[20.,-10.]])
    expected=incident_spectrum(train,xy)
    actual=evaluator[2].incident_spectrum(train,xy)
    assert_allclose(actual,expected,rtol=1e-11,atol=1e-11)


def test_source_application_matches_cpu_and_batches_memory(evaluator,monkeypatch):
    green,mesh,gpu=evaluator
    points=np.array([[0.,0.,0.],[1.,1.,0.],[2.,-3.,0.]])
    ks=np.array([.025,1.,2.,50.])
    source=np.random.default_rng(2).normal(size=(4,2,6)).astype(complex)+.4j
    expected=np.einsum('fpj,fjm->fpm',gpu.evaluate(points,ks),source)
    actual=gpu.apply_sources(points,ks,source)
    assert_allclose(actual,expected,rtol=1e-11,atol=1e-12)
    source *= 2
    assert_allclose(gpu.apply_sources(points,ks,source),2*expected,rtol=1e-11,atol=1e-12)
    # Bound the device allocation to force both spatial and frequency batches.
    actual_device=gpu.device
    monkeypatch.setattr(gpu,'device',type('SmallDevice',(),dict(max_mem_alloc_size=64,global_mem_size=512))())
    try:
        batched=gpu.apply_sources(points,ks,source)
    finally:
        gpu.device=actual_device
    assert_allclose(batched,2*expected,rtol=1e-11,atol=1e-12)


def test_close_releases_resources_and_is_idempotent(evaluator):
    from topwater.solve.gpu import OpenCLPotential
    green,mesh,_=evaluator
    gpu=OpenCLPotential(green,mesh,device_name=os.environ.get('TOPWATER_GPU_DEVICE'))
    gpu.evaluate(np.zeros((1,3)),np.array([1.]))
    gpu.close()
    gpu.close()
    assert gpu.context is None
    with pytest.raises(RuntimeError,match='closed'):
        gpu.evaluate(np.zeros((1,3)),np.array([1.]))


@pytest.mark.parametrize('points,ks', [
    (np.array([[0.,0.,.1]]),np.array([1.])),
    (np.array([[np.nan,0.,0.]]),np.array([1.])),
    (np.zeros((0,3)),np.array([1.])),
    (np.zeros((1,3)),np.array([0.])),
    (np.zeros((1,3)),np.array([np.nan])),
])
def test_unsupported_inputs_fail_explicitly(evaluator,points,ks):
    with pytest.raises(ValueError):
        evaluator[2].evaluate(points,ks)
