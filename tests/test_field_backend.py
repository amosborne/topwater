"""Reconstruction device selection and ownership without GPU hardware."""
from types import SimpleNamespace

import numpy as np
import pytest

from topwater.solve.field_backend import select_field_backend
from topwater.solve.field_backend import GPUUnavailableError
from topwater.solve.gpu import backend as gpu_module
from topwater.solve.types import MODES


def test_cpu_selection_does_not_construct_gpu(monkeypatch):
    from topwater.solve.cpu import backend as cpu
    backend=SimpleNamespace(evaluator=SimpleNamespace(threads=8))
    monkeypatch.setattr(cpu, 'CPUFieldBackend', lambda _, **kw: backend)
    monkeypatch.setattr(gpu_module,'GPUFieldBackend',lambda *args,**kw: pytest.fail('GPU constructed'))
    assert select_field_backend(object(),'cpu') is backend


@pytest.mark.parametrize('device,raises', [('auto',False),('gpu',True),('missing',True)])
def test_initialization_fallback_is_explicit(monkeypatch,device,raises,capsys):
    from topwater.solve.cpu import backend as cpu
    backend=SimpleNamespace(evaluator=SimpleNamespace(threads=8))
    monkeypatch.setattr(cpu, 'CPUFieldBackend', lambda _, **kw: backend)
    def unavailable(*args,**kwargs):
        raise GPUUnavailableError('driver unavailable')
    monkeypatch.setattr(gpu_module,'GPUFieldBackend',unavailable)
    if raises:
        with pytest.raises(GPUUnavailableError,match='driver unavailable'):
            select_field_backend(object(),device)
    else:
        assert select_field_backend(object(),device) is backend
        assert 'device=cpu reason=driver unavailable' in capsys.readouterr().err


def test_unexpected_initialization_bug_is_not_hidden(monkeypatch):
    def broken(*args,**kwargs):
        raise TypeError('implementation error')
    monkeypatch.setattr(gpu_module,'GPUFieldBackend',broken)
    with pytest.raises(TypeError,match='implementation error'):
        select_field_backend(object(),'auto')


def test_device_selection_returns_owned_backend(monkeypatch,capsys):
    backend=SimpleNamespace(device_name='test GPU',driver_version='1')
    def create(hydro,*,device_name):
        assert device_name=='test'
        return backend
    monkeypatch.setattr(gpu_module,'GPUFieldBackend',create)
    assert select_field_backend(object(),'test') is backend
    assert 'device=gpu name=test GPU' in capsys.readouterr().err


def test_adapter_reuses_source_bank_without_modifying_solver(monkeypatch):
    source=np.arange(12).reshape(2,6)+1j
    results={(w,mode):SimpleNamespace(sources=source[:,j])
             for w in (1.,2.) for j,mode in enumerate(MODES)}
    mesh=SimpleNamespace(nb_faces=2)
    original=lambda *args: None
    hydro=SimpleNamespace(depth=np.inf,body=SimpleNamespace(mesh_including_lid=mesh),
        solver=SimpleNamespace(engine=SimpleNamespace(green_function=object()),
                               compute_free_surface_elevation=original),
        radiation_results=results,gravity=9.81)
    class Evaluator:
        device=SimpleNamespace(name='test',driver_version='1')
        def __init__(self,*args,**kwargs):self.closed=False
        def apply_sources(self,points,ks,sources):
            return np.broadcast_to(sources.sum(axis=1)[:,None,:],(len(ks),len(points),sources.shape[2]))
        def apply_sources_transfer(self,points,ks,sources):
            from topwater.solve.gpu.potential import DeviceTransfer
            values=self.apply_sources(points,ks,sources)
            return DeviceTransfer(values,values.shape,None)
        def close(self):self.closed=True
    monkeypatch.setattr(gpu_module,'OpenCLPotential',Evaluator)
    monkeypatch.setattr('topwater.solve.gpu.fourier.GPUFourierWorkspace',
                        lambda _:SimpleNamespace(close=lambda:None))
    backend=gpu_module.GPUFieldBackend(hydro,cache_sources=False)
    actual=backend.radiation_transfer(hydro,np.zeros((3,2)),np.array([1.,2.]))
    expected=np.broadcast_to(-source.sum(axis=0)/9.81,(2,3,6))
    np.testing.assert_allclose(actual.buffer*actual.scale,expected)
    cached=backend._sources
    backend.radiation_transfer(hydro,np.zeros((1,2)),np.array([1.,2.]))
    assert backend._sources is cached
    assert hydro.solver.compute_free_surface_elevation is original
    backend.close()
    assert backend.evaluator.closed


def test_fft_initialization_failure_releases_evaluator(monkeypatch):
    closed = []
    evaluator = SimpleNamespace(close=lambda: closed.append(True))
    hydro = SimpleNamespace(depth=np.inf, body=SimpleNamespace(mesh_including_lid=object()),
        solver=SimpleNamespace(engine=SimpleNamespace(green_function=object())))
    monkeypatch.setattr(gpu_module, 'OpenCLPotential', lambda *args, **kwargs: evaluator)
    def fail(*args):
        raise RuntimeError('injected FFT initialization failure')
    monkeypatch.setattr('topwater.solve.gpu.fourier.GPUFourierWorkspace', fail)
    with pytest.raises(GPUUnavailableError, match='injected FFT'):
        gpu_module.GPUFieldBackend(hydro)
    assert closed == [True]


def test_backend_close_attempts_both_resources_and_preserves_first_error():
    backend = gpu_module.GPUFieldBackend.__new__(gpu_module.GPUFieldBackend)
    closed = []
    def fail(name):
        closed.append(name)
        raise RuntimeError(name)
    backend.fourier = SimpleNamespace(close=lambda: fail('fourier'))
    backend.evaluator = SimpleNamespace(close=lambda: fail('evaluator'))
    with pytest.raises(RuntimeError, match='fourier') as caught:
        backend.close()
    assert closed == ['fourier', 'evaluator']
    assert any('evaluator' in note for note in caught.value.__notes__)


@pytest.mark.parametrize('selector,expected,error', [
    ('test gpu', 0, None),
    ('TEST GPU PRO', 1, None),
    ('test', None, 'Ambiguous GPU name'),
    ('missing', None, 'available GPUs: Test GPU, Test GPU Pro'),
])
def test_gpu_name_selection_prefers_exact_and_reports_available_names(monkeypatch, selector, expected, error):
    import sys
    from topwater.solve import capytaine as adapter
    from topwater.solve.gpu.potential import OpenCLPotential
    devices = [SimpleNamespace(name=name, type=1, double_fp_config=1,
        host_unified_memory=False, max_compute_units=1, global_mem_size=1)
        for name in ('Test GPU', 'Test GPU Pro')]
    cl = SimpleNamespace(device_type=SimpleNamespace(GPU=1),
        # Drivers can enumerate the same model through multiple ICD entries.
        get_platforms=lambda: [SimpleNamespace(get_devices=lambda: devices*2)],
        Context=lambda devices: object(), CommandQueue=lambda *args, **kw: object(),
        command_queue_properties=SimpleNamespace(PROFILING_ENABLE=1))
    monkeypatch.setitem(sys.modules, 'pyopencl', cl)
    monkeypatch.setattr(adapter, 'potential_inputs', lambda *args: object())
    monkeypatch.setattr(OpenCLPotential, '_initialize', lambda *args: None)
    if error:
        with pytest.raises(RuntimeError, match=error):
            OpenCLPotential(None, None, device_name=selector)
    else:
        evaluator = OpenCLPotential(None, None, device_name=selector)
        assert evaluator.device is devices[expected]
