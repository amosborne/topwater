import os
from pathlib import Path
import numpy as np
import pytest

from topwater.solve.fourier import offline_radiation_response
from topwater.solve.settings import OfflineRadiationSettings

pytestmark=pytest.mark.skipif(os.environ.get('TOPWATER_TEST_GPU')!='1',reason='explicit GPU opt-in required')


@pytest.fixture(scope='module')
def workspace():
    import capytaine as cpt
    from capytaine.green_functions.delhommeau import Delhommeau
    from topwater.solve.gpu.potential import OpenCLPotential
    from topwater.solve.gpu.fourier import GPUFourierWorkspace
    mesh=cpt.Mesh(vertices=np.array([[-1.,-1.,-.2],[1.,-1.,-.2],[1.,1.,-.2],[-1.,1.,-.2]]),
                  faces=np.array([[0,1,2,3]]),auto_clean=False,auto_check=False)
    green=Delhommeau(tabulation_cache_dir=str(Path.home()/'.cache/topwater/capytaine-3.0.0'))
    evaluator=OpenCLPotential(green,mesh,device_name=os.environ.get('TOPWATER_GPU_DEVICE'))
    value=GPUFourierWorkspace(evaluator)
    yield value
    value.close();evaluator.close()


@pytest.mark.parametrize('n',[17,64,129,13001])
def test_gpu_arbitrary_length_real_fourier_preserves_full_complex_response(workspace,n):
    t=np.arange(n)*.01
    random=np.random.default_rng(n)
    v=random.normal(size=(n,6));v[0]=0
    a=np.zeros_like(v)
    w=2*np.pi*np.fft.rfftfreq(len(t),.01)
    H=random.normal(size=(len(w),3,6))+1j*random.normal(size=(len(w),3,6))
    output=t[::8]
    settings=OfflineRadiationSettings()
    expected,evidence=offline_radiation_response(w,H,t,v,a,settings=settings,output_time_s=output)
    actual,gpu=offline_radiation_response(w,H,t,v,a,settings=settings,output_time_s=output,workspace=workspace)
    np.testing.assert_allclose(actual,expected,atol=1e-11,rtol=1e-11)
    np.testing.assert_allclose(gpu['raw_initial_radiation_m'],evidence['raw_initial_radiation_m'],atol=1e-11)


def test_resident_transfer_virtual_zero_and_scaling(workspace):
    from topwater.solve.gpu.potential import DeviceTransfer
    n=129;t=np.arange(n)*.01
    v=np.zeros((n,6));v[:,0]=t*np.sin(t)
    a=np.zeros_like(v);w=2*np.pi*np.fft.rfftfreq(len(t),.01)
    H=np.zeros((len(w),2,6),complex);H[1:,:,0]=np.exp(.03j*w[1:,None])
    device=DeviceTransfer(workspace.upload(H[1:]),H[1:].shape,workspace.context).scaled(.7).with_zero()
    expected,_=offline_radiation_response(w,.7*H,t,v,a)
    actual,_=offline_radiation_response(w,device,t,v,a,workspace=workspace)
    np.testing.assert_allclose(actual,expected,atol=1e-11,rtol=1e-11)


def test_gpu_preserves_historical_padded_pre_response(workspace):
    from .test_offline_radiation import compact_case
    t,v,a,w,H,settings=compact_case()
    expected,evidence=offline_radiation_response(w,H,t,v,a,settings=settings)
    actual,gpu=offline_radiation_response(w,H,t,v,a,settings=settings,workspace=workspace)
    np.testing.assert_allclose(actual,expected,atol=1e-11,rtol=1e-11)
    np.testing.assert_allclose(gpu['pre_response_max_abs_m'],evidence['pre_response_max_abs_m'],atol=1e-11)
