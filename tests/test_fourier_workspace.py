from dataclasses import replace
import numpy as np

from topwater.solve.cpu.fourier import CPUFourierWorkspace as FourierWorkspace
from topwater.solve.fourier import offline_radiation_response
from .test_offline_radiation import compact_case
from topwater.solve.transfer import Transfer


def test_shared_workspace_matches_reference_and_invalidates_motion():
    t,v,a,w,H,settings=compact_case()
    workspace=FourierWorkspace()
    for multiplier in (1.,1.,2.):
        expected,evidence=offline_radiation_response(w,H,t,v*multiplier,a,settings=settings)
        actual,reused=offline_radiation_response(w,H,t,v*multiplier,a,settings=settings,workspace=workspace)
        np.testing.assert_allclose(actual,expected,atol=1e-13,rtol=1e-13)
        np.testing.assert_allclose(reused['raw_initial_radiation_m'],evidence['raw_initial_radiation_m'],atol=1e-13)


def test_workspace_keeps_complex_phase_on_native_odd_clock():
    from topwater.solve.settings import OfflineRadiationSettings
    t=np.arange(129)*.01
    v=np.zeros((len(t),6));v[:,0]=np.sin(2*np.pi*t)*t
    a=np.zeros_like(v)
    w=2*np.pi*np.fft.rfftfreq(len(t),.01)
    H=np.zeros((len(w),3,6),complex);H[:,:,0]=np.exp(.03j*w[:,None])
    settings=OfflineRadiationSettings()
    expected,_=offline_radiation_response(w,H,t,v,a,settings=settings)
    actual,_=offline_radiation_response(w,H,t,v,a,settings=settings,workspace=FourierWorkspace())
    np.testing.assert_allclose(actual,expected,atol=1e-13,rtol=1e-13)


def test_cpu_lazy_zero_and_scale_preserve_transfer_and_reuse_adjacent_lengths():
    random=np.random.default_rng(93)
    workspace=FourierWorkspace()
    for n in (128,129):
        time=np.arange(n)*.01
        w=2*np.pi*np.fft.rfftfreq(n,.01)
        values=random.normal(size=(len(w)-1,3,6))+1j*random.normal(size=(len(w)-1,3,6))
        transfer=Transfer(values,values.shape).scaled(-1/9.81).with_zero()
        explicit=np.concatenate([np.zeros((1,3,6)),values*(-1/9.81)])
        velocity=random.normal(size=(n,6));velocity[0]=0
        actual,_=offline_radiation_response(w,transfer,time,velocity,np.zeros_like(velocity),workspace=workspace)
        # Independent full-bin matrix expression, including even-size Nyquist.
        spectrum=np.einsum('fpm,fm->fp',explicit.conj(),np.fft.rfft(velocity,axis=0))
        np.testing.assert_allclose(actual,np.fft.irfft(spectrum,n=n,axis=0),atol=1e-13,rtol=1e-13)
