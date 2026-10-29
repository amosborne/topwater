"""Potential-only CPU prototype equivalence without GPU dependencies."""
from pathlib import Path
from types import SimpleNamespace

import capytaine as cpt
import numpy as np
from numpy.testing import assert_allclose
from numba import get_num_threads

from capytaine.green_functions.delhommeau import Delhommeau
from topwater.solve.cpu.potential import CPUPotential


def test_potential_only_cpu_preserves_near_far_and_thread_state():
    mesh=cpt.Mesh(vertices=np.array([[-1.,-1.,-.4],[1.,-1.,-.4],[1.,1.,-.4],[-1.,1.,-.4],
                                    [-1.,-1.,0.],[1.,-1.,0.],[1.,1.,0.],[-1.,1.,0.]]),
                  faces=np.array([[0,1,2,3],[4,5,6,7]]),auto_clean=False,auto_check=False)
    green=Delhommeau(tabulation_cache_dir=str(Path.home()/'.cache/topwater/capytaine-3.0.0'))
    points=np.array([[0.,0.,0.],[1.,1.,0.],[1.,0.,0.],[.123,.456,0.],
                     [20.,20.,0.],[1e3,0.,0.]])
    ks=np.array([1e-8,.025,.5,2.,50.,235.])
    expected=np.array([green.evaluate(points,mesh,wavenumber=float(k))[0] for k in ks])
    before=get_num_threads()
    evaluator=CPUPotential(green,mesh,threads=min(2,before))
    actual=evaluator.evaluate(points,ks)
    assert get_num_threads()==before
    assert_allclose(actual,expected,rtol=1e-10,atol=2e-13)
    assert_allclose(evaluator.evaluate(points[:2],ks),actual[:,:2],rtol=0,atol=0)
    source=np.random.default_rng(7).normal(size=(len(ks),2,6)).astype(complex)+.2j
    expected_fields=np.einsum('fpj,fjm->fpm',expected,source)
    evaluator.reserve_radiation(len(points),len(ks))
    fields=evaluator.apply_sources_transfer(points,ks,source)
    assert_allclose(fields.buffer,expected_fields,rtol=1e-10,atol=2e-13)
    source*=2
    assert_allclose(evaluator.apply_sources(points,ks,source),2*expected_fields,rtol=1e-10,atol=2e-13)
    assert get_num_threads()==before

    from topwater.solve.surface import incident_spectrum
    random=np.random.default_rng(21)
    coefficients=random.normal(size=(4,2,5))+1j*random.normal(size=(4,2,5))
    coefficients[:,:,0]=0
    coefficients[3]=0
    train=SimpleNamespace(theta=np.array([0.,.23,.23,2.]),k=np.array([0.,.1,2.,20.,200.]),
        omega=np.array([0.,1.,2.,3.,4.]),C=coefficients)
    xy=np.array([[-1.,0.],[0.,0.],[.1,.2],[20.,-10.]])
    assert_allclose(evaluator.incident_spectrum(train,xy),incident_spectrum(train,xy),
                    rtol=1e-11,atol=1e-11)
    assert get_num_threads()==before
