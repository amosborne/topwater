import numpy as np
import pytest
from .common_output_sources import source_quality
from .frozen_source_quadrature import evaluate_trace


def test_source_gate_charges_both_adjacent_quadrature_uncertainties():
    A=np.zeros((3,6,6));A[:,0,0]=[0,8e-7,1.4e-6]
    U=np.zeros_like(A);U[:,0,0]=2e-7
    _,uncharged=source_quality(A,np.zeros_like(A),np.ones((6,6)))
    data,charged=source_quality(A,U,np.ones((6,6)))
    assert uncharged['passed'] and not charged['passed']
    assert data['source_previous_ratio'][0,0]>1
    np.testing.assert_allclose(data['previous_quadrature_charge'][0,0],4e-7)


def test_source_gate_reports_signed_increments_and_rejects_missing_uncertainty():
    A=np.zeros((3,6,6));A[:,1,5]=[2e-7,1e-7,0]
    data,gate=source_quality(A,np.zeros_like(A),np.ones((6,6)))
    assert gate['passed'] and data['source_previous_increment'][1,5]<0
    with pytest.raises(ValueError):source_quality(A,np.zeros((2,6,6)),np.ones((6,6)))
    with pytest.raises(ValueError):source_quality(A,-np.ones_like(A),np.ones((6,6)))


def test_shared_trace_evaluator_checks_deadline_around_every_native_block():
    calls=[]
    class Engine:
        def build_S_matrix(self,points,mesh,**kwargs):
            calls.append(('native',len(points)))
            return np.ones((len(points),2))
    def budget():calls.append(('budget',0))
    phi=evaluate_trace(Engine(),None,np.zeros((5,3)),np.ones((2,6)),2,budget)
    np.testing.assert_array_equal(phi,2)
    assert calls==[('budget',0),('native',2),('budget',0),('budget',0),('native',2),('budget',0),('budget',0),('native',1),('budget',0)]
