"""Weak damping must use consistent projection, including variable sigma."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.sparse.linalg import spsolve
from .experiments.free_surface.surface_damping import matrices


def test_linear_sigma_integrates_exactly_and_is_not_nodal_multiplication():
    points=np.array([[0.,0.,0.],[2.,0.,0.],[0.,1.,0.]])
    mass,weighted=matrices(points,np.array([[0,1,2]]),np.arange(3),np.array([0.,1.,2.]))
    # Integral lambda_i*lambda_j*lambda_k on the unit-area triangle:
    # 1/10 if all equal, 1/30 if two equal, 1/60 if all distinct.
    expected=np.array([[.1,1/15,1/12],[1/15,1/6,.1],[1/12,.1,7/30]])
    assert_allclose(weighted.toarray(),expected,atol=2e-15)
    psi=np.array([.7,-.3,.2])
    projected=spsolve(mass,weighted@psi)
    assert not np.allclose(projected,np.array([0.,1.,2.])*psi)
    assert_allclose(np.ones(3)@weighted@np.ones(3),1.)
