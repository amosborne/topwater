"""Source convention tests independent of the expensive external reference."""
import numpy as np
import pytest
from .periodic_sources import combine,integrate_pressure,require_close


def test_complex_source_combination_preserves_phase_and_column_order():
    basis=np.arange(21).reshape(3,7)*(1+2j)
    xi=np.array([1,2j,3,-4j,5,6j]);zeta=-.02j
    expected=sum(basis[:,j]*xi[j] for j in range(6))+basis[:,6]*zeta
    np.testing.assert_allclose(combine(basis,zeta,xi),expected)
    with pytest.raises(ValueError):combine(basis[:,:6],zeta,xi)


def test_pressure_force_sign_and_world_cg_moment():
    triangles=np.array([[[0.,0.,0.],[2.,0.,0.],[0.,3.,0.]]])
    center=np.array([1.,2.,3.]);pressure=np.array([2.+3j])
    force=np.array([0,0,-3])*pressure[0]
    moment=np.cross(triangles[0].mean(axis=0)-center,force)
    np.testing.assert_allclose(integrate_pressure(triangles,center,pressure),np.r_[force,moment])


def test_reproduction_fails_on_phase_sign_and_nonfinite():
    expected=np.array([1+2j,0.])
    assert require_close(expected,expected,1.)==0
    for actual in (expected.conj(),np.array([np.nan,0]),np.ones(1)):
        with pytest.raises(ValueError):require_close(actual,expected,1.)
