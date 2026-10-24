"""Source convention tests independent of the expensive external reference."""
import numpy as np
import pytest
from .periodic_sources import combine,integrate_pressure,require_close,rigid_normal_modes


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


def test_geometry_modes_match_pinned_capytaine_motion_objects_without_green_solve():
    pytest.importorskip('capytaine',reason='optional correlation API check; recovery command requires dependency')
    from tests.capytaine_reference import body_from_triangles,MODES
    triangles=np.array([[[0.,0.,-1.],[2.,0.,-1.],[0.,3.,-1.]],
                        [[0.,0.,-1.],[0.,3.,-1.],[0.,0.,-2.]]])
    center=np.array([.3,-.4,.2])
    body=body_from_triangles(triangles,center,2500.)
    expected=np.stack([np.einsum('ni,ni->n',body.dofs[m].evaluate_motion(body.mesh),
                        body.mesh.faces_normals) for m in MODES],axis=1)
    np.testing.assert_allclose(rigid_normal_modes(triangles,center),expected,atol=1e-15)
