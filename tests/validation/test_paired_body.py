import numpy as np
import pytest
from .paired_body import symmetric_family,scaled_check,entry_diagnostics,controls,capytaine_radius_check
from .reference_mesh import verify_surface
from .common import midpoint_children


def test_mirrored_nested_family_preserves_canonical_physical_facets():
    from tests.reference_boat import equilibrium
    _,loads=equilibrium();source=loads.wetted_triangles
    full,half,maximum=symmetric_family(source,refinements=(0,1,2))
    assert maximum<=1e-12
    assert [len(x) for x in full]==[len(full[0])*4**n for n in range(3)]
    for n,triangles in enumerate(full):
        verify_surface(source,triangles)
        # Compare unordered triangle vertices to avoid cyclic FFF choices.
        def canonical(values):return np.array([t[np.lexsort((t[:,2],t[:,1],t[:,0]))] for t in values])
        first=canonical(triangles[:len(half[n])]);second=canonical(triangles[len(half[n]):]*np.array([1,-1,1]))
        np.testing.assert_array_equal(first,second)
        if n:np.testing.assert_array_equal(half[n],midpoint_children(half[n-1]))


def test_asymmetric_physical_surface_cannot_be_silently_reflected():
    from tests.reference_boat import equilibrium
    _,loads=equilibrium();source=loads.wetted_triangles.copy();source[:,:,1]+=.1
    with pytest.raises((ValueError,AssertionError)):symmetric_family(source,refinements=(0,))


def test_infinity_trace_normalization_cancels_symbolic_factor():
    pytest.importorskip('capytaine',reason='optional installed-Capytaine API check')
    from capytaine.tools.symbolic_multiplication import SymbolicMultiplication
    omega=SymbolicMultiplication('∞',1.)
    chi=np.array([1.,2.,3.]);phi=(-1j*omega)*chi
    np.testing.assert_array_equal(np.asarray(phi/(-1j*omega)),chi)
    pressure=1j*omega*1025*phi
    np.testing.assert_array_equal(np.asarray(pressure/(omega*omega)),1025*chi)


def test_pinned_capytaine_radius_api_without_green_solve():
    pytest.importorskip('capytaine',reason='optional installed-Capytaine geometry API check')
    from tests.reference_boat import equilibrium
    from tests.capytaine_reference import body_from_triangles
    state,loads=equilibrium();family,_,_=symmetric_family(loads.wetted_triangles,refinements=(0,))
    body=body_from_triangles(family[0],state.position,2500.)
    assert capytaine_radius_check(body,family[0])<=1e-12


def test_trace_check_rejects_missing_and_nonfinite_values():
    assert scaled_check(np.ones((2,6)),np.ones((2,6)),1.,1e-10)==0
    for value in (np.full((2,6),np.nan),np.ones(3),-np.ones((2,6))):
        with pytest.raises(ValueError):scaled_check(value,np.ones((2,6)),1.,1e-10)


def test_all_entry_masks_expose_an_offdiagonal_failure():
    c=controls('validation/paired-body-controls.json')
    A=np.zeros((4,6,6));C=A.copy();A[1,0,4]=.001;A[2,0,4]=.0015;A[3,0,4]=.0018
    arrays,failures=entry_diagnostics(C,A,np.ones((6,6)),c)
    assert failures['reference_latest']==[[0,4]]
    assert failures['reference_previous']==[[0,4]]
    assert failures['candidate_comparison']==[[0,4]]
    assert not arrays['candidate_mesh_failed'].any()
