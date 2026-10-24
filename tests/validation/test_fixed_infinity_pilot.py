import numpy as np
import pytest
from .fixed_infinity_pilot import validate_infinite_basis


def test_infinite_basis_rejects_missing_mode_finite_seventh_column_and_complex_trace():
    valid=np.ones((3,6),complex);validate_infinite_basis(valid,valid,valid,3,6)
    with pytest.raises(ValueError):validate_infinite_basis(valid,valid,valid,3,5)
    with pytest.raises(ValueError):validate_infinite_basis(valid,np.ones((3,7)),valid,3,6)
    with pytest.raises(ValueError):validate_infinite_basis(valid+1j,valid,valid,3,6)
    damaged=valid.copy();damaged[0,0]=np.nan
    with pytest.raises(ValueError):validate_infinite_basis(damaged,valid,valid,3,6)
