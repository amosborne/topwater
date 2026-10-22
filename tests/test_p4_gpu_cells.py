"""Explicit optional-GPU physical cell controls, no global fluid solve."""
import numpy as np
import pytest

pytest.importorskip('cupy')

from .experiments.free_surface import cuda_operator  # Existing Windows CUDA DLL setup.
from .experiments.free_surface.p4_geometry import integrals, element_stiffness, InvalidGeometryError
from .experiments.free_surface.p4_gpu_cells import cell_terms
from .test_p4_geometry import curved_cell


def test_gpu_curved_physical_energy_force_and_stiffness():
    mesh, coordinates = curved_cell(.7)
    x = coordinates[mesh.elements[0]]
    phi = x[:, 0]*x[:, 1]
    result = cell_terms(x[None], phi[None], chunk_size=1)
    reference = integrals(x, phi)
    assert result['kinetic_energy'][0] == pytest.approx(500*(1/30+.7/90), abs=3e-10)
    np.testing.assert_allclose(result['kinetic_gradient'][0], reference['kinetic_gradient'], rtol=3e-12, atol=3e-10)
    np.testing.assert_allclose(result['element_stiffness'][0], element_stiffness(x), rtol=3e-12, atol=3e-12)
    inverted = x.copy();inverted[:, 0] *= -1
    with pytest.raises(InvalidGeometryError):
        cell_terms(inverted[None], phi[None], chunk_size=1)
