"""Same-physics kernel candidate checked against oracle and physical work."""
import numpy as np
import pytest

from .test_p4_geometry import curved_cell
from .experiments.free_surface.p4_geometry import (integrals as reference,
    InvalidGeometryError, kinetic_integrals as integrals, metrics)


def test_curved_nonuniform_flow_keeps_energy_and_full_shape_force():
    mesh, coordinates = curved_cell(.7)
    nodes = mesh.elements[0];x = coordinates[nodes]
    phi = np.sin(x[:, 0]+.3*x[:, 1]-.2*x[:, 2])
    old = reference(x, phi);new = integrals(x, phi)
    np.testing.assert_allclose(new['kinetic_energy'], old['kinetic_energy'], rtol=3e-14)
    np.testing.assert_allclose(new['kinetic_gradient'], old['kinetic_gradient'], rtol=2e-12, atol=2e-11)
    direction = np.column_stack((x[:, 0]*x[:, 1], x[:, 1]*x[:, 2], x[:, 0]*x[:, 2]))
    expected = np.sum(new['kinetic_gradient']*direction)
    errors = []
    for h in (1e-3, 5e-4):
        fd = (reference(x+h*direction, phi)['kinetic_energy']
              -reference(x-h*direction, phi)['kinetic_energy'])/(2*h)
        errors.append(abs(fd-expected))
    assert errors[1] < .4*errors[0]+1e-9
    assert errors[1] < 1e-7*max(1., abs(expected))


def test_physical_affine_flow_and_invalid_cell_guard_are_preserved():
    mesh, coordinates = curved_cell(.7);x = coordinates[mesh.elements[0]]
    expected = 500*(1/6+.7/24)
    for axis in range(3):
        assert integrals(x, x[:, axis])['kinetic_energy'] == pytest.approx(expected, abs=3e-10)
    inverted = x.copy();inverted[:, 0] *= -1
    with pytest.raises(InvalidGeometryError):
        metrics(inverted, 8)
    collapsed = x.copy();collapsed[:, 2] = 0
    with pytest.raises(InvalidGeometryError):
        integrals(collapsed, x[:, 2])
