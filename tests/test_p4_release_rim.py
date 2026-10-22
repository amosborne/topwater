"""Physical polynomial overshoot control for the conservative P4 rim guard."""
from types import SimpleNamespace

import numpy as np
import pytest

from .experiments.free_surface.p4_tetra import indices
from .experiments.free_surface.canoe_contact import ContactTopologyError
from .experiments.free_surface.p4_canoe_release import _rim_bound


def test_global_polynomial_bound_catches_height_missed_by_every_nodal_sample():
    # A real graph z=lambda0*lambda1*lambda2 on the unit projected triangle.
    # AM-GM gives its exact physical maximum1/27 at the centroid; no P4 node
    # is the centroid. Degree elevation3->4 gives Bernstein max1/12.
    barycentric = indices(3)/4
    height = np.prod(barycentric, axis=1)
    coordinates = np.column_stack((barycentric[:, 1:], height))
    assert height.max() == 1/32
    dense = np.array([[i/60, j/60, 1-(i+j)/60]
                      for i in range(61) for j in range(61-i)])
    assert np.prod(dense, axis=1).max() == pytest.approx(1/27, abs=2e-17)
    # Explicit straight flat physical opening rim, with its projection inside
    # this graph. This supplies geometry only; no PDE or full-canoe mock.
    rim = np.array([[.2, .2, .1], [.5, .2, .1], [.4, .4, .1]])
    owner = SimpleNamespace(chart=SimpleNamespace(linear=SimpleNamespace(
        rim_vertices=rim, cg=np.zeros(3))))
    args = (owner, coordinates, np.zeros(3), np.array([1., 0., 0., 0.]),
            np.arange(15)[None], (np.array([-1., -1.]), np.array([2., 2.])))
    clearance, rim_z, upper, margin = _rim_bound(*args)
    assert 1/12 <= upper < 1/12+2e-15
    assert rim_z == .1 and margin > 0
    assert 0 < clearance <= .1-1/12
    # Nodal max would falsely certify .034>1/32 even though the graph's
    # actual maximum exceeds the rim. The conservative guard must stop.
    rim[:, 2] = .034
    assert height.max() < rim[0, 2] < 1/27
    with pytest.raises(ContactTopologyError, match='not a flooding classification'):
        _rim_bound(*args)
    # A bound can also be inconclusive for a safe rim: never misclassify that
    # numerical stop as physical flooding or silently loosen the bound.
    rim[:, 2] = .05
    assert 1/27 < rim[0, 2] < 1/12
    with pytest.raises(ContactTopologyError, match='cannot certify'):
        _rim_bound(*args)
    rim[:, 2] = .1
    rim[0, 0] = 2.
    with pytest.raises(ContactTopologyError, match='outside the certified outer tank box'):
        _rim_bound(*args)
