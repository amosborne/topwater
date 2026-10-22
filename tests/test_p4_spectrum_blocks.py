"""Curved-cell restriction check for the actual-canoe spectrum capture.

The all-free mask is an algebraic edge case, not a physical water cell.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from tests.experiments.free_surface.p4_canoe_spectrum import surface_blocks
from tests.experiments.free_surface.p4_geometry import stiffness
from tests.experiments.free_surface.p4_tetra import indices


@pytest.mark.parametrize('free', (np.arange(0, 35, 4), np.arange(35)))
def test_curved_surface_blocks_equal_full_stiffness_restriction(free):
    reference = indices(4)[:, 1:]/4
    points = reference@np.array([[1.1, .1, -.03], [.05, .8, .04], [0., .02, 1.2]])
    points[:, 2] += .08*reference[:, 0]**2+.03*reference[:, 0]*reference[:, 1]
    # A disconnected cell with no free nodes must not enter either block.
    points = np.vstack((points, points+[3., 1., 0.]))
    mesh = SimpleNamespace(points=points, elements=np.arange(70).reshape(2, 35))
    owner = SimpleNamespace(mesh=mesh)
    inside = np.setdiff1d(np.arange(70), free)
    kif, kff, touched = surface_blocks(owner, points, inside, free)
    full = stiffness(mesh, points)
    assert touched == 1
    for candidate, expected in ((kif, full[inside][:, free]), (kff, full[free][:, free])):
        np.testing.assert_allclose(candidate.toarray(), expected.toarray(), rtol=2e-13, atol=2e-14)
