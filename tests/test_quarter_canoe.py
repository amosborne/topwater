"""Topology and quadrature ownership controls, without TetGen or a PDE solve."""
import numpy as np

from tests.experiments.free_surface import p4_tetra
from tests.experiments.free_surface.quarter_canoe import POCKET, make_plc, mirror
from tests.experiments.free_surface.quarter_canoe_control import rule, prolongation


def test_actual_canonical_quarter_plc_has_closed_winding_and_parent_ownership():
    plc = make_plc([0., POCKET[0], POCKET[0]+1.],
                   [0., POCKET[1], POCKET[1]+1.], [-4., -1., 0.], subdivisions=1)
    assert set(plc.markers) == {1, 2, 3, 4, 5}
    assert len(plc.body_parent) == np.sum(plc.markers == 1)
    p = plc.points[plc.faces]
    volume = np.einsum('fi,fi->', np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0]), p.mean(axis=1))/6
    assert 0 < volume < POCKET.prod()
    assert np.all(plc.points >= [0., 0., -1.])


def test_reflected_positive_tets_mirror_the_integration_points_and_parity_trace():
    # One cell touching both symmetry planes is enough to expose shared-node
    # ownership and orientation-permuted Duffy rules, including edge modes.
    points = np.array([[0., 0., -1.], [1., 0., -1.], [0., 1., -1.], [0., 0., 0.]])
    quarter = p4_tetra.make_mesh(points, np.array([[0, 1, 2, 3]]))
    full, maps, swaps, signs = mirror(quarter)
    q, _, weights = rule(5)
    for i, (sign, swap) in enumerate(zip(signs, swaps)):
        f, _, fw = rule(5, bool(swap))
        np.testing.assert_allclose(f@full.points[full.tets[i]], (q@points)*[*sign, 1], atol=2e-16)
        np.testing.assert_array_equal(fw, weights)
        np.testing.assert_array_equal(full.points[maps[i]], quarter.points*[*sign, 1])
    for parity in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
        t, unknown = prolongation(quarter, full, maps, signs, parity)
        value = t@np.ones(len(unknown))
        for axis in (0, 1):
            if parity[axis] == -1:
                np.testing.assert_array_equal(value[full.points[:, axis] == 0], 0.)
