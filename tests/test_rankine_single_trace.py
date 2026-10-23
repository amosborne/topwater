"""Owned S traces against independent edge integration and one-sided bulk limits."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.integrate import quad
from scipy.spatial.transform import Rotation

from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels, panel_influence, single_layer_trace


def triangles():
    return np.array([[[.1, -.2, 0.], [1.7, .1, 0.], [.3, 1.1, 0.]]])


@pytest.mark.parametrize('sign', [-1, 1])
@pytest.mark.parametrize('scale', [.3, 3.7])
def test_self_pv_edge_quadrature_and_oriented_jump(sign, scale):
    rotation = Rotation.from_rotvec([.2, -.4, .3]).as_matrix()
    panels = TrianglePanels(scale*triangles()@rotation.T+[.7, -.4, 1.2], sign)
    target = panels.centroids[0]
    pv = np.zeros(3)
    for edge in range(3):
        a, b = panels.triangles[0, edge], panels.triangles[0, (edge+1) % 3]
        length = np.linalg.norm(b-a)
        value, error = quad(lambda t: length/(4*np.pi*np.linalg.norm(target-(a+t*(b-a)))),
                            0., 1., epsabs=2e-14, epsrel=2e-13)
        assert error < 2e-13
        pv -= value*panels.edge_outward[0, edge]
    trace = single_layer_trace(panels, [0])
    assert trace.gradient_double_layer is None
    assert trace.double_layer[0, 0] == 0
    assert_allclose(trace.gradient_single_layer[0, 0], pv+.5*panels.normals[0], atol=3e-14, rtol=3e-13)
    scalar = panel_influence(panels.centroids, panels, self_indices=[0])
    assert_allclose(trace.single_layer, scalar.single_layer, atol=0, rtol=0)
    errors = []
    for step in [2e-4*scale, 1e-4*scale]:
        point = target-step*panels.normals[0]
        bulk = panel_influence(point[None], panels, gradients=True)
        errors.append(np.linalg.norm(bulk.gradient_single_layer[0, 0]-trace.gradient_single_layer[0, 0]))
    assert errors[1] < .51*errors[0]


def test_indexed_blocks_offdiagonal_and_existing_guards():
    panels = TrianglePanels(np.concatenate([triangles()+[4*j, 0., .3*j] for j in range(3)]), [-1, 1, -1])
    controls = LaplaceControls(target_block_size=2)
    indices = np.array([2, 0, 2])
    blocks = [single_layer_trace(panels, row, controls=controls) for row in [indices[:2], indices[2:]]]
    for result, row in zip(blocks, [indices[:2], indices[2:]], strict=True):
        for i, source in enumerate(row):
            # Subtract the one owned jump; all other source contributions have no jump.
            for other in range(3):
                if other == source:
                    continue
                single_panel = TrianglePanels(panels.triangles[other:other+1], panels.normal_sign[other])
                reference = panel_influence(panels.centroids[source:source+1], single_panel, gradients=True)
                assert_allclose(result.gradient_single_layer[i, other], reference.gradient_single_layer[0, 0], atol=2e-14)
    assert_allclose(blocks[0].gradient_single_layer[0], blocks[1].gradient_single_layer[0], atol=0, rtol=0)
    with pytest.raises(ValueError, match='target_block_size'):
        single_layer_trace(panels, indices, controls=controls)
    for bad in [[-1], [3], [.5], [True]]:
        with pytest.raises(ValueError, match='indices'):
            single_layer_trace(panels, bad)
    with pytest.raises(ValueError, match='trace limits'):
        panel_influence(panels.centroids[:1], panels, gradients=True, self_indices=[0])
    with pytest.raises(ValueError, match='self_indices'):
        panel_influence(panels.centroids[:1], panels, gradients=True)
    duplicate = TrianglePanels(np.repeat(triangles(), 2, axis=0))
    with pytest.raises(ValueError, match='self_indices'):
        single_layer_trace(duplicate, [0])
    for point, match in [(triangles()[0, 0], 'vertices'),
                         (triangles()[0, :2].mean(axis=0), 'edges')]:
        with pytest.raises(ValueError, match=match):
            panel_influence(point[None], TrianglePanels(triangles()), gradients=True)
