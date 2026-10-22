import json

import numpy as np
from numpy.testing import assert_allclose
import pytest

from topwater.frames import rotation_matrix
from topwater import HullMesh
from topwater.bem import LaplaceOperators, assemble_laplace
from topwater.numerics import LaplaceControls
from topwater.rankine import QuadratureError, TrianglePanels, _fallback, panel_influence
from .rankine_reference import panel_reference
from .reference_hull import mesh_arrays


TRIANGLE = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)


def test_equilateral_self_integral_and_explicit_jump_exclusion():
    triangle = np.array([[0, 0, 0], [1, 0, 0], [0.5, np.sqrt(3)/2, 0]])
    panels = TrianglePanels(triangle[None])
    result = panel_influence(panels.centroids, panels, self_indices=[0])
    expected = np.sqrt(3)*np.log(2+np.sqrt(3))/(4*np.pi)
    assert_allclose(result.single_layer[0, 0], expected, rtol=1e-13)
    assert result.double_layer[0, 0] == 0
    reference = panel_reference(panels.centroids[0], triangle)
    assert_allclose(result.single_layer[0, 0], reference[0], rtol=1e-13)


@pytest.mark.parametrize("height", [10, 1, 1e-1, 1e-3, 1e-6, -1e-1, -1e-6])
@pytest.mark.parametrize("projection", [[1/3, 1/3], [0.5, 1e-3]])
def test_self_near_far_and_bulk_gradients_against_independent_quadrature(height, projection):
    target = np.r_[projection, height]
    panels = TrianglePanels(TRIANGLE[None])
    result = panel_influence(target[None], panels, gradients=True)
    s, k, gs, gk = panel_reference(target, TRIANGLE)
    diameter = panels.diameters[0]
    assert_allclose(result.single_layer[0, 0]/diameter, s/diameter, rtol=1e-8, atol=1e-10)
    assert_allclose(result.double_layer[0, 0], k, rtol=1e-8, atol=1e-10)
    assert_allclose(result.gradient_single_layer[0, 0], gs, rtol=1e-7, atol=1e-9)
    assert_allclose(result.gradient_double_layer[0, 0]*diameter, gk*diameter, rtol=1e-7, atol=1e-9)


def test_normal_reversal_is_not_a_source_target_derivative_swap():
    target = [[0.2, 0.3, 0.4]]
    forward = panel_influence(target, TrianglePanels(TRIANGLE[None]), gradients=True)
    reverse = panel_influence(target, TrianglePanels(TRIANGLE[None], -1), gradients=True)
    assert_allclose(reverse.single_layer, forward.single_layer)
    assert_allclose(reverse.gradient_single_layer, forward.gradient_single_layer)
    assert_allclose(reverse.double_layer, -forward.double_layer)
    assert_allclose(reverse.gradient_double_layer, -forward.gradient_double_layer)
    rewound = panel_influence(target, TrianglePanels(TRIANGLE[None, ::-1]), gradients=True)
    assert_allclose(rewound.single_layer, forward.single_layer)
    assert_allclose(rewound.double_layer, -forward.double_layer)


def test_rigid_covariance_and_length_scaling():
    q = np.array([1, 2, 3, 4])/np.sqrt(30)
    rotation = rotation_matrix(q)
    translation = np.array([4, -5, 6])
    target = np.array([[0.2, 0.3, 0.4]])
    reference = panel_influence(target, TrianglePanels(TRIANGLE[None]), gradients=True)
    for scale in (1e-3, 1, 1e3):
        panels = TrianglePanels((scale*TRIANGLE@rotation.T+translation)[None])
        result = panel_influence(scale*target@rotation.T+translation, panels, gradients=True)
        assert_allclose(result.single_layer/scale, reference.single_layer, rtol=1e-9, atol=1e-10)
        assert_allclose(result.double_layer, reference.double_layer, rtol=1e-9, atol=1e-10)
        assert_allclose(result.gradient_single_layer@rotation, reference.gradient_single_layer,
                        rtol=1e-7, atol=1e-9)
        assert_allclose(scale*result.gradient_double_layer@rotation, reference.gradient_double_layer,
                        rtol=1e-7, atol=1e-9)


def test_boundary_and_oversized_target_requests_fail_explicitly():
    panels = TrianglePanels(TRIANGLE[None])
    for target in ([[0, 0, 0]], [[0.5, 0, 0]], panels.centroids):
        with pytest.raises(ValueError):
            panel_influence(target, panels)
    with pytest.raises(ValueError, match="boundary gradients"):
        panel_influence(panels.centroids, panels, self_indices=[0], gradients=True)
    with pytest.raises(ValueError, match="bounded blocks"):
        panel_influence(np.ones((129, 3)), panels)
    with pytest.raises(ValueError, match="centroid"):
        panel_influence([[0.2, 0.2, 1]], panels, self_indices=[0])


def test_cached_geometry_owns_arrays_and_rejects_degeneracy():
    original = TRIANGLE[None].copy()
    panels = TrianglePanels(original)
    original[:] = 0
    assert_allclose(panels.triangles[0], TRIANGLE)
    assert not panels.areas.flags.writeable
    with pytest.raises(ValueError, match="nondegenerate"):
        TrianglePanels(np.zeros((1, 3, 3)))


def test_coplanar_exterior_normal_derivative_is_nonzero():
    panels = TrianglePanels(TRIANGLE[None])
    target = np.array([2., 2., 0.])
    result = panel_influence(target[None], panels, gradients=True)
    reference = panel_reference(target, TRIANGLE)
    assert_allclose(result.gradient_double_layer[0, 0], reference[3], atol=1e-12)
    assert result.gradient_double_layer[0, 0, 2] > 0.003
    step = np.array([0, 0, 1e-5])
    difference = (panel_influence((target+step)[None], panels).double_layer
                  - panel_influence((target-step)[None], panels).double_layer)/(2e-5)
    assert_allclose(difference[0, 0], reference[3][2], rtol=1e-9)


def test_cancellation_fallback_and_explicit_failure_for_unresolved_gradients():
    panels = TrianglePanels(TRIANGLE[None])
    controls = LaplaceControls(quadrature_rtol=1e-12, quadrature_atol=1e-30)
    result = panel_influence([[1e8, 3e7, 4e7]], panels, controls=controls, gradients=True)
    assert result.fallback_pairs == 1
    # At this separation the leading monopole gives an independent 1/r check.
    delta = np.array([1e8, 3e7, 4e7])-panels.centroids[0]
    assert_allclose(result.single_layer[0, 0], .5/(4*np.pi*np.linalg.norm(delta)), rtol=1e-14)
    # Solid-angle limits at a right-angle corner and one height inside an edge.
    # Generic adaptive quadrature can miss the 1e-12-wide peak here.
    for target, expected in (([0, 0, 1e-12], .125), ([.5, 1e-12, 1e-12], .375)):
        actual, _ = _fallback(np.array(target), TRIANGLE, np.array([0, 0, 1]), LaplaceControls(), False)
        assert_allclose(actual[1], expected, atol=1e-11)
    result = panel_influence([[.5, 1e-12, 1e-12]], panels)
    assert result.fallback_pairs == 1
    with pytest.raises(QuadratureError, match="close panel"):
        panel_influence([[.5, -1e-12, 0]], panels, gradients=True)
    with pytest.raises(QuadratureError, match="max_refinement"):
        _fallback(np.array([3, 3, 3]), TRIANGLE, np.array([0, 0, 1]),
                  LaplaceControls(max_refinement=0), False)


def test_closed_canonical_hull_interior_exterior_identities_and_actual_result_roundtrip(tmp_path):
    vertices, faces = mesh_arrays()
    mesh = HullMesh(vertices, faces)
    panels = TrianglePanels.from_mesh(mesh)
    operators = assemble_laplace(panels)
    assert_allclose(operators.double_layer.sum(axis=1), -0.5, atol=1e-8)
    assert_allclose(np.diag(operators.double_layer), 0, atol=0)
    reversed_operators = assemble_laplace(TrianglePanels.from_mesh(mesh, -1))
    assert_allclose(reversed_operators.double_layer.sum(axis=1), 0.5, atol=1e-8)
    assert_allclose(reversed_operators.single_layer, operators.single_layer, atol=0)
    # Interior/exterior evaluation identities are independent of a linear solve.
    points = np.array([[0, 0, -0.175], [0, 0, 1]])
    field = panel_influence(points, panels, gradients=True)
    assert_allclose(field.double_layer.sum(axis=1), [-1, 0], atol=1e-8)
    assert_allclose(field.gradient_double_layer.sum(axis=1), 0, atol=1e-8)
    path = tmp_path / "operators.npz"
    operators.save(path)
    restored = LaplaceOperators.load(path)
    assert_allclose(restored.single_layer, operators.single_layer, rtol=0, atol=0)
    assert_allclose(restored.double_layer, operators.double_layer, rtol=0, atol=0)
    assert restored.diagnostics == operators.diagnostics
    assert operators.diagnostics["fallback_pairs"] == 0


@pytest.mark.parametrize("corruption", ["trace", "count", "shape", "dtype", "nonfinite", "budget"])
def test_saved_operator_rejects_corrupt_headers_metadata_and_budget(tmp_path, corruption):
    operators = assemble_laplace(TrianglePanels(TRIANGLE[None]))
    metadata = operators.diagnostics.copy()
    single = operators.single_layer.copy()
    double = operators.double_layer.copy()
    controls = LaplaceControls()
    if corruption == "trace":
        metadata["double_layer_trace"] = "includes_jump"
    elif corruption == "count":
        metadata["panel_count"] = 2
    elif corruption == "shape":
        single = np.zeros((2, 3))
    elif corruption == "dtype":
        single = single.astype(np.float32)
    elif corruption == "nonfinite":
        single[0, 0] = np.nan
    elif corruption == "budget":
        # Actual headers exceed budget despite metadata claiming a single panel.
        single = double = np.ones((3, 3))
        controls = LaplaceControls(max_operator_bytes=64)
    path = tmp_path / "corrupt.npz"
    np.savez(path, single_layer=single, double_layer=double, metadata=json.dumps(metadata))
    with pytest.raises(MemoryError if corruption == "budget" else ValueError):
        LaplaceOperators.load(path, controls=controls)
