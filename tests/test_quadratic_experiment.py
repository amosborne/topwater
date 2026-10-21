"""Independent source and physical checks for the flat P2/P1 experiment."""

import numpy as np
from numpy.testing import assert_allclose
import pytest
from threadpoolctl import threadpool_limits

from topwater.rankine import TrianglePanels
from .reference_hull import mesh_arrays
from .laplace_reference import triangle_samples
from .experiments.free_surface.quadratic import p2_shape, run, wave
from .experiments.free_surface.quadratic_kernel import source


def source_data(face):
    vertices, faces = mesh_arrays(3)
    panels = TrianglePanels(vertices[faces[[face]]], 1)
    opposite = np.roll(panels.triangles, -2, axis=1) - np.roll(
        panels.triangles, -1, axis=1
    )
    gradient = np.cross(panels.geometric_normals[:, None], opposite) / (
        2 * panels.areas[:, None, None]
    )
    far, weights = triangle_samples(panels.triangles, 8)
    bary, _ = triangle_samples(np.eye(3)[None], 8)

    def evaluate(x, self_pair=False):
        single, double = np.empty(3), np.empty(6)
        source(
            x,
            panels.triangles[0],
            panels.geometric_normals[0],
            panels.edge_outward[0],
            panels.tangents[0],
            panels.edge_lengths[0],
            panels.areas[0],
            panels.centroids[0],
            gradient[0],
            self_pair,
            far[0],
            weights,
            bary[0],
            single,
            double,
            np.empty((3, 3)),
            np.empty(3),
        )
        return single, double

    return panels, evaluate


@pytest.mark.parametrize("face", [11, 47])
@pytest.mark.parametrize("distance", [0.1, 1.0, 7.0, 20.0, 10000.0])
def test_canonical_source_moments_against_independent_quadrature(face, distance):
    panels, evaluate = source_data(face)
    point = panels.centroids[0] + distance * panels.diameters[0] * (
        panels.normals[0] + 0.3 * panels.tangents[0, 0]
    )
    single, double = evaluate(point)
    samples, weights = triangle_samples(panels.triangles, 96)
    bary, _ = triangle_samples(np.eye(3)[None], 96)
    delta = point - samples[0]
    radius = np.linalg.norm(delta, axis=1)
    g = panels.areas[0] * weights / (4 * np.pi * radius)
    dg = g * (delta @ panels.normals[0]) / radius**2
    # Constant-kernel family scaling avoids relative division by P2 vertex
    # integrals whose leading far-field term is genuinely zero.
    assert np.max(np.abs(single - g @ bary[0])) / g.sum() < 1e-9
    assert np.max(np.abs(double - dg @ p2_shape(bary[0]))) / abs(dg.sum()) < 5e-9


@pytest.mark.parametrize("face", [11, 47])
def test_self_integral_and_both_polynomial_trace_limits(face):
    panels, evaluate = source_data(face)
    location = np.array([0.2, 0.3, 0.5])
    point = location @ panels.triangles[0]
    single, double = evaluate(point, True)
    nodes, weights = np.polynomial.legendre.leggauss(256)
    t, weights = (nodes + 1) / 2, weights / 2
    expected = np.zeros(3)
    # Independent radial fan integration cancels the weak 1/r singularity;
    # the remaining line integrals have no numerical source moment recurrence.
    for edge in range(3):
        a, b = panels.triangles[0, edge], panels.triangles[0, (edge + 1) % 3]
        ray = (1 - t[:, None]) * (a - point) + t[:, None] * (b - point)
        edge_shape = (1 - t[:, None]) * np.eye(3)[edge] + t[:, None] * np.eye(3)[
            (edge + 1) % 3
        ]
        expected += (
            np.linalg.norm(np.cross(a - point, b - point))
            / (4 * np.pi)
            * np.sum(
                weights[:, None]
                * (edge_shape + location)
                / 2
                / np.linalg.norm(ray, axis=1)[:, None],
                axis=0,
            )
        )
    assert_allclose(single, expected, rtol=1e-9, atol=0)
    assert_allclose(double, 0, atol=0, rtol=0)
    for sign in (-1, 1):
        _, value = evaluate(
            point + sign * 1e-9 * panels.diameters[0] * panels.normals[0]
        )
        assert_allclose(
            value, sign * p2_shape(location[None])[0] / 2, atol=1e-7, rtol=0
        )


def test_complete_quadratic_space_energy_and_quadrature_refinement():
    with threadpool_limits(1):
        coarse, _ = run((8, 4, 4), 3)
        refined, _ = run((8, 4, 4), 6)
    for row in (coarse, refined):
        assert row["trace_dofs"] == 642
        assert row["flux_dofs"] == 960
        assert row["h_constant_max"] < 1e-11
        assert row["constant_q"] < 1e-10
        assert_allclose(row["net_flux"], 0, atol=1e-11, rtol=0)
        assert abs(row["lowest_eigenvalues"][0]) < 1e-10
        assert row["lowest_eigenvalues"][1] == pytest.approx(
            0.25 * np.tanh(1.0), rel=0.001
        )
    # All eight harmonic degree <=2 traces are represented exactly. Their
    # cross-energy matrix is integrated independently in the actual volume.
    assert (
        refined["harmonic_degree_two_energy_matrix_error"]
        < coarse["harmonic_degree_two_energy_matrix_error"] / 8
    )
    assert refined["harmonic_degree_two_energy_matrix_error"] < 1e-4


def test_five_period_wave_at_matched_surface_node_count():
    with threadpool_limits(1):
        linear = wave((16, 8, 8), degree=1)
        quadratic = wave((8, 4, 4), degree=2)
    assert linear["surface_dofs"] == quadratic["surface_dofs"] == 153
    # Actual DtN plus RK23, compared with the continuous analytic wave over
    # the whole surface and five periods. This is a method comparison, not a
    # one-percent gate or nonlinear/curved-geometry qualification.
    assert (
        quadratic["continuous_elevation_relative_l2_time"]
        < linear["continuous_elevation_relative_l2_time"] * 0.6
    )
    assert quadratic["continuous_elevation_relative_l2_time"] < 0.04
    assert quadratic["energy_drift"] < 1e-6
    assert quadratic["volume_drift"] < 1e-12
