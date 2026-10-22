"""Actual moving-graph BEM + weak Zakharov RHS; contained architecture probe.

All potential values are graph nodes. Kinematics uses projected-area flux.
Dynamics integrates facet gradients and interpolated graph flux with P1 tests.
No nodal-to-P0 potential transfer or gradient averaging is used.
"""

import json
import argparse
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.linalg import solve, cho_factor, cho_solve
from threadpoolctl import threadpool_limits
from .domain import tank
from .operator import assemble_energy
from topwater import HullMesh
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes


@dataclass(frozen=True)
class Tank:
    mesh: HullMesh
    surface_nodes: np.ndarray
    wall_nodes: np.ndarray
    surface_faces: np.ndarray
    projected_areas: np.ndarray
    horizontal_gradients: np.ndarray
    projected_mass: np.ndarray
    projected_factor: tuple[np.ndarray, bool]
    shape_values: np.ndarray
    quadrature_weights: np.ndarray
    column_x: np.ndarray
    column_y: np.ndarray
    surface_shape: tuple[int, int]


def make_tank(resolution):
    mesh, free, node_map = tank(*resolution)
    fs = node_map.ravel()
    other = np.setdiff1d(np.arange(len(mesh.vertices)), fs)
    local = np.full(len(mesh.vertices), -1, dtype=int)
    local[fs] = np.arange(len(fs))
    faces = local[mesh.faces[free]]
    panels = TrianglePanels(mesh.vertices[mesh.faces[free]], 1)
    opposite = np.roll(panels.triangles, -2, axis=1) - np.roll(
        panels.triangles, -1, axis=1
    )
    gradient = np.cross(panels.geometric_normals[:, None, :], opposite) / (
        2 * panels.areas[:, None, None]
    )
    mass = np.zeros((len(fs), len(fs)))
    for face, area in zip(faces, panels.areas):
        mass[np.ix_(face, face)] += area * (np.ones((3, 3)) + np.eye(3)) / 12
    barycentric, weights = triangle_samples(np.eye(3)[None], 3)
    x = np.unique(mesh.vertices[:, 0])
    y = np.unique(mesh.vertices[:, 1])
    ix = np.searchsorted(x, mesh.vertices[:, 0])
    iy = np.searchsorted(y, mesh.vertices[:, 1])
    return Tank(
        mesh,
        fs,
        other,
        faces,
        panels.areas,
        gradient[:, :, :2],
        mass,
        cho_factor(mass),
        barycentric[0],
        weights,
        ix,
        iy,
        node_map.shape,
    )


def rhs(data, eta, psi, order=3):
    mesh, fs, other = data.mesh, data.surface_nodes, data.wall_nodes
    faces, areas, gradient = (
        data.surface_faces,
        data.projected_areas,
        data.horizontal_gradients,
    )
    mass, shape, weights = (
        data.projected_mass,
        data.shape_values,
        data.quadrature_weights,
    )
    ix, iy = data.column_x, data.column_y
    depth = 4.0
    g = 9.81
    if (
        not np.all(np.isfinite(eta))
        or not np.all(np.isfinite(psi))
        or np.min(eta) <= -depth
    ):
        raise ValueError("nonfinite state or graph touching the tank bottom")
    columns = eta.reshape(data.surface_shape)[ix, iy]
    vertices = mesh.vertices.copy()
    vertices[:, 2] += (vertices[:, 2] + depth) / depth * columns
    # Connectivity was validated once. Only panel geometry changes per stage.
    panels = TrianglePanels(vertices[mesh.faces], 1)
    energy, _ = assemble_energy(mesh.faces, panels, order)
    # Natural zero side/bottom flux. No constant-potential mean subtraction.
    wall_phi = -solve(
        energy[np.ix_(other, other)], energy[np.ix_(other, fs)] @ psi, assume_a="pos"
    )
    flux_load = energy[np.ix_(fs, fs)] @ psi + energy[np.ix_(fs, other)] @ wall_phi
    eta_t = cho_solve(data.projected_factor, flux_load)
    load = bernoulli_load(data, eta, psi, eta_t)
    psi_t = -g * eta + cho_solve(data.projected_factor, load)
    fluid_energy = 0.5 * (psi @ flux_load + g * eta @ mass @ eta)
    return eta_t, psi_t, float(np.ones(len(fs)) @ mass @ eta_t), float(fluid_energy)


def bernoulli_load(data, eta, psi, eta_t):
    """P1 weak nonlinear surface term on the fixed horizontal projection."""
    faces, areas, gradient = data.surface_faces, data.projected_areas, data.horizontal_gradients
    shape, weights = data.shape_values, data.quadrature_weights
    slope = np.einsum("fcj,fc->fj", gradient, eta[faces])
    psi_gradient = np.einsum("fcj,fc->fj", gradient, psi[faces])
    graph_flux = eta_t[faces] @ shape.T
    dot = np.sum(slope * psi_gradient, axis=1)
    nonlinear = (
        -0.5 * np.sum(psi_gradient**2, axis=1)[:, None]
        + 0.5
        * (graph_flux + dot[:, None]) ** 2
        / (1 + np.sum(slope**2, axis=1))[:, None]
    )
    load = np.zeros(len(eta))
    for corner in range(3):
        value = areas * np.einsum("q,q,fq->f", weights, shape[:, corner], nonlinear)
        np.add.at(load, faces[:, corner], value)
    return load


def run(resolution, amplitude, k=1.0):
    start = perf_counter()
    data = make_tank(resolution)
    result = check_quadratic(data, rhs, amplitude, k)
    result['resolution'] = resolution
    result['seconds'] = perf_counter() - start
    result['peak_rss_bytes'] = peak_rss_bytes()
    return result


def check_quadratic(data, evaluate, amplitude, k):
    """Independent finite-depth wave terms for a supplied evolving-domain solve."""
    mesh, fs = data.mesh, data.surface_nodes
    x = mesh.vertices[fs, 0] + 2 * np.pi
    g = 9.81
    k1 = k * np.tanh(4 * k)
    k2 = 2 * k * np.tanh(8 * k)
    omega = np.sqrt(g * k1)
    theta = np.pi / 4
    eta = amplitude * np.cos(k * x) * np.cos(theta)
    psi = -amplitude * g / omega * np.cos(k * x) * np.sin(theta)
    plus = evaluate(data, eta, psi)
    minus = evaluate(data, -eta, -psi)
    kin = (plus[0] + minus[0]) / (2 * amplitude**2)
    dyn = (plus[1] + minus[1]) / (2 * amplitude**2)
    exact_kin = (
        g
        / omega
        * (k1 * k2 / 2 - k * k)
        * np.cos(theta)
        * np.sin(theta)
        * np.cos(2 * k * x)
    )
    exact_dyn = (
        g
        * g
        / (4 * omega**2)
        * ((k1 * k1 - k * k) + (k1 * k1 + k * k) * np.cos(2 * k * x))
        * np.sin(theta) ** 2
    )
    samples, weights = triangle_samples(mesh.vertices[fs[data.surface_faces]], 16)
    sample_shape, _ = triangle_samples(np.eye(3)[None], 16)
    xx = samples[:, :, 0] + 2 * np.pi
    surface_kin = (
        g
        / omega
        * (k1 * k2 / 2 - k * k)
        * np.cos(theta)
        * np.sin(theta)
        * np.cos(2 * k * xx)
    )
    surface_dyn = (
        g
        * g
        / (4 * omega**2)
        * ((k1 * k1 - k * k) + (k1 * k1 + k * k) * np.cos(2 * k * xx))
        * np.sin(theta) ** 2
    )
    numerical_kin = kin[data.surface_faces] @ sample_shape[0].T
    numerical_dyn = dyn[data.surface_faces] @ sample_shape[0].T

    def surface_error(numerical, exact):
        return float(
            np.sqrt(
                np.einsum(
                    "f,q,fq->", data.projected_areas, weights, (numerical - exact) ** 2
                )
                / np.einsum("f,q,fq->", data.projected_areas, weights, exact**2)
            )
        )

    interior = (np.abs(mesh.vertices[fs, 0]) <= np.pi + 1e-12) & (
        np.abs(mesh.vertices[fs, 1]) <= np.pi / 2 + 1e-12
    )
    # The deep-water kinematic quadratic coefficient is tiny. Normalize by
    # natural k^2*g/omega scale, not by that cancellation.
    result = dict(
        qualified=False,
        amplitude=amplitude,
        wavenumber=k,
        quadratic_kinematic_normalized_rms=float(
            np.sqrt(np.mean((kin[interior] - exact_kin[interior]) ** 2))
            / (k * k * g / omega)
        ),
        quadratic_kinematic_relative_rms=float(
            np.sqrt(
                np.mean((kin[interior] - exact_kin[interior]) ** 2)
                / np.mean(exact_kin[interior] ** 2)
            )
        ),
        quadratic_dynamic_relative_rms=float(
            np.sqrt(
                np.mean((dyn[interior] - exact_dyn[interior]) ** 2)
                / np.mean(exact_dyn[interior] ** 2)
            )
        ),
        comparison_location="nodes in a fixed physical interior region",
        continuous_quadratic_kinematic_relative_rms=surface_error(
            numerical_kin, surface_kin
        ),
        continuous_quadratic_dynamic_relative_rms=surface_error(
            numerical_dyn, surface_dyn
        ),
        integrated_elevation_rates=[plus[2], minus[2]],
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fine",
        action="store_true",
        help="include the 5120-face whole-field finite-depth check",
    )
    args = parser.parse_args()
    cases = [
        ((8, 4, 4), 0.02, 1.0),
        ((16, 8, 8), 0.02, 1.0),
        ((16, 8, 8), 0.01, 1.0),
        ((8, 4, 4), 0.02, 0.25),
        ((16, 8, 8), 0.02, 0.25),
    ]
    if args.fine:
        cases.append(((32, 16, 16), 0.02, 0.25))
    rows = []
    with threadpool_limits(1):
        for resolution, amplitude, k in cases:
            row = run(resolution, amplitude, k)
            rows.append(row)
            print(json.dumps(row), flush=True)
    Path("artifacts").mkdir(exist_ok=True)
    Path("artifacts/free-surface-quadratic.json").write_text(
        json.dumps(rows, indent=2) + "\n"
    )
