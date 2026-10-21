"""Contained conforming P2/P1 flat-wave comparison; not a runtime backend."""

import argparse
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.linalg import solve, eigh
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix
from threadpoolctl import threadpool_limits
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes
from .domain import tank, trace
from .kernel import symmetrize_single_layer
from .quadratic_kernel import assemble_pairs


def p2_shape(bary):
    return np.column_stack(
        (
            bary * (2 * bary - 1),
            4 * bary[:, 0] * bary[:, 1],
            4 * bary[:, 1] * bary[:, 2],
            4 * bary[:, 2] * bary[:, 0],
        )
    )


def assemble(mesh, panels, order=3):
    """Closed bounded-domain energy using geometric outward normals (+1)."""
    nf = len(panels)
    edges = np.sort(mesh.faces[:, [(0, 1), (1, 2), (2, 0)]].reshape(-1, 2), axis=1)
    edges, index = np.unique(edges, axis=0, return_inverse=True)
    dofs = np.column_stack((mesh.faces, len(mesh.vertices) + index.reshape(-1, 3)))
    coords = np.vstack((mesh.vertices, mesh.vertices[edges].mean(axis=1)))
    nv = len(coords)
    opposite = np.roll(panels.triangles, -2, axis=1) - np.roll(
        panels.triangles, -1, axis=1
    )
    gradient = np.cross(panels.geometric_normals[:, None], opposite) / (
        2 * panels.areas[:, None, None]
    )
    points, weights = triangle_samples(panels.triangles, order)
    bary, _ = triangle_samples(np.eye(3)[None], order)
    far_points, far_weights = triangle_samples(panels.triangles, 8)
    far_shape, _ = triangle_samples(np.eye(3)[None], 8)
    v, h = assemble_pairs(
        panels.triangles,
        panels.geometric_normals,
        panels.edge_outward,
        panels.tangents,
        panels.edge_lengths,
        panels.areas,
        panels.centroids,
        gradient,
        dofs,
        points,
        weights,
        bary[0],
        far_points,
        far_weights,
        far_shape[0],
        nv,
    )
    asym = float(symmetrize_single_layer(v))
    # Nodal coefficients of the exact linear curl of each P2 basis function.
    grads = np.empty((nf, 3, 6, 3))
    for a in range(3):
        for i in range(3):
            grads[:, a, i] = (4 * (a == i) - 1) * gradient[:, i]
        for j, (i, k) in enumerate(((0, 1), (1, 2), (2, 0))):
            grads[:, a, 3 + j] = 4 * (
                (a == i) * gradient[:, k] + (a == k) * gradient[:, i]
            )
    curl = np.cross(panels.geometric_normals[:, None, None], grads)
    energy = np.zeros((nv, nv))
    rows = np.repeat(np.arange(3 * nf), 6)
    cols = np.broadcast_to(dofs[:, None], (nf, 3, 6)).ravel()
    for axis in range(3):
        c = coo_matrix(
            (curl[:, :, :, axis].ravel(), (rows, cols)), shape=(3 * nf, nv)
        ).tocsr()
        energy += c.T @ (v @ c)
    auxiliary = solve(v.T, h, assume_a="pos", overwrite_a=True)
    energy += h.T @ auxiliary
    return (
        energy,
        coords,
        dofs,
        dict(
            asymmetry=asym,
            h_constant_max=float(np.max(np.abs(h @ np.ones(nv)))),
            flux_dofs=3 * nf,
            trace_dofs=nv,
            operator_bytes=dict(
                single_layer_factor=v.nbytes, mixed_trace=h.nbytes, energy=energy.nbytes
            ),
        ),
    )


def run(resolution=(8, 4, 4), order=3, degree=2):
    start = perf_counter()
    mesh, free, _ = tank(*resolution)
    panels = TrianglePanels.from_mesh(mesh)
    if degree == 2:
        energy, coords, dofs, diagnostics = assemble(mesh, panels, order)
    else:
        from .operator import assemble_energy

        energy, diagnostics = assemble_energy(mesh.faces, panels, order)
        coords, dofs = mesh.vertices, mesh.faces
    fs = np.unique(dofs[free])
    other = np.setdiff1d(np.arange(len(coords)), fs)
    st = energy[np.ix_(fs, fs)] - energy[np.ix_(fs, other)] @ solve(
        energy[np.ix_(other, other)], energy[np.ix_(other, fs)], assume_a="pos"
    )
    lookup = np.full(len(coords), -1, dtype=int)
    lookup[fs] = np.arange(len(fs))
    surface_faces = lookup[dofs[free]]
    samples, weights = triangle_samples(panels.triangles[free], 16)
    bary, _ = triangle_samples(np.eye(3)[None], 16)
    shape = p2_shape(bary[0]) if degree == 2 else bary[0]
    mass = np.zeros_like(st)
    for face, area in zip(surface_faces, panels.areas[free]):
        mass[np.ix_(face, face)] += area * shape.T @ (weights[:, None] * shape)
    potential = trace(coords[fs])
    q = solve(mass, st @ potential, assume_a="pos")
    lam = np.array([0.0, np.tanh(4.0), np.sqrt(1.25) * np.tanh(4 * np.sqrt(1.25))])
    values = shape @ q[surface_faces]
    exact = (trace(samples.reshape(-1, 3)) * lam).reshape(values.shape)
    error = np.sqrt(
        np.einsum(
            "f,q,fqk->k",
            panels.areas[free],
            weights,
            (values[:, :, 1:] - exact[:, :, 1:]) ** 2,
        )
        / np.einsum("f,q,fqk->k", panels.areas[free], weights, exact[:, :, 1:] ** 2)
    )
    eig = eigh(st, mass, eigvals_only=True)
    volume = 32 * np.pi**2

    def fields(points):
        x, y, z = points.T
        value = np.column_stack(
            (x, y, z, x * x - y * y, x * y, x * z, y * z, z * z - x * x)
        )
        grad = np.zeros((len(points), 8, 3))
        grad[:, :3] = np.eye(3)
        grad[:, 3, 0] = 2 * x
        grad[:, 3, 1] = -2 * y
        grad[:, 4, 0] = y
        grad[:, 4, 1] = x
        grad[:, 5, 0] = z
        grad[:, 5, 2] = x
        grad[:, 6, 1] = z
        grad[:, 6, 2] = y
        grad[:, 7, 0] = -2 * x
        grad[:, 7, 2] = 2 * z
        return value, grad

    u, w = np.polynomial.legendre.leggauss(4)
    x, y, z = np.meshgrid(2 * np.pi * u, np.pi * u, 2 * (u - 1), indexing="ij")
    volume_weights = np.einsum("i,j,k->ijk", w, w, w).ravel() * (4 * np.pi**2)
    _, g = fields(np.column_stack((x.ravel(), y.ravel(), z.ravel())))
    expected = np.einsum("q,qik,qjk->ij", volume_weights, g, g)
    polynomial, _ = fields(coords)
    actual = polynomial.T @ energy @ polynomial
    patch_scale = np.sqrt(np.diag(expected)[:, None] * np.diag(expected)[None, :])
    result = dict(
        qualified=False,
        degree=degree,
        resolution=resolution,
        order=order,
        faces=len(mesh.faces),
        surface_dofs=len(fs),
        seconds=perf_counter() - start,
        **diagnostics,
        constant_q=float(np.max(np.abs(q[:, 0]))),
        net_flux=(np.ones(len(fs)) @ mass @ q).tolist(),
        continuous_q_relative_l2=error.tolist(),
        lowest_eigenvalues=eig[:5].tolist(),
        linear_energy_error=float(
            np.max(np.abs(coords.T @ energy @ coords - volume * np.eye(3))) / volume
        ),
        harmonic_degree_two_energy_matrix_error=float(
            np.max(np.abs(actual - expected) / patch_scale)
        ),
        harmonic_degree_two_diagonal_errors=(
            (np.diag(actual) - np.diag(expected)) / np.diag(expected)
        ).tolist(),
        peak_rss_bytes=peak_rss_bytes(),
    )
    return result, (
        coords[fs],
        surface_faces,
        mass,
        st,
        panels.areas[free],
        samples,
        weights,
        shape,
    )


def wave(resolution=(8, 4, 4), order=3, degree=2):
    start = perf_counter()
    check, data = run(resolution, order, degree)
    coords, surface_faces, mass, s, areas, points, weights, shape = data
    dtn = solve(mass, s, assume_a="pos")
    n = len(coords)
    g = 9.81
    a = 0.001
    initial = a * np.cos(coords[:, 0] + 2 * np.pi)
    omega = np.sqrt(g * np.tanh(4.0))
    times = np.linspace(0, 10 * np.pi / omega, 1001)

    def rhs(t, y):
        return np.r_[dtn @ y[n:], -g * y[:n]]

    solution = solve_ivp(
        rhs,
        (times[0], times[-1]),
        np.r_[initial, np.zeros(n)],
        method="RK23",
        rtol=1e-9,
        atol=1e-12,
        t_eval=times,
    )
    assert solution.success
    eta, psi = solution.y[:n], solution.y[n:]
    exact = a * np.cos(points[:, :, 0] + 2 * np.pi)
    c = np.einsum("f,q,fq->", areas, weights, exact**2)
    b = np.zeros(n)
    for j in range(shape.shape[1]):
        load = np.einsum("f,q,q,fq->f", areas, weights, shape[:, j], exact)
        np.add.at(b, surface_faces[:, j], load)
    cosine = np.cos(omega * times)
    error = np.sum(eta * (mass @ eta), axis=0) - 2 * cosine * (b @ eta) + cosine**2 * c
    continuous = float(
        np.sqrt(np.trapezoid(error, times) / np.trapezoid(cosine**2 * c, times))
    )
    energy = 0.5 * (
        g * np.sum(eta * (mass @ eta), axis=0) + np.sum(psi * (s @ psi), axis=0)
    )
    volume = np.ones(n) @ mass @ eta
    return dict(
        **check,
        periods=5,
        rtol=1e-9,
        rhs=solution.nfev,
        continuous_elevation_relative_l2_time=continuous,
        energy_drift=float(np.max(np.abs(energy / energy[0] - 1))),
        volume_drift=float(np.max(np.abs(volume - volume[0]))),
        total_seconds=perf_counter() - start,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fine",
        action="store_true",
        help="include the more expensive matched fine pair",
    )
    args = parser.parse_args()
    cases = [((16, 8, 8), 1), ((8, 4, 4), 2)]
    if args.fine:
        cases.extend([((32, 16, 16), 1), ((16, 8, 8), 2)])
    rows = []
    Path("artifacts").mkdir(exist_ok=True)
    with threadpool_limits(1):
        for resolution, degree in cases:
            row = wave(resolution, degree=degree)
            rows.append(row)
            Path("artifacts/free-surface-p2-comparison.json").write_text(
                json.dumps(rows, indent=2) + "\n", encoding="utf-8"
            )
            print(json.dumps(row), flush=True)
