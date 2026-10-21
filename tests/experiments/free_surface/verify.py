"""Independent flat-tank checks, not nonlinear/open-water qualification."""
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.linalg import solve, eigh
from threadpoolctl import threadpool_limits
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from tests.laplace_reference import triangle_samples
from .domain import tank, trace
from .operator import assemble_energy


def run(resolution, order=3):
    start = perf_counter()
    mesh, free, node_map = tank(*resolution)
    panels = TrianglePanels.from_mesh(mesh)
    energy, diagnostics = assemble_energy(mesh.faces, panels, order)
    fs = node_map.ravel()
    other = np.setdiff1d(np.arange(len(mesh.vertices)), fs)
    elimination = solve(energy[np.ix_(other, other)],
                        energy[np.ix_(other, fs)], assume_a='pos')
    steklov = (energy[np.ix_(fs, fs)]
               - energy[np.ix_(fs, other)]@elimination)
    lookup = np.full(len(mesh.vertices), -1, dtype=int)
    lookup[fs] = np.arange(len(fs))
    mass = np.zeros((len(fs), len(fs)))
    for face, area in zip(lookup[mesh.faces[free]], panels.areas[free]):
        mass[np.ix_(face, face)] += area*(np.ones((3, 3))+np.eye(3))/12
    potential = trace(mesh.vertices[fs])
    q = solve(mass, steklov@potential, assume_a='pos')
    expected_eigenvalues = np.array([
        0., np.tanh(4.), np.sqrt(1.25)*np.tanh(4*np.sqrt(1.25)),
    ])
    exact = potential*expected_eigenvalues
    interior = ((np.abs(mesh.vertices[fs, 0]) <= np.pi+1e-12)
                & (np.abs(mesh.vertices[fs, 1]) <= np.pi/2+1e-12))
    error = np.sqrt(np.sum((q[interior, 1:]-exact[interior, 1:])**2, axis=0)
                    / np.sum(exact[interior, 1:]**2, axis=0))
    samples, weights = triangle_samples(panels.triangles[free], 16)
    shape, _ = triangle_samples(np.eye(3)[None], 16)
    values = np.einsum('qc,fcr->fqr', shape[0], q[lookup[mesh.faces[free]]])
    surface_exact = (trace(samples.reshape(-1, 3))*expected_eigenvalues).reshape(values.shape)
    surface_error = np.sqrt(np.einsum(
        'f,q,fqr->r', panels.areas[free], weights, (values[:, :, 1:]-surface_exact[:, :, 1:])**2,
    ) / np.einsum('f,q,fqr->r', panels.areas[free], weights, surface_exact[:, :, 1:]**2))
    spectrum = eigh(steklov, mass, eigvals_only=True)
    expected_linear_energy = np.eye(3)*(4*np.pi*2*np.pi*4)
    linear = mesh.vertices
    patch_error = float(np.max(np.abs(
        linear.T@energy@linear-expected_linear_energy,
    ))/expected_linear_energy[0, 0])
    # Counterexample to the earlier vertex-average -> P0 trace proposal.
    # A nonconstant three-color nodal field disappears on every top triangle.
    i, j = np.indices(node_map.shape)
    pattern = np.zeros(len(mesh.vertices))
    pattern[node_map] = np.array([1., -.5, -.5])[(i+j) % 3]
    lost_trace = pattern[mesh.faces[free]].mean(axis=1)
    return dict(
        resolution=resolution, panels=len(panels), nodes=len(mesh.vertices),
        order=order, seconds=perf_counter()-start,
        peak_rss_bytes=peak_rss_bytes(), **diagnostics,
        linear_energy_relative_error=patch_error,
        constant_flux_max=float(np.max(np.abs(q[:, 0]))),
        interior_relative_rms=error.tolist(),
        continuous_surface_velocity_relative_rms=surface_error.tolist(),
        modal_eigenvalues=(np.sum(potential*(steklov@potential), axis=0)
                           / np.sum(potential*(mass@potential), axis=0)).tolist(),
        expected_eigenvalues=expected_eigenvalues.tolist(),
        lowest_eigenvalues=spectrum[:5].tolist(),
        net_top_flux=(np.ones(len(fs))@mass@q).tolist(),
        rejected_p0_transfer_nonconstant_trace_max=float(np.max(np.abs(lost_trace))),
    )


if __name__ == '__main__':
    with threadpool_limits(1):
        rows = [run(resolution) for resolution in ((8, 4, 4), (16, 8, 8))]
    report = dict(results=rows, qualified=False, missing=[
        'tolerance-controlled singular/near/far Galerkin integration',
        'nonlinear energy and long-duration space/time convergence',
        'open/deep-water domain and radiation treatment',
        'wave-body coupling and finite-frequency Capytaine correlation',
    ])
    Path('artifacts').mkdir(exist_ok=True)
    Path('artifacts/free-surface-verification.json').write_text(
        json.dumps(report, indent=2)+'\n', encoding='utf-8',
    )
    print(json.dumps(report, indent=2))
