"""Independent canonical-hull field accuracy report; exits 1 for unmet gates.

This command deliberately exposes discretization failures even when the algebraic
solver tests pass: uv run --locked python -m tests.verify_laplace
"""

import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from topwater import HullMesh
from topwater.bem import assemble_laplace
from topwater.boundary import BoundarySystem
from topwater.rankine import TrianglePanels, panel_influence
from .benchmark_rankine import benchmark_environment
from .laplace_reference import harmonic_field, point_source
from .reference_hull import mesh_arrays


def field(points, case):
    if case == 'harmonic':
        return harmonic_field(points)
    sources = dict(external_source=[3.5, 1., .7], internal_source=[0, 0, -.175])
    return point_source(points, sources[case])


def verify():
    rows = []
    with threadpool_limits(1):
        for bounded in (True, False):
            for refinement in (2, 3, 4):
                v, f = mesh_arrays(refinement)
                panels = TrianglePanels.from_mesh(HullMesh(v, f), 1 if bounded else -1)
                operators = assemble_laplace(panels)
                points = np.array([[-.5, 0, -.15], [0, .1, -.15], [.5, -.1, -.1]] if bounded else
                                  [[-2.5, .2, .3], [2, .8, -.1], [0, 1, .5]])
                influence = panel_influence(points, panels, gradients=True)
                for case in (('harmonic', 'external_source') if bounded else ('internal_source',)):
                    phi, gradient = field(panels.centroids, case)
                    q = np.einsum('nrj,nj->nr', gradient, panels.normals)
                    bulk_phi, bulk_gradient = field(points, case)
                    for kind in ('dirichlet', 'neumann', 'mixed'):
                        # Centroid point-source flux is not exactly compatible.
                        # Do not repair it or relax the Neumann compatibility gate.
                        if case == 'external_source' and kind == 'neumann':
                            continue
                        mask = (np.ones(len(panels), bool) if kind == 'dirichlet' else
                                np.zeros(len(panels), bool) if kind == 'neumann' else
                                panels.centroids[:, 0] > 0)
                        system = BoundarySystem(operators, dirichlet=mask, areas=panels.areas,
                                                length_scale=3.5, bounded=bounded)
                        solution = system.solve(potential=phi, normal_velocity=q)
                        expected_phi, expected_bulk_phi = phi.copy(), bulk_phi.copy()
                        if bounded and kind == 'neumann':
                            mean = panels.areas @ phi/panels.areas.sum()
                            expected_phi -= mean
                            expected_bulk_phi -= mean
                        evaluated_phi = (influence.single_layer @ solution.normal_velocity
                                         - influence.double_layer @ solution.potential)
                        evaluated_gradient = (
                            np.einsum('mnj,nr->mrj', influence.gradient_single_layer,
                                      solution.normal_velocity)
                            - np.einsum('mnj,nr->mrj', influence.gradient_double_layer,
                                        solution.potential))
                        def rms(values):
                            return np.sqrt(panels.areas @ (values**2)/panels.areas.sum())
                        ep = rms(solution.potential-expected_phi)
                        eq = 3.5*rms(solution.normal_velocity-q)
                        eu = 3.5*np.sqrt(np.mean(np.sum((evaluated_gradient-bulk_gradient)**2, axis=2), axis=0))
                        rows.append(dict(panel_count=len(panels), bounded=bounded, case=case,
                                         boundary=kind, potential_scale=1., length_scale=3.5,
                                         potential_error=ep.tolist(), normal_velocity_error=eq.tolist(),
                                         bulk_velocity_error=eu.tolist(),
                                         bulk_potential_error=np.sqrt(np.mean((evaluated_phi-expected_bulk_phi)**2,
                                                                              axis=0)).tolist(),
                                         gate_passed=bool(max(np.max(ep), np.max(eq), np.max(eu)) <= .01),
                                         solve_diagnostics=solution.diagnostics))
    final = [row for row in rows if row['panel_count'] == 2048]
    convergence = []
    for case, boundary in sorted({(row['case'], row['boundary']) for row in rows}):
        series = [row for row in rows if row['case'] == case and row['boundary'] == boundary]
        metrics = {}
        for name in ('potential_error', 'normal_velocity_error', 'bulk_velocity_error'):
            values = np.array([row[name] for row in series])
            # Prescribed fields are exact; do not demand a decrease below roundoff.
            metrics[name] = (np.all(np.diff(values, axis=0) < 0, axis=0)
                             | np.all(values < 1e-12, axis=0)).tolist()
        convergence.append(dict(case=case, boundary=boundary, decreasing=metrics))
    missing = ['bounded point-source pure Neumann with panel-average flux',
               'surface velocity recovery']
    decreasing = all(all(all(flags) for flags in row['decreasing'].values()) for row in convergence)
    magnitude_passed = all(row['gate_passed'] for row in final)
    return dict(benchmark_environment(), results=rows,
                finest_gate_passed=magnitude_passed, refinement_convergence=convergence,
                missing_cases=missing, stage_gate_passed=magnitude_passed and decreasing and not missing)


if __name__ == '__main__':
    report = verify()
    path = Path('artifacts/laplace-verification.json')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+"\n", encoding='utf-8')
    for row in report['results']:
        if row['panel_count'] == 2048:
            print(row['case'], row['boundary'], 'gate', row['gate_passed'],
                  'phi', row['potential_error'], 'q', row['normal_velocity_error'],
                  'velocity', row['bulk_velocity_error'])
    raise SystemExit(0 if report['stage_gate_passed'] else 1)
