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
from .laplace_reference import harmonic_field, point_source, triangle_samples
from .analytic_geometry import sphere_mesh


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
                v, f = sphere_mesh(refinement)
                panels = TrianglePanels.from_mesh(HullMesh(v, f), 1 if bounded else -1)
                operators = assemble_laplace(panels)
                points = np.array([[-.5, 0, -.15], [0, .1, -.15], [.5, -.1, -.1]] if bounded else
                                  [[-2.5, .2, .3], [2, .8, -.1], [0, 1, .5]])
                influence = panel_influence(points, panels, gradients=True)
                surface_points, surface_weights = triangle_samples(panels.triangles)
                for case in (('harmonic', 'external_source') if bounded else ('internal_source',)):
                    phi, gradient = field(panels.centroids, case)
                    q = np.einsum('nrj,nj->nr', gradient, panels.normals)
                    bulk_phi, bulk_gradient = field(points, case)
                    surface_phi, surface_gradient = field(surface_points.reshape(-1, 3), case)
                    surface_phi = surface_phi.reshape(len(panels), len(surface_weights), -1)
                    surface_q = np.einsum('fqrj,fj->fqr',
                        surface_gradient.reshape(len(panels), len(surface_weights), -1, 3),
                        panels.normals)
                    average_q = np.einsum('q,fqr->fr', surface_weights, surface_q)
                    def surface_rms(values):
                        return np.sqrt(np.einsum('f,q,fqr->r', panels.areas,
                            surface_weights, values**2)/panels.areas.sum())
                    for kind in ('dirichlet', 'neumann', 'mixed'):
                        mask = (np.ones(len(panels), bool) if kind == 'dirichlet' else
                                np.zeros(len(panels), bool) if kind == 'neumann' else
                                panels.centroids[:, 0] > 0)
                        system = BoundarySystem(operators, dirichlet=mask, areas=panels.areas,
                                                length_scale=3.5, bounded=bounded)
                        # Integrate the actual prescribed flux, without projecting
                        # away its net value or weakening the compatibility gate.
                        known_q = average_q if bounded and kind == 'neumann' else q
                        solution = system.solve(potential=phi, normal_velocity=known_q)
                        expected_phi, expected_bulk_phi = phi.copy(), bulk_phi.copy()
                        expected_surface_phi = surface_phi.copy()
                        if bounded and kind == 'neumann':
                            mean = np.einsum('f,q,fqr->r', panels.areas, surface_weights,
                                             surface_phi)/panels.areas.sum()
                            expected_phi -= mean
                            expected_bulk_phi -= mean
                            expected_surface_phi -= mean
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
                        surface_ep = surface_rms(solution.potential[:, None]-expected_surface_phi)
                        surface_eq = surface_rms(solution.normal_velocity[:, None]-surface_q)
                        relative_surface_ep = surface_ep/surface_rms(expected_surface_phi)
                        relative_surface_eq = surface_eq/surface_rms(surface_q)
                        relative_eu = (eu/3.5
                            /np.sqrt(np.mean(np.sum(bulk_gradient**2, axis=2), axis=0)))
                        rows.append(dict(panel_count=len(panels), bounded=bounded, case=case,
                                         boundary=kind, potential_scale=1., length_scale=3.5,
                                         potential_error=ep.tolist(), normal_velocity_error=eq.tolist(),
                                         bulk_velocity_error=eu.tolist(),
                                         surface_potential_error=surface_ep.tolist(),
                                         surface_normal_velocity_error=(3.5*surface_eq).tolist(),
                                         relative_surface_potential_error=relative_surface_ep.tolist(),
                                         relative_surface_normal_velocity_error=relative_surface_eq.tolist(),
                                         relative_bulk_velocity_error=relative_eu.tolist(),
                                         reference_net_flux=(panels.areas@average_q).tolist(),
                                         bulk_potential_error=np.sqrt(np.mean((evaluated_phi-expected_bulk_phi)**2,
                                                                              axis=0)).tolist(),
                                         gate_passed=bool(max(np.max(ep), np.max(eq), np.max(eu),
                                             np.max(surface_ep), 3.5*np.max(surface_eq)) <= .01),
                                         solve_diagnostics=solution.diagnostics))
    final = [row for row in rows if row['panel_count'] == 2048]
    convergence = []
    for case, boundary in sorted({(row['case'], row['boundary']) for row in rows}):
        series = [row for row in rows if row['case'] == case and row['boundary'] == boundary]
        metrics = {}
        for name in ('potential_error', 'normal_velocity_error', 'bulk_velocity_error',
                     'surface_potential_error', 'surface_normal_velocity_error',
                     'original_bie_residual'):
            values = np.array([row['solve_diagnostics'][name] if name == 'original_bie_residual'
                               else row[name] for row in series])
            # Prescribed fields are exact; do not demand a decrease below roundoff.
            metrics[name] = (np.all(np.diff(values, axis=0) < 0, axis=0)
                             | np.all(values < 1e-12, axis=0)).tolist()
        convergence.append(dict(case=case, boundary=boundary, decreasing=metrics))
    missing = ['surface velocity recovery']
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
