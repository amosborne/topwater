from dataclasses import FrozenInstanceError

import numpy as np
from numpy.testing import assert_allclose
import pytest
from threadpoolctl import threadpool_info, threadpool_limits

from topwater import HullMesh
from topwater.bem import LaplaceOperators, assemble_laplace
from topwater.boundary import BoundarySystem
from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels, panel_influence
from .laplace_reference import harmonic_field, point_source, triangle_samples
from .reference_hull import mesh_arrays


@pytest.fixture(scope="module")
def canonical():
    vertices, faces = mesh_arrays()
    panels = TrianglePanels.from_mesh(HullMesh(vertices, faces))
    return panels, assemble_laplace(panels)


def harmonic_fields(panels):
    potential, gradient = harmonic_field(panels.centroids)
    flux = np.einsum("nrj,nj->nr", gradient, panels.normals)
    return potential, flux


@pytest.mark.parametrize("kind", ["dirichlet", "neumann", "mixed"])
def test_reused_factors_batch_fields_scaling_and_independent_harmonic_convergence(kind):
    errors = []
    with threadpool_limits(1):
        for refinement in (2, 3, 4):
            v, f = mesh_arrays(refinement)
            p = TrianglePanels.from_mesh(HullMesh(v, f))
            operators = assemble_laplace(p)
            phi, q = harmonic_fields(p)
            mask = (np.ones(len(p), bool) if kind == "dirichlet" else
                    np.zeros(len(p), bool) if kind == "neumann" else p.centroids[:, 0] > 0)
            system = BoundarySystem(operators, dirichlet=mask, areas=p.areas,
                                    length_scale=3.5, bounded=True)
            solution = system.solve(potential=phi, normal_velocity=q)
            expected = phi.copy()
            if kind == "neumann":
                expected -= p.areas @ phi/p.areas.sum()
                assert_allclose(solution.diagnostics['gauge_residual'], 0, atol=1e-14)
                assert max(solution.diagnostics['compatibility_residual']) < 1e-12
            phi_error = np.sqrt(np.sum(p.areas[:, None]*(solution.potential-expected)**2,
                                      axis=0)/p.areas.sum())
            q_error = 3.5*np.sqrt(np.sum(p.areas[:, None]*(solution.normal_velocity-q)**2,
                                        axis=0)/p.areas.sum())
            errors.append(np.maximum(phi_error, q_error))
            assert max(solution.diagnostics['backward_error']) < 1e-12
            # Repeated solves reuse precisely the same factors and match batch columns.
            for column in range(2):
                separate = system.solve(potential=phi[:, column], normal_velocity=q[:, column])
                assert_allclose(separate.potential, solution.potential[:, column], atol=1e-12)
                assert_allclose(separate.normal_velocity, solution.normal_velocity[:, column], atol=1e-12)
    assert np.all(np.diff(errors, axis=0) < 0)
    if kind == "neumann":
        assert np.max(errors[-1]) < .01
    # Dirichlet/mixed flux converges but does NOT yet meet the 1% stage gate.
    # The verification report exposes that limitation; do not certify Stage 1.


def test_configuration_ownership_and_physical_length_scale_invariance(canonical):
    panels, operators = canonical
    phi, q = harmonic_fields(panels)
    mask = panels.centroids[:, 0] > 0
    solutions = []
    for length in (.035, 3.5, 350.):
        system = BoundarySystem(operators, dirichlet=mask, areas=panels.areas,
                                length_scale=length, bounded=True)
        with pytest.raises(FrozenInstanceError):
            system.length_scale = 1
        solutions.append(system.solve(potential=phi, normal_velocity=q))
    for solution in solutions[1:]:
        assert_allclose(solution.potential, solutions[0].potential, atol=1e-12)
        assert_allclose(solution.normal_velocity, solutions[0].normal_velocity, atol=1e-12)
        assert not solution.potential.flags.writeable
        assert not solution.normal_velocity.flags.writeable


def test_neumann_compatibility_zero_problem_and_validation(canonical):
    panels, operators = canonical
    system = BoundarySystem(operators, dirichlet=np.zeros(len(panels), bool),
                            areas=panels.areas, length_scale=3.5, bounded=True)
    with pytest.raises(ValueError, match="net flux"):
        system.solve(normal_velocity=np.ones(len(panels)))
    with pytest.raises(ValueError, match="missing prescribed"):
        system.solve(potential=np.zeros(len(panels)))
    with pytest.raises(ValueError, match="finite"):
        system.solve(normal_velocity=np.full(len(panels), np.nan))
    solution = system.solve(normal_velocity=np.zeros(len(panels)))
    assert_allclose(solution.potential, 0, atol=0)
    assert solution.diagnostics['backward_error'] == [0]
    assert solution.diagnostics['original_bie_residual'] == [0]


def test_external_source_neumann_uses_converged_flux_integrals_without_projection(canonical):
    panels, operators = canonical
    source = [3.5, 1., .7]
    averages = []
    for order in (8, 12):
        points, weights = triangle_samples(panels.triangles, order)
        _, gradient = point_source(points.reshape(-1, 3), source)
        q = np.einsum('fqrj,fj->fqr', gradient.reshape(len(panels), len(weights), -1, 3),
                      panels.normals)
        averages.append(np.einsum('q,fqr->fr', weights, q))
    assert_allclose(averages[0], averages[1], atol=1e-14, rtol=1e-12)
    assert_allclose(panels.areas@averages[1], 0, atol=1e-14)
    _, gradient = point_source(panels.centroids, source)
    centroid_q = np.einsum('frj,fj->fr', gradient, panels.normals)
    assert abs((panels.areas@centroid_q)[0]) > 1e-7
    system = BoundarySystem(operators, dirichlet=np.zeros(len(panels), bool),
                            areas=panels.areas, length_scale=3.5, bounded=True)
    result = system.solve(normal_velocity=averages[1])
    assert max(result.diagnostics['compatibility_residual']) < 1e-12
    assert max(result.diagnostics['backward_error']) < 1e-12
    # The bordered solve is accurate. Its original BIE defect remains a
    # discretization error, reported separately rather than called roundoff.
    assert np.all(np.isfinite(result.diagnostics['original_bie_residual']))


def test_surface_quadrature_exposes_constant_panel_interpolation_error(canonical):
    panels, _ = canonical
    points, weights = triangle_samples(panels.triangles)
    phi, _ = harmonic_field(panels.centroids)
    exact, _ = harmonic_field(points.reshape(-1, 3))
    error = phi[:, None]-exact.reshape(len(panels), len(weights), -1)
    surface_mse = np.einsum('q,fqr->fr', weights, error**2)
    # A linear field's exact triangle variance is the mean squared vertex
    # deviation divided by four; its value at the centroid alone has zero error.
    vertex_deviation = (panels.triangles[:, :, 0]-panels.centroids[:, None, 0])/3.5
    expected = np.sum(vertex_deviation**2, axis=1)/12
    assert_allclose(surface_mse[:, 0], expected, atol=1e-16, rtol=1e-12)
    assert np.all(surface_mse[:, 0] > 0)


def test_surface_error_quadrature_resolves_the_near_internal_source(canonical):
    panels, _ = canonical
    source = [0, 0, -.175]
    phi, gradient = point_source(panels.centroids, source)
    q = np.einsum('frj,fj->fr', gradient, panels.normals)
    errors = []
    for order in (16, 24):
        points, weights = triangle_samples(panels.triangles, order)
        exact, gradient = point_source(points.reshape(-1, 3), source)
        exact_q = np.einsum('fqrj,fj->fqr',
            gradient.reshape(len(panels), len(weights), -1, 3), panels.normals)
        metrics = []
        for value, reference in ((phi, exact.reshape(exact_q.shape)), (q, exact_q)):
            metrics.append(np.sqrt(np.einsum('f,q,fqr->r', panels.areas, weights,
                (value[:, None]-reference)**2)/panels.areas.sum()))
        errors.append(metrics)
    assert_allclose(errors[0], errors[1], atol=1e-8, rtol=1e-9)


def test_factorization_limits_are_scoped_and_explicit(canonical):
    panels, operators = canonical
    with threadpool_limits(2, user_api='blas'):
        system = BoundarySystem(operators, dirichlet=np.ones(len(panels), bool),
                                areas=panels.areas, length_scale=3.5, bounded=True)
        assert all(pool['num_threads'] == 1 for pool in system.diagnostics['factorization_blas'])
        assert all(pool['num_threads'] == 2 for pool in threadpool_info() if pool['user_api'] == 'blas')
        ambient = BoundarySystem(operators, dirichlet=np.ones(len(panels), bool),
                                 areas=panels.areas, length_scale=3.5, bounded=True,
                                 factorization_threads=None)
        assert all(pool['num_threads'] == 2 for pool in ambient.diagnostics['factorization_blas'])


@pytest.mark.parametrize("small", [0, 1e-14])
def test_singular_and_ill_conditioned_systems_are_rejected(small):
    single = np.diag([1., small])
    double = np.zeros((2, 2))
    single.setflags(write=False)
    double.setflags(write=False)
    operators = LaplaceOperators(single, double, {})
    with pytest.raises(np.linalg.LinAlgError):
        BoundarySystem(operators, dirichlet=np.ones(2, bool), areas=np.ones(2),
                       length_scale=1, bounded=True)


def test_rigid_transform_and_uniform_length_change_preserve_neumann_solution(canonical):
    panels, operators = canonical
    _, q = harmonic_fields(panels)
    baseline = BoundarySystem(operators, dirichlet=np.zeros(len(panels), bool),
                              areas=panels.areas, length_scale=3.5, bounded=True).solve(normal_velocity=q)
    rotation = rotation_matrix(np.array([1, 2, 3, 4])/np.sqrt(30))
    scale = 2.5
    moved = TrianglePanels(scale*panels.triangles@rotation.T+[3, -4, 2])
    transformed = BoundarySystem(assemble_laplace(moved), dirichlet=np.zeros(len(moved), bool),
                                 areas=moved.areas, length_scale=3.5*scale,
                                 bounded=True).solve(normal_velocity=q/scale)
    assert_allclose(transformed.potential, baseline.potential, rtol=1e-10, atol=1e-11)


def test_decaying_exterior_source_allows_net_flux_and_converges_in_bulk_velocity():
    points = np.array([[-2.5, .2, .3], [2, .8, -.1], [0, 1, .5]])
    source = np.array([0, 0, -.175])
    _, expected_velocity = point_source(points, source)
    errors = []
    flux_errors = []
    with threadpool_limits(1):
        for refinement in (2, 3, 4):
            v, f = mesh_arrays(refinement)
            panels = TrianglePanels.from_mesh(HullMesh(v, f), -1)
            _, gradient = point_source(panels.centroids, source)
            q = np.einsum('nrj,nj->nr', gradient, panels.normals)
            system = BoundarySystem(assemble_laplace(panels), dirichlet=np.zeros(len(panels), bool),
                                    areas=panels.areas, length_scale=3.5, bounded=False)
            solution = system.solve(normal_velocity=q)
            net_flux = solution.diagnostics['net_flux'][0]
            assert net_flux > 0
            flux_errors.append(abs(net_flux-4*np.pi))
            assert solution.diagnostics['compatibility_residual'] is None
            influence = panel_influence(points, panels, gradients=True)
            velocity = (np.einsum('mnj,nr->mrj', influence.gradient_single_layer, solution.normal_velocity)
                        - np.einsum('mnj,nr->mrj', influence.gradient_double_layer, solution.potential))
            errors.append(3.5*np.sqrt(np.mean(np.sum((velocity-expected_velocity)**2, axis=2))))
    assert np.all(np.diff(errors) < 0)
    assert np.all(np.diff(flux_errors) < 0)
    assert flux_errors[-1]/(4*np.pi) < .01
    assert errors[-1] < .01
