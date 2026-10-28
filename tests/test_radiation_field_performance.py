"""The optimized field path must match Capytaine's public evaluator."""
from types import SimpleNamespace

import numpy as np
import pytest


def test_shared_source_mesh_preserves_native_field_matrices_and_complex_transfers(monkeypatch):
    cpt = pytest.importorskip('capytaine')
    from capytaine.green_functions.delhommeau import Delhommeau
    from threadpoolctl import threadpool_info, threadpool_limits
    from topwater._capytaine_surface import MODES, radiation_velocity_transfer
    from topwater._wecsim_setup import floating_equilibrium
    from .reference_boat import reference_body

    body = reference_body()
    state, h = floating_equilibrium(body)
    points, ids = np.unique(h.wetted_triangles.reshape(-1, 3), axis=0, return_inverse=True)
    mesh = cpt.Mesh(vertices=points, faces=ids.reshape(-1, 3), auto_clean=False, auto_check=False)
    floating = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=state.position),
                                center_of_mass=state.position)
    # Both paths use the very same evaluator/settings. No large tabulation is
    # needed to test this identity, and no optional MATLAB engine is involved.
    solver = cpt.BEMSolver(green_function=Delhommeau(tabulation_nr=0, tabulation_nz=0,
                                                   tabulation_cache_dir=None))
    omega = np.array([.8, 2., 8.])
    xy = np.array([[-6., -3.], [0., 3.], [6., 3.]])
    results = {}
    for w in omega:
        for mode in MODES:
            results[(w, mode)] = solver.solve(cpt.RadiationProblem(body=floating, omega=w,
                rho=1025., g=9.81, water_depth=np.inf, radiating_dof=mode),
                keep_details=True, n_threads=1)
    hydro = SimpleNamespace(body=floating, solver=solver, radiation_results=results,
                            radiation_omega=omega, gravity=9.81, density=1025., depth=np.inf)
    field = np.c_[xy, np.zeros(len(xy))]
    source = floating.mesh_including_lid
    expected = []
    # Production already used one thread around the complete reconstruction.
    with threadpool_limits(1):
        for w in omega:
            result = results[(w, MODES[0])]
            matrix = solver.engine.build_S_matrix(field, source, free_surface=result.free_surface,
                water_depth=result.water_depth, wavenumber=result.encounter_wavenumber)
            sources = np.column_stack([results[(w, mode)].sources for mode in MODES])
            expected.append(-(matrix @ sources)/9.81)
    calls = []
    merged = type(source).merged
    def count_merge(self, *args, **kwargs):
        calls.append(True)
        return merged(self, *args, **kwargs)
    monkeypatch.setattr(type(source), 'merged', count_merge)
    before_threads = [(p['prefix'], p['num_threads']) for p in threadpool_info()]
    solve = solver.solve
    thread_checks = []
    def controlled_solve(*args, **kwargs):
        assert kwargs.get('n_threads') is None  # No repeated discovery inside Capytaine.
        assert all(p['num_threads'] == 1 for p in threadpool_info())
        thread_checks.append(True)
        return solve(*args, **kwargs)
    monkeypatch.setattr(solver, 'solve', controlled_solve)
    del results[(omega[-1], MODES[-1])]
    actual = radiation_velocity_transfer(hydro, xy, shared_green=True)
    np.testing.assert_array_equal(actual, expected)
    assert len(calls) == 1  # One immutable source mesh for every frequency/mode.
    assert len(thread_checks) == 1
    assert [(p['prefix'], p['num_threads']) for p in threadpool_info()] == before_threads
    public = radiation_velocity_transfer(hydro, xy, shared_green=False)
    np.testing.assert_allclose(actual, public, rtol=2e-13, atol=2e-14)
