"""Independent plane-wave checks for exact native packet forcing and source reuse."""
from types import SimpleNamespace
from dataclasses import replace

import capytaine as cpt
import numpy as np
import pytest
from numpy.testing import assert_allclose
from capytaine.bem.airy_waves import airy_waves_pressure, airy_waves_velocity, froude_krylov_force

from topwater.solve.preparation import prepare_case
from topwater.inputs import WaveComponents
from topwater.solve.diffraction import packet_diffraction
from topwater.solve.hydrodynamics import prepare_hydrodynamics
from topwater.solve.settings import SolveOptions
from topwater.solve.types import MODES
from topwater.solve.surface import combined_diffraction
from .test_workflow_data import boat


@pytest.fixture
def packet_case(tmp_path):
    setup = prepare_case(boat().to_body(), WaveComponents([], [], [], []), 1., rest_s=1., ramp_s=1.)
    random = np.random.default_rng(71)
    coefficients = random.normal(size=(3, 2, 4))+1j*random.normal(size=(3, 2, 4))
    coefficients[:, :, 0] = 0.
    setup = replace(setup, train=SimpleNamespace(omega=np.array([0., .37, 1.13, 2.41]),
        theta=np.array([.123, .789, 1.2]), C=coefficients))
    options = SolveOptions(coefficient_step_rad_s=.5, radiation_bandwidth_rad_s=1.)
    hydro = prepare_hydrodynamics(setup, tmp_path/'hydrodynamics', options)
    return setup, hydro


def test_native_preparation_omits_reference_exports_and_diffraction_grid(packet_case):
    setup, hydro = packet_case
    assert list(hydro.directory.iterdir()) == []
    assert not hasattr(hydro, 'excitation_headings')
    assert not hasattr(hydro, 'excitation_omega')
    assert not hasattr(hydro, 'netcdf_path')
    assert len(hydro.radiation_results) == 6*(len(hydro.radiation_omega)+1)
    assert all(mode in MODES for _, mode in hydro.radiation_results)


@pytest.mark.parametrize('index', [1, 2, 3])
def test_coherent_boundary_pressure_and_force_match_independent_plane_solves(packet_case, index):
    setup, hydro = packet_case
    response = packet_diffraction(hydro, setup.train)
    points = hydro.body.mesh_including_lid.faces_centers
    physical = hydro.body.mesh.faces_centers
    velocity = np.zeros((len(points), 3), complex)
    pressure = np.zeros(len(physical), complex)
    force = np.zeros(6, complex)
    for j, theta in enumerate(setup.train.theta):
        for branch in (0, 1):
            problem = cpt.DiffractionProblem(body=hydro.body, rho=hydro.density,
                g=hydro.gravity, water_depth=hydro.depth, omega=setup.train.omega[index],
                wave_direction=float((theta+branch*np.pi) % (2*np.pi)))
            result = hydro.solver.solve(problem, keep_details=True)
            coefficient = setup.train.C[j, branch, index]
            velocity += coefficient*airy_waves_velocity(points, problem)
            pressure += coefficient*airy_waves_pressure(physical, problem)
            fk = froude_krylov_force(problem)
            force += coefficient*np.array([result.forces[mode]+fk[mode] for mode in MODES])
    assert_allclose(response.incident(index, points, velocity=True), velocity, rtol=1e-12, atol=1e-12)
    assert_allclose(response.incident(index, physical), pressure, rtol=1e-12, atol=1e-9)
    assert_allclose(response.force(index), force, rtol=1e-10, atol=1e-8)


def test_reconstruction_reuses_motion_sources_without_recomputing_forces(packet_case, monkeypatch):
    setup, hydro = packet_case
    response = packet_diffraction(hydro, setup.train)
    for index in range(1, len(setup.train.omega)):
        response.force(index)
    monkeypatch.setattr(hydro.solver, 'solve', lambda *a, **kw: pytest.fail('repeated BEM solve'))
    monkeypatch.setattr(response, 'force', lambda *a: pytest.fail('reconstruction computed force'))
    backend = SimpleNamespace(free_surface_elevation=lambda xy, result: np.full(len(xy), result.omega))
    _, scattered = combined_diffraction(hydro, setup.train, [[4., 4.]],
        return_incident=False, field_backend=backend)
    assert_allclose(scattered[:, 0], setup.train.omega)


def test_calm_packet_has_zero_force_without_diffraction_solves(packet_case, monkeypatch):
    setup, hydro = packet_case
    setup.train.C.fill(0.)
    response = packet_diffraction(hydro, setup.train)
    monkeypatch.setattr(hydro.solver, 'solve', lambda *a, **kw: pytest.fail('calm BEM solve'))
    for index in range(1, len(setup.train.omega)):
        assert_allclose(response.force(index), 0., atol=0.)


def test_reference_backend_still_exports_its_original_tables(packet_case, tmp_path):
    from topwater.solve.wecsim import exchange as producer
    setup, _ = packet_case
    h = setup.hydrostatic
    hydro = producer.build_hydrodynamics(setup.wet_triangles_m, setup.mean_cg_m,
        setup.mass_matrix, tmp_path/'reference', hydrostatic_stiffness=h.stiffness,
        cb=h.center_of_buoyancy, volume=h.displaced_volume,
        radiation_omega=np.array([.1, .2, .3]), excitation_omega=np.array([.1, .2]),
        headings=np.array([0., np.pi/4, np.pi/2, 2*np.pi]), field_omega=np.array([.1, .2]),
        green_cache_dir=tmp_path/'green')
    assert hydro.netcdf_path.is_file()
    assert hydro.native_mat_path.is_file()
    assert hydro.excitation.shape == (2, 4, 6)
    assert len(hydro.diffraction_results) == 9
