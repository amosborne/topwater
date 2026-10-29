"""Shared surface adapters and public external-model entry point."""
from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose

from topwater.solve import simulate
from topwater.solve import SolveOptions
from topwater.solve import workflow
from topwater.solve.native import backend as motion_backend
from topwater.solve.native.dynamics import Motion
from .test_workflow_data import boat, result
from topwater.workflow_data import Sea


@pytest.mark.parametrize('failure', [None, 'close', 'reconstruct', 'both'])
def test_native_public_workflow_preserves_surface_adapter_and_world_offset(monkeypatch, tmp_path, failure):
    from topwater.solve import surface as producer
    expected = result()
    offset = .03
    setup = SimpleNamespace(wet_triangles_m=np.zeros((1, 3, 3)), mean_cg_m=np.zeros(3),
        mass_matrix=np.eye(6), world_vertical_offset_m=offset,
        hydrostatic=SimpleNamespace(stiffness=np.eye(6), center_of_buoyancy=np.zeros(3), displaced_volume=1.))
    setup.train = SimpleNamespace(omega=np.array([0., 1.]), certify_incident_startup=lambda _: None)
    motion = Motion(expected.time_s*.1, np.zeros((len(expected.time_s), 6)),
        np.zeros((len(expected.time_s), 6)), np.zeros((len(expected.time_s), 6)),
        expected.position_m, expected.orientation_wxyz, np.zeros((len(expected.time_s), 6)),
        np.zeros((len(expected.time_s), 0)), {})
    influence = object()
    options = SolveOptions(time_step_s=.05, packet_tolerance_m=.005,
        cpu_workers=2, batch_points=8)
    def prepare(*args, **kwargs):
        assert kwargs['time_step_s'] == .05
        assert kwargs['envelope_tolerance_m'] == .005
        return setup
    monkeypatch.setattr(workflow, 'prepare_case', prepare)
    closed = []
    def close():
        closed.append(True)
        if failure in ('close', 'both'):
            raise RuntimeError('injected close failure')
    backend = SimpleNamespace(close=close)
    monkeypatch.setattr(workflow, 'select_field_backend', lambda *args, **kw: backend)
    monkeypatch.setattr(workflow.tempfile, 'mkdtemp', lambda **kw: str(tmp_path/'scratch'))
    (tmp_path/'scratch').mkdir()
    hydro = object()
    monkeypatch.setattr(motion_backend, 'prepare_hydrodynamics', lambda *args, **kw: hydro)
    excitation = object()
    monkeypatch.setattr(motion_backend.Excitation, 'from_case', lambda *args: excitation)
    def integrate(case, coefficients, *, influence, **kwargs):
        assert case is setup and coefficients is hydro
        assert influence is supplied
        assert kwargs == dict(dt_s=.05, excitation=excitation)
        return motion
    supplied = influence
    monkeypatch.setattr(motion_backend, 'integrate_motion', integrate)
    def reconstruct(coefficients, case, computational, surface_time, *, field_backend, surface_options):
        assert field_backend is backend
        assert surface_options is options
        if failure in ('reconstruct', 'both'):
            raise RuntimeError('injected reconstruct failure')
        assert coefficients is hydro and case is setup
        assert_allclose(computational.position_m, expected.position_m-[0., 0., offset])
        # Motion/radiation inputs are passed through unchanged.
        assert computational.velocity_native is motion.velocity_native
        vertices = expected.surface_vertices_m.copy()
        vertices[:, 2] -= offset
        return (vertices, expected.surface_faces,
                expected.surface_vertex_offsets, expected.surface_face_offsets)
    monkeypatch.setattr(producer, 'reconstruct_total_surface', reconstruct)
    if failure is not None:
        expected_error = 'reconstruct' if failure in ('reconstruct', 'both') else 'close'
        with pytest.raises(RuntimeError, match=f'injected {expected_error} failure') as caught:
            simulate(boat(), Sea([.1], [1.], [0.], [0.]), 1., influence=influence, options=options)
        if failure == 'both':
            assert any('close failure' in note for note in caught.value.__notes__)
        assert closed == [True]
        assert not (tmp_path/'scratch').exists()
        return
    actual = simulate(boat(), Sea([.1], [1.], [0.], [0.]), 1., influence=influence, options=options)
    assert_allclose(actual.position_m, expected.position_m)
    assert_allclose(actual.orientation_wxyz, expected.orientation_wxyz)
    assert_allclose(actual.surface_vertices_m, expected.surface_vertices_m)
    np.testing.assert_array_equal(actual.surface_faces, expected.surface_faces)
    assert not (tmp_path/'scratch').exists()


@pytest.mark.parametrize('end', [.01, .08, .09, 130., 130.01, 130.04])
def test_surface_clock_selects_native_samples_and_preserves_endpoint(end):
    from topwater.solve.sampling import surface_sample_times
    from topwater.solve.fourier import offline_radiation_response
    time = np.arange(round(end/.01)+1)*.01
    output = surface_sample_times(time)
    assert output[0] == 0. and output[-1] == time[-1]
    assert np.isin(output, time).all()
    velocity = np.zeros((len(time), 6))
    omega=2*np.pi*np.fft.rfftfreq(len(time),.01)
    response, _ = offline_radiation_response(
        omega, np.zeros((len(omega), 1, 6), complex), time, velocity, velocity,
        output_time_s=output)
    assert response.shape == (len(output), 1)
