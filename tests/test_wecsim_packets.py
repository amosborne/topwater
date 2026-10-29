"""Packet-to-native phase, direction, coherence and static-reference contracts."""
from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose

from topwater.solve.preparation import IncomingTrain
from topwater.solve.wecsim.spectral import _coherent_waves
from topwater.inputs import WaveComponents
from topwater.solve.packets import minimum_cutoff_packet


def test_phase_is_direct_at_dwell_midpoint_with_existing_negative_time_sign():
    phase, heading, amplitude = np.pi/2, .73, .2
    waves = WaveComponents([amplitude], [2*np.pi*.1], [heading], [phase])
    train = IncomingTrain(waves)
    packet = minimum_cutoff_packet(.1, -phase, 10., 5., 100., amplitude_m=amplitude)
    assert_allclose(train.elevation([[0., 0.]], packet.target.time_s)[:, 0],
                    packet.elevation_m, rtol=0, atol=1e-13)
    # phi=+pi/2 in cos(-omega*(t-mid)+phi) has a POSITIVE midpoint slope.
    assert train.elevation([[0., 0.]], [65.01])[0, 0] > train.elevation([[0., 0.]], [64.99])[0, 0]
    omega, _, _, directions, _, coefficients = _coherent_waves(train)
    assert_allclose(omega, train.omega[1:])
    # At oblique spatial targets, every actual sideband uses its own dispersion.
    distance = (2*np.pi/(2*np.pi*.1)**2*9.81)/4
    xy = np.outer([-distance, 0., distance], [np.cos(heading), np.sin(heading)])
    times = np.array([12.3, 65., 103.7])
    expected = np.full((len(times), len(xy)), train.dc_mean_m)
    for theta, column in zip(directions, coefficients.T):
        s = xy @ [np.cos(theta), np.sin(theta)]
        expected += (np.exp(-1j*times[:, None]*omega) @
                     (column[:, None]*np.exp(1j*omega[:, None]**2/9.81*s))).real
    assert_allclose(train.elevation(xy, times), expected, rtol=0, atol=1e-13)


def test_aggregate_budget_preserves_original_carriers_and_real_dc():
    waves = WaveComponents([.1, .2], [2*np.pi*.1, 2*np.pi*.2], [.2, 1.3], [.4, -.7])
    originals = [x.copy() for x in (waves.amplitudes, waves.angular_frequencies, waves.directions, waves.phases)]
    train = IncomingTrain(waves)
    assert_allclose([row['allocated_tolerance_m'] for row in train.packet_diagnostics], [.01/3, .02/3])
    assert train.certified_origin_exceedance_m <= .01
    assert_allclose(train.dc_mean_m, sum(train.dc_by_carrier_m), rtol=0, atol=0)
    assert np.all(train.C[:, :, 0] == 0)
    assert np.all(train.C[:, 1] == 0)
    directions, records = train.directional_elevation_records(np.array([0., 65., 130.]))
    assert_allclose(records.sum(axis=1), train.elevation([[0., 0.]], [0., 65., 130.])[:, 0], atol=1e-14)
    for actual, original in zip((waves.amplitudes, waves.angular_frequencies, waves.directions, waves.phases), originals):
        np.testing.assert_array_equal(actual, original)


def test_native_coherence_preserves_cancellation_and_does_not_add_zero_headings():
    C = np.zeros((3, 2, 3), complex)
    C[0, 0, 1:] = [1.+2j, -.3j]
    C[1, 0, 1:] = -C[0, 0, 1:]
    C[2, 0, 1:] = [.1, .2j]
    train = SimpleNamespace(omega=np.array([0., .1, .2]), theta=np.array([.2, .2, .7]),
                            C=C, carriers=SimpleNamespace(amplitudes=np.ones(3)))
    *_, directions, degrees, coefficients = _coherent_waves(train)
    assert_allclose(directions, [.2, .7])
    assert_allclose(coefficients[:, 0], 0., atol=0.)
    assert_allclose(coefficients[:, 1], C[2, 0, 1:])


def test_empty_sea_has_valid_exactly_zero_native_spectrum():
    train = SimpleNamespace(omega=np.array([0., .1, .2]), theta=np.array([]),
        C=np.zeros((0, 2, 3), complex), carriers=SimpleNamespace(amplitudes=np.array([])))
    *_, directions, degrees, coefficients = _coherent_waves(train)
    assert len(directions) == 2
    assert_allclose(coefficients, np.zeros((2, 2)))


def test_nyquist_singleton_is_not_doubled_when_packet_reaches_grid_endpoint(monkeypatch):
    import topwater.solve.preparation as setup_module
    # A controlled packet isolates adapter normalization; no bandwidth search.
    packet = SimpleNamespace(positive_components=4, coefficients=np.array([.02, .1j, .2, .3j, .04]),
        target=SimpleNamespace(time_s=np.arange(8)*.01), dc_m=.02, cutoff_hz=50.,
        diagnostics={'certified_envelope_exceedance_bound_m': 0.})
    monkeypatch.setattr(setup_module, 'minimum_cutoff_packet', lambda *args, **kwargs: packet)
    waves = WaveComponents([.2], [2*np.pi*.1], [.3], [.4])
    train = IncomingTrain(waves, rest_s=0., ramp_s=.02, dwell_s=.04)
    assert_allclose(train.C[0, 0, -1], .04)
    assert_allclose(train.C[0, 0, 1:-1], np.conj(2*packet.coefficients[1:-1]))


def test_static_mean_coordinate_shift_retains_full_restoring_couplings():
    # All six loads must see the hydrostatic heave column, including moments.
    stiffness = np.arange(36, dtype=float).reshape(6, 6)
    stiffness = stiffness+stiffness.T+np.eye(6)*100
    physical_q = np.array([.1, -.2, .03, .02, -.04, .01])
    for mean in [-.02, 0., .013]:
        native_q = physical_q-np.array([0., 0., mean, 0., 0., 0.])
        assert_allclose(-stiffness@native_q,
                        -stiffness@physical_q+mean*stiffness[:, 2], atol=1e-14)


@pytest.mark.parametrize('solver', ['native', 'wecsim'])
def test_solve_surface_wrapper_keeps_geometry_computational_and_exports_world_once(monkeypatch, tmp_path, solver):
    from contextlib import contextmanager
    from topwater.solve import workflow
    from topwater.solve.native import backend as native_backend
    from topwater.solve.wecsim import backend as wecsim_backend
    from topwater.solve.types import Motion, PreparedMotion
    motion_backend = native_backend if solver == 'native' else wecsim_backend
    d = .007
    positions = np.array([[0., 0., .2], [.1, 0., .21]])
    stamp = np.array([0., .1])
    motion = Motion(stamp, np.zeros((2, 6)), np.zeros((2, 6)), np.zeros((2, 6)),
        positions.copy(), np.tile([1., 0., 0., 0.], (2, 1)), np.zeros((2, 6)))
    setup = SimpleNamespace(world_vertical_offset_m=d,
        train=SimpleNamespace(certify_incident_startup=lambda setup: None))
    water = np.array([[0., 0., -.001], [1., 0., 0.], [0., 1., .001]])
    faces = np.array([[0, 1, 2]])
    hydro = object()
    selected = []
    closed = []
    field_backend = SimpleNamespace(close=lambda: closed.append(True))
    def field(coefficients, device, **kwargs):
        assert coefficients is hydro
        selected.append(device)
        return field_backend
    def reconstruct(coefficients, actual_setup, computational, times, **kwargs):
        assert coefficients is hydro and actual_setup is setup
        assert_allclose(computational.position_m, positions-[0., 0., d])
        assert_allclose(motion.position_m, positions)
        assert (tmp_path/'surface-world-offset.json').is_file()
        assert kwargs['field_backend'] is field_backend
        return water.copy(), faces, np.array([0, 3]), np.array([0, 1])
    @contextmanager
    def scratch(name):
        assert name == solver
        yield tmp_path
    monkeypatch.setattr(motion_backend, 'preflight', lambda options, influence: None)
    monkeypatch.setattr(motion_backend, 'prepare', lambda *args: PreparedMotion(hydro))
    monkeypatch.setattr(motion_backend, 'integrate', lambda *args: motion)
    monkeypatch.setattr(workflow, 'prepare_case', lambda *args, **kwargs: setup)
    monkeypatch.setattr(workflow, 'waves_from_sea', lambda sea: object())
    monkeypatch.setattr(workflow, 'scratch_directory', scratch)
    monkeypatch.setattr(workflow, 'select_field_backend', field)
    monkeypatch.setattr(workflow.producer, 'reconstruct_total_surface', reconstruct)
    monkeypatch.setattr(workflow, 'Result', lambda **kwargs: SimpleNamespace(**kwargs))
    boat = SimpleNamespace(to_body=lambda: object(), vertices_m=water, faces=faces, cg_m=np.zeros(3))
    result = workflow.solve(boat, object(), 100., solver=solver)
    assert_allclose(result.position_m, positions)
    assert_allclose(result.surface_vertices_m, water+[0., 0., d])
    assert_allclose(water[:, 2], [-.001, 0., .001])
    assert selected == ['auto' if solver == 'native' else 'cpu']
    assert closed == [True]
