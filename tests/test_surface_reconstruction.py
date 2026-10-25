"""Analytical field-convention fixtures; no Capytaine or MATLAB jobs.

Execution still belongs to the coordinator's serial numerical allocation.
"""
import numpy as np
import pytest
from types import SimpleNamespace

from topwater._capytaine_surface import (
    SurfaceKernel, convolve_radiation, recover_transfer, synthesize_spectrum,
    require_surface_readiness, cosine_kernel, exterior_field_points,
    exterior_trace_points, incident_spectrum,
    WaterGraphPartition, _regular_grid, _grid_p1_values, _delta_moments,
    _total_waterline, _waterline_distance_xy,
)


def exponential_kernel(*, accepted=True):
    t = np.arange(20001)*.001
    values = np.zeros((len(t), 1, 6))
    values[:, 0, 0] = 3*np.exp(-2*t)
    direct = np.zeros((1, 6)); direct[0, 0] = .4
    return SurfaceKernel(t, values, direct, accepted)


def test_exponential_plus_direct_has_correct_complex_phase():
    kernel = exponential_kernel()
    omega = np.array([.07, .31, 1.13, 3.71, 7.19])
    actual = recover_transfer(kernel, omega)[:, 0, 0]
    expected = .4 + 3/(2-1j*omega)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-6)
    assert np.all(actual.imag > 0), 'exp(-iwt) uses the positive-imaginary causal transform'


def test_cosine_inversion_has_bounded_finite_band_error_and_refines():
    # h=A exp(-a*t), H=D+A/(a-i*w).  The infinite-frequency tail is bounded
    # uniformly in lag by (2*A/pi)*(pi/2-atan(W/a)).  The trapezoid bound follows
    # from |(R cos(w*t))''| <= 10*A/a^3+4*A*t/a^2+A*t^2/a,
    # where R=A*a/(a^2+w^2).  Thus the test budgets actual truncation and
    # quadrature error instead of comparing a finite band to an exact identity.
    A, a, D = 3., 2., .4
    lag = np.linspace(0., .25, 21)
    errors = []
    for maximum, spacing in [(64., .02), (256., .01)]:
        omega = np.arange(round(maximum/spacing)+1)*spacing
        transfer = np.zeros((len(omega), 1, 6), dtype=complex)
        transfer[:, 0, 0] = D + A/(a-1j*omega)
        direct = np.zeros((1, 6)); direct[0, 0] = D
        candidate = cosine_kernel(omega, transfer, lag, direct_velocity=direct)
        error = float(np.max(abs(candidate.values[:, 0, 0]-A*np.exp(-a*lag))))
        tail = 2*A/np.pi*(np.pi/2-np.arctan(maximum/a))
        derivative_bound = 10*A/a**3 + 4*A*lag[-1]/a**2 + A*lag[-1]**2/a
        quadrature = 2/np.pi*maximum*spacing**2/12*derivative_bound
        assert error <= tail+quadrature+1e-12
        np.testing.assert_array_equal(candidate.direct_velocity, direct)
        assert not candidate.accepted
        errors.append(error)
    assert errors[1] < errors[0]/3


def test_mean_footprint_guard_rejects_interior_and_on_waterline():
    vertices = np.array([[-2., -1., 0.], [2., -1., 0.], [2., 1., 0.],
                         [-2., 1., 0.], [0., 0., -1.]])
    hydro = SimpleNamespace(body=SimpleNamespace(mesh=SimpleNamespace(vertices=vertices)))
    exterior = np.array([[2.1, 0.], [-3., .5], [0., 1.1]])
    np.testing.assert_array_equal(exterior_field_points(hydro, exterior), exterior)
    for point in ([[0., 0.]], [[2., 0.]], [[-2., -1.]]):
        with pytest.raises(ValueError, match='strictly outside'):
            exterior_field_points(hydro, point)


def test_exterior_trace_maps_only_disturbance_queries():
    vertices = np.array([[-2., -1., 0.], [2., -1., 0.], [2., 1., 0.],
                         [-2., 1., 0.], [0., 0., -1.]])
    hydro = SimpleNamespace(body=SimpleNamespace(mesh=SimpleNamespace(vertices=vertices)))
    original = np.array([[1.9, 0.], [3., 0.]])
    mapped, continued = exterior_trace_points(hydro, original, .025)
    np.testing.assert_allclose(mapped, [[2.025, 0.], [3., 0.]], atol=1e-14)
    np.testing.assert_array_equal(continued, [True, False])
    train = SimpleNamespace(omega=np.array([0., 1.]), k=np.array([0., .2]),
                            theta=np.array([0.]), C=np.array([[[0., 1.], [0., 0.]]]))
    np.testing.assert_allclose(incident_spectrum(train, original)[1], np.exp(1j*.2*original[:, 0]))
    assert abs(incident_spectrum(train, mapped)[1, 0]-incident_spectrum(train, original)[1, 0]) > .01


def test_identified_acceleration_feedthrough_is_explicit_and_uses_actual_acceleration():
    kernel = exponential_kernel()
    kernel.direct_acceleration = np.zeros((1, 6)); kernel.direct_acceleration[0, 0] = .07
    omega = np.array([.17, 1.13, 3.71])
    np.testing.assert_allclose(recover_transfer(kernel, omega)[:, 0, 0],
                               .4+3/(2-1j*omega)-1j*omega*.07, atol=2e-6)
    t = kernel.lag_s[:101]
    velocity = np.zeros((len(t), 6))
    acceleration = np.zeros_like(velocity); acceleration[:, 0] = 2.
    with pytest.raises(ValueError, match='Actual modal accelerations'):
        convolve_radiation(kernel, t, velocity)
    np.testing.assert_allclose(convolve_radiation(kernel, t, velocity, acceleration)[:, 0], .14, atol=1e-14)


def test_causal_convolution_has_no_future_motion_or_startup_history():
    kernel = exponential_kernel()
    t = kernel.lag_s[:4001]
    w = 1.7
    v = np.zeros((len(t), 6)); v[:, 0] = np.sin(w*t)
    actual = convolve_radiation(kernel, t, v)[:, 0]
    expected = .4*np.sin(w*t) + 3*(2*np.sin(w*t)-w*np.cos(w*t)+w*np.exp(-2*t))/(4+w*w)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-6)
    changed = v.copy(); changed[2001:, 0] += 100
    np.testing.assert_allclose(convolve_radiation(kernel, t, changed)[:2001, 0], actual[:2001], rtol=0, atol=1e-11)
    assert actual[0] == pytest.approx(0, abs=1e-14)


def test_unaccepted_or_short_kernel_cannot_publish_radiation():
    kernel = exponential_kernel(accepted=False)
    t = kernel.lag_s[:11]
    with pytest.raises(ValueError, match='review'):
        convolve_radiation(kernel, t, np.zeros((len(t), 6)))
    kernel.accepted = True
    kernel.lag_s = kernel.lag_s[:5]
    kernel.values = kernel.values[:5]
    with pytest.raises(ValueError, match='entire requested history'):
        convolve_radiation(kernel, t, np.zeros((len(t), 6)))


def test_weighted_spectrum_has_no_hidden_factor_or_conjugation():
    omega = np.array([0., .3, .9])
    coefficients = np.array([[0.], [2+3j], [-.7+.2j]])
    t = np.array([0., .17, 1.2])
    expected = sum((coefficients[i, 0]*np.exp(-1j*w*t)).real for i, w in enumerate(omega))
    np.testing.assert_allclose(synthesize_spectrum(omega, coefficients, t)[:, 0], expected, atol=1e-14)


def test_total_surface_capability_checks_selected_policy(monkeypatch):
    import topwater._capytaine_surface as producer
    require_surface_readiness()
    monkeypatch.setattr(producer, '_SURFACE_POLICY_ACCEPTED', False)
    with pytest.raises(ValueError, match='pilot'):
        require_surface_readiness()


def test_wet_set_area_is_geometric_and_owns_shared_vertical_prism_edges_once():
    xy, faces, _ = _regular_grid([[-2., -2.], [2., 2.]], 2.)
    graph = WaterGraphPartition(xy, faces)
    # A vertical physical rectangle lies exactly on the grid's x=0 boundary.
    hull = np.array([[0., -1., -1.], [0., 1., -1.], [0., 1., 1.], [0., -1., 1.]])
    physical = np.array([[0, 1, 2], [0, 2, 3]])
    fraction, wet = graph.wet_set_area(hull, physical, np.full(len(xy), .5), [0., 0., 1.], 0., 2.)
    assert fraction == pytest.approx(.5, abs=1e-9)
    assert wet == pytest.approx(3., abs=1e-9)
    # Current pose moved upward .25m; transport the original material wet plane.
    moved = hull+[0., 0., .25]
    fraction, wet = graph.wet_set_area(moved, physical, np.zeros(len(xy)), [0., 0., 1.], -.25, 2.)
    assert fraction == pytest.approx(.25, abs=1e-9)
    assert wet == pytest.approx(1.5, abs=1e-9)


def test_wet_set_area_uses_symmetric_difference_not_only_wetted_area_change():
    xy, faces, _ = _regular_grid([[-2., -2.], [2., 2.]], 2.)
    graph = WaterGraphPartition(xy, faces)
    hull = np.array([[0., -1., -1.], [0., 1., -1.], [0., 1., 1.], [0., -1., 1.]])
    physical = np.array([[0, 1, 2], [0, 2, 3]])
    fraction, wet = graph.wet_set_area(hull, physical, np.zeros(len(xy)), [0., 1., 0.], 0., 2.)
    assert wet == pytest.approx(2., abs=1e-9)
    assert fraction == pytest.approx(1., abs=1e-9)


def test_wet_set_area_fails_on_incomplete_water_graph_coverage():
    xy, faces, _ = _regular_grid([[-.5, -.5], [.5, .5]], .5)
    graph = WaterGraphPartition(xy, faces)
    hull = np.array([[0., -1., -1.], [0., 1., -1.], [0., 1., 1.], [0., -1., 1.]])
    with pytest.raises(ValueError, match='partition'):
        graph.wet_set_area(hull, [[0, 1, 2], [0, 2, 3]], np.zeros(len(xy)), [0., 0., 1.], 0., 2.)


def test_p1_cut_interpolation_and_projected_delta_moments():
    xy, faces, axes = _regular_grid([[0., 0.], [1., 1.]], 1.)
    value = xy[:, 0]+2*xy[:, 1]
    query = np.array([[.2, .1], [.1, .9]])
    np.testing.assert_allclose(_grid_p1_values(axes, value, query), query[:, 0]+2*query[:, 1])
    # Unit patch with a quarter-square physical hole, preserving an actual
    # total-waterline inner boundary rather than using the hull projection.
    xy = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.],
                   [.25, .25], [.75, .25], [.75, .75], [.25, .75]])
    faces = np.array([[0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
                      [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])
    hull = np.array([[x, y, z] for z in (1., 6.) for x, y in
                     ((.25, .25), (.75, .25), (.75, .75), (.25, .75))])
    value = xy[:, 0]+2*xy[:, 1]
    vertices = np.c_[xy, value+2]
    moments = _delta_moments(vertices, faces, value, hull, 10., [[0., 0.], [1., 1.]])
    for result in moments.values():
        np.testing.assert_allclose(result[:5], [.75, 1.5, 3., 2., 133/64], atol=1e-14)
        assert result[7] >= 2


def test_near_band_uses_actual_waterline_not_flared_hull_outline():
    xy = np.array([[-3., -3.], [3., -3.], [3., 3.], [-3., 3.],
                   [-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])
    faces = np.array([[0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
                      [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])
    vertices = np.c_[xy, np.zeros(len(xy))]
    hull = np.array([[s*x, s*y, z] for s, z in ((.5, -1.), (1.5, 1.))
                     for x, y in ((-1., -1.), (1., -1.), (1., 1.), (-1., 1.))])
    line = _total_waterline(vertices, faces, hull, [[-3., -3.], [3., 3.]])
    assert len(line) == 4
    # This point is inside the flared hull's horizontal outline, but is .4m
    # from its actual total-waterline contact and outside a .2m waterline band.
    distance = _waterline_distance_xy(np.array([[1.4, 0.]]), line)
    np.testing.assert_allclose(distance, [.4], atol=1e-14)
    assert distance[0] > .2
    with pytest.raises(ValueError, match='no resolved hull-contact waterline'):
        _total_waterline(vertices[:4], np.array([[0, 1, 2], [0, 2, 3]]), hull,
                         [[-3., -3.], [3., 3.]])


@pytest.fixture
def calm_surface_case(monkeypatch, tmp_path):
    """Geometry/integration control with analytical zero fields, no BEM job."""
    from scipy.spatial.transform import Rotation
    from topwater.workflow_data import Boat
    import topwater._capytaine_surface as producer

    vertices = np.array([[-1., -1., 0.], [1., -1., 0.], [1., 1., 0.], [-1., 1., 0.],
                         [-1., -1., 1.], [1., -1., 1.], [1., 1., 1.], [-1., 1., 1.]])
    physical = np.array([[0, 2, 1], [0, 3, 2], [0, 1, 5], [0, 5, 4],
                         [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])
    boat = Boat(vertices, physical, 1000., np.array([0., 0., .5]), np.diag([400., 400., 600.]))
    body = boat.to_body()
    embedding = Rotation.from_euler('z', .7).as_matrix()
    stamp = np.array([0., .02])
    orientation = np.tile(Rotation.from_matrix(embedding).as_quat(scalar_first=True), (2, 1))
    train = SimpleNamespace(omega=np.array([0., 1.]), k=np.array([0., 1/9.81]),
        theta=np.empty(0), C=np.zeros((0, 2, 2), dtype=complex),
        carriers=SimpleNamespace(amplitudes=np.empty(0)),
        excitation=lambda *args: (np.zeros((2, 6)), np.zeros(6)))
    setup = SimpleNamespace(body=body, physical_faces=body.mesh.faces[body.mesh.physical_faces],
        world_from_native=embedding, equilibrium=SimpleNamespace(orientation=np.array([1., 0., 0., 0.])),
        mean_cg_m=np.zeros(3), hydrostatic=SimpleNamespace(wetted_area=8.), train=train)
    motion = SimpleNamespace(time_s=stamp, position_m=np.zeros((2, 3)), orientation_wxyz=orientation,
                             velocity_native=np.zeros((2, 6)))
    waterline = np.array([[-1., -1., 0.], [1., -1., 0.], [1., 1., 0.], [-1., 1., 0.], [0., 0., -.5]])
    hydro = SimpleNamespace(directory=tmp_path/'hydrodynamics',
        body=SimpleNamespace(mesh=SimpleNamespace(vertices=waterline)), depth=np.inf,
        excitation_omega=np.array([.1, 1.]), excitation_headings=np.array([0., 2*np.pi]),
        excitation=np.zeros((2, 2, 6), dtype=complex), K=np.zeros((6, 6)))
    monkeypatch.setattr(producer, 'radiation_velocity_transfer',
        lambda h, xy, omega, **kw: np.zeros((len(omega), len(xy), 6), dtype=complex))
    monkeypatch.setattr(producer, 'combined_diffraction',
        lambda h, tr, xy, **kw: (np.zeros((2, len(xy)), dtype=complex),
                                 np.zeros((2, len(xy)), dtype=complex), np.zeros((2, 6), dtype=complex)))
    return SimpleNamespace(producer=producer, hydro=hydro, setup=setup, motion=motion, stamp=stamp,
        vertices=vertices, physical=physical, boat=boat, embedding=embedding, orientation=orientation)


def test_production_calm_geometry_packs_contact_and_covers_oblique_world_view(calm_surface_case, capsys):
    from topwater.workflow_data import Result
    case = calm_surface_case
    producer, hydro, setup, motion, stamp = case.producer, case.hydro, case.setup, case.motion, case.stamp
    vertices, physical, boat = case.vertices, case.physical, case.boat
    embedding, orientation = case.embedding, case.orientation
    water, faces, vertex_offsets, face_offsets = producer.reconstruct_total_surface(hydro, setup, motion, stamp)
    result = Result(stamp, motion.position_m, orientation, vertices, physical, boat.cg_m, stamp,
                    water, faces, vertex_offsets, face_offsets)
    for i in range(2):
        frame, connectivity = result.surface_frame(i)
        np.testing.assert_allclose(frame[:, 2], 0., atol=1e-14)
        posed = (vertices-boat.cg_m) @ embedding.T
        center = (posed[:, :2].min(axis=0)+posed[:, :2].max(axis=0))/2
        half = np.max(posed[:, :2].max(axis=0)-posed[:, :2].min(axis=0))/2*1.45
        assert np.all(frame[:, :2].min(axis=0) <= center-half)
        assert np.all(frame[:, :2].max(axis=0) >= center+half)
        assert connectivity.max() < len(frame)
    diagnostic = capsys.readouterr().err
    assert 'degenerate incident RMS' in diagnostic
    assert 'no total disturbed-fluid velocity bound' in diagnostic


def test_failed_frame_preserves_primary_context_and_replays_without_physics(calm_surface_case, monkeypatch, tmp_path, capsys):
    import json
    from topwater.workflow_data import Result
    case = calm_surface_case
    producer = case.producer
    original_metrics = producer._delta_moments
    original_error = ValueError('sentinel contact diagnostic failure')

    def failing_metrics(*args, **kwargs):
        raise original_error

    monkeypatch.setattr(producer, '_delta_moments', failing_metrics)
    with pytest.raises(producer.SurfaceExportError, match='frame 0, time 0s') as caught:
        producer.reconstruct_total_surface(case.hydro, case.setup, case.motion, case.stamp)
    error = caught.value
    assert error.__cause__ is original_error
    assert error.stage == 'waterline-delta-metrics'
    assert error.frame_index == 0 and error.time_s == 0
    cache = error.recovery_directory/'completed-field.npz'
    assert cache.is_file()
    report = json.loads((error.recovery_directory/'failure.json').read_text())
    assert report['stage'] == error.stage
    assert 'sentinel contact diagnostic failure' in report['traceback']
    assert 'sentinel contact diagnostic failure' in capsys.readouterr().err
    with np.load(error.recovery_directory/'failure-frame.npz', allow_pickle=False) as saved:
        assert saved['water_vertices_world_m'].shape[1] == 3
        assert saved['hull_faces'].shape[1] == 3
        assert saved['reference_normal_native'].shape == (3,)
    monkeypatch.setattr(producer, '_delta_moments', original_metrics)

    def no_physics(*args, **kwargs):
        raise AssertionError('A geometry replay must not execute a field producer')

    monkeypatch.setattr(producer, 'radiation_velocity_transfer', no_physics)
    monkeypatch.setattr(producer, 'combined_diffraction', no_physics)
    monkeypatch.setattr(producer, 'incident_spectrum', no_physics)
    output = tmp_path/'replayed-result.npz'
    producer.replay_surface_cache(cache, diagnostic_directory=tmp_path/'replay-diagnostics', output_path=output)
    result = Result.load(output)
    np.testing.assert_array_equal(result.position_m, case.motion.position_m)
    np.testing.assert_array_equal(result.orientation_wxyz, case.motion.orientation_wxyz)
    np.testing.assert_array_equal(result.surface_time_s, case.stamp)
    np.testing.assert_array_equal(np.unique(result.hull_faces), np.arange(len(result.hull_vertices_m)))
    np.testing.assert_array_equal(result.cg_m, case.setup.body.cg)
    with np.load(cache, allow_pickle=False) as saved:
        # Artificial planar closure introduces a center used only by its caps.
        # Keep it in the recovery geometry, but not the public physical mesh.
        assert len(saved['body_vertices_m']) > len(result.hull_vertices_m)
        np.testing.assert_array_equal(saved['body_faces'], case.setup.body.mesh.faces)
        np.testing.assert_array_equal(saved['opening_rim'], case.setup.body.mesh.opening_rim)
        np.testing.assert_array_equal(saved['artificial_faces'], case.setup.body.mesh.artificial_faces)
        np.testing.assert_allclose(result.hull_vertices_m[result.hull_faces],
                                   saved['body_vertices_m'][saved['physical_faces']], atol=0, rtol=0)


@pytest.mark.parametrize('bad_time', ['missing', 'empty', 'nonnumeric'])
def test_malformed_completed_time_keeps_original_validation_cause(bad_time, tmp_path):
    import json
    from topwater._capytaine_surface import SurfaceExportError, replay_surface_cache
    data = dict(cache_version=np.array(1), producer_source_sha256=np.array('diagnostic-fixture'))
    if bad_time == 'empty':
        data['surface_time_s'] = np.empty(0)
    elif bad_time == 'nonnumeric':
        data['surface_time_s'] = np.array(['invalid-time'])
    cache = tmp_path/'malformed-field.npz'
    np.savez(cache, **data)
    directory = tmp_path/'malformed-diagnostics'
    with pytest.raises(SurfaceExportError, match='completed-field-validation-or-summary') as caught:
        replay_surface_cache(cache, diagnostic_directory=directory)
    cause = caught.value.__cause__
    assert isinstance(cause, KeyError if bad_time == 'missing' else ValueError)
    assert caught.value.time_s is None
    report = json.loads((directory/'failure.json').read_text())
    assert report['time_s'] is None
    assert report['error_type'] == type(cause).__name__
    assert report['error_message'] == str(cause)
