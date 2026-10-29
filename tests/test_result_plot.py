"""Saved six-motion conventions and failure-safe PNG publication."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from topwater.solve.plot import _rms_color_maximum
from topwater.solve.plot import _sea_spectrum
from topwater.solve.plot import motion_histories
from topwater.solve.plot import save_motion_plot
from topwater.workflow_data import Result
from topwater.workflow_data import Sea


def saved_motion_result(times=None):
    times = np.linspace(0., 60., 3001) if times is None else np.asarray(times, dtype=float)
    hull = np.array([[-1., -1., -1.], [1., -1., -1.], [1., 1., -1.], [-1., 1., -1.],
                     [-1., -1., 1.], [1., -1., 1.], [1., 1., 1.], [-1., 1., 1.]])
    faces = np.array([[0, 2, 1], [0, 3, 2], [0, 1, 5], [0, 5, 4],
                      [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6],
                      [3, 0, 4], [3, 4, 7]])
    water = np.array([[-5., -5., 0.], [5., -5., 0.], [0., 5., 0.]])
    return Result(time_s=times, position_m=np.zeros((len(times), 3)),
                  orientation_wxyz=np.tile([1., 0., 0., 0.], (len(times), 1)),
                  hull_vertices_m=hull, hull_faces=faces, cg_m=np.zeros(3),
                  surface_time_s=times[[0, -1]], surface_vertices_m=np.tile(water, (2, 1)),
                  surface_faces=np.array([[0, 1, 2], [0, 1, 2]]),
                  surface_vertex_offsets=np.array([0, 3, 6]),
                  surface_face_offsets=np.array([0, 1, 2]))


def test_original_sea_spectrum_bins_power_in_world_from_directions_without_phase():
    sea = Sea(np.array([.3, .4, .12, .8]),
        2*np.pi*np.array([.112, .113, .262, .612]),
        np.array([-np.pi, -np.pi+np.deg2rad(5), -np.pi/2, 0.]), np.zeros(4))
    angle, frequency, bins, height, heading = _sea_spectrum(sea)
    np.testing.assert_allclose(np.diff(angle), np.deg2rad(15), rtol=0., atol=1e-15)
    np.testing.assert_allclose(np.diff(frequency), .025, rtol=0., atol=1e-16)
    assert frequency[-1] == pytest.approx(.625)  # Include carriers above .5 Hz.
    expected = np.zeros((25, 24))
    expected[4, 0], expected[10, 6], expected[24, 12] = np.sqrt(.5)*np.array([.5, .12, .8])
    np.testing.assert_allclose(bins, expected, rtol=0., atol=1e-15)
    assert np.sum(bins**2) == pytest.approx(.5*np.sum(sea.amplitudes_m**2))
    assert height == pytest.approx(4*np.sqrt(.5*np.sum(sea.amplitudes_m**2)))
    assert height == pytest.approx(4*np.sqrt(np.sum(bins**2)))
    shifted = replace(sea, phases_rad=np.array([.7, -.4, 1., 2.]))
    repeated = _sea_spectrum(shifted)
    np.testing.assert_array_equal(repeated[2], bins)
    assert repeated[3:] == (height, heading)


def test_calm_and_ambiguous_default_sea_headings():
    calm = Sea(np.array([]), np.array([]), np.array([]), np.array([]))
    _, frequency, bins, height, heading = _sea_spectrum(calm)
    assert frequency[-1] == .5 and bins.shape == (20, 24)
    assert np.count_nonzero(bins) == 0 and height == heading == 0.
    opposed = Sea(np.ones(2), np.ones(2), np.array([0., np.pi]), np.zeros(2))
    with pytest.raises(ValueError, match='no unique default boat heading'):
        _sea_spectrum(opposed)
    with pytest.raises(TypeError, match='original Sea'):
        _sea_spectrum(object())


@pytest.mark.parametrize('peak, expected', [
    (0., 1.), (.011, .011), (.115420427, .12), (.1, .1), (.415, .42), (.51, .51), (2.51, 2.6)])
def test_colorbar_maximum_is_rounded_up_without_altering_bin_values(peak, expected):
    values = np.array([[peak/3, peak]])
    before = values.copy()
    maximum = _rms_color_maximum(values)
    assert maximum == pytest.approx(expected)
    assert maximum >= peak
    np.testing.assert_array_equal(values, before)


@pytest.mark.parametrize('propagation, heading_offset, peak, expected_ticks', [
    (.37, 0., .415, [0., .1, .2, .3, .4, .42]),
    (-1.2, .23, .115420427, [0., .02, .04, .06, .08, .1, .12])])
def test_combined_png_preserves_motion_axes_and_marks_actual_heading_relative_to_default(
        tmp_path, monkeypatch, propagation, heading_offset, peak, expected_ticks):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.image import imread
    actual = saved_motion_result([0., .04, .08])
    default = propagation+np.pi
    orientation = Rotation.from_euler('z', default+heading_offset).as_quat(scalar_first=True)
    actual = replace(actual, orientation_wxyz=np.tile(orientation, (3, 1)))
    sea = Sea(np.array([peak/np.sqrt(.5)]), np.array([2*np.pi*.112]),
              np.array([propagation]), np.zeros(1))
    observed = []
    print_png = FigureCanvasAgg.print_png

    def inspect(canvas, *args, **kwargs):
        figure = canvas.figure
        assert len(figure.axes) == 8  # Six histories, polar spectrum and colorbar.
        motion = figure.axes[:6]
        polar, bar = figure.axes[6:]
        assert polar.get_label() == 'sea spectrum'
        assert polar.get_theta_offset() == pytest.approx(np.pi/2-default)
        assert polar.get_theta_direction() == 1
        downward = np.deg2rad(polar.get_rlabel_position())+polar.get_theta_offset()
        assert np.mod(downward, 2*np.pi) == pytest.approx(3*np.pi/2)
        mesh = polar.collections[0]
        assert mesh.get_cmap().name == 'cividis'
        values = np.ma.asarray(mesh.get_array())
        zero = np.asarray(values.data) == 0.
        assert np.any(zero) and not np.any(np.ma.getmaskarray(values))
        zero_rgba = mesh.to_rgba(values)[zero]
        np.testing.assert_allclose(zero_rgba,
            np.tile(mesh.get_cmap()(0.), (int(np.count_nonzero(zero)), 1)), rtol=0., atol=1e-15)
        assert not np.allclose(mesh.get_cmap()(0.)[:3], [1., 1., 1.])
        assert polar.get_facecolor()[:3] == (1., 1., 1.)
        assert all(label.get_color() == 'black' for label in polar.get_yticklabels())
        np.testing.assert_allclose(polar.collections[0].get_coordinates()[0, :, 0],
                                   np.linspace(0., 2*np.pi, 25), atol=1e-15)
        assert all(label.get_text() == '' for label in polar.get_xticklabels())
        assert len(polar.get_xticks()) == 0
        assert not any(line.get_visible() for line in polar.get_xgridlines())
        assert all(line.get_visible() for line in polar.get_ygridlines())
        assert [label.get_text() for label in polar.get_yticklabels()] == ['0.1 Hz', '0.2', '0.3', '0.4']
        line = next(line for line in polar.lines if line.get_label() == 'initial boat heading')
        # The default is computed from Sea, not taken from this actual bow.
        screen = np.mod(line.get_xdata()+polar.get_theta_offset(), 2*np.pi)
        np.testing.assert_allclose(screen, np.pi/2+heading_offset, atol=1e-14)
        assert line.get_linestyle() == ':'
        assert line.get_color() == 'black'
        assert bar.get_xlabel() == 'Bin RMS elevation (m)'
        assert not figure.texts  # No figure supertitle/footer or polar title/subtitles.
        norm = polar.collections[0].norm
        assert norm.vmin == 0.
        assert norm.vmax == pytest.approx(expected_ticks[-1])
        assert np.max(polar.collections[0].get_array()) == pytest.approx(peak)
        np.testing.assert_allclose(bar.get_xlim(), [norm.vmin, norm.vmax], atol=1e-15)
        ticks = bar.get_xticks()
        np.testing.assert_allclose(ticks, expected_ticks, atol=1e-15)
        assert np.all((ticks >= 0.) & (ticks <= expected_ticks[-1]))
        assert min(polar.get_position().x0, bar.get_position().x0) > max(
            axis.get_position().x1 for axis in motion)
        assert bar.get_position().y1 < polar.get_position().y0
        assert bar.get_position().width > bar.get_position().height
        right_center = 14.75/17.5
        assert (polar.get_position().x0+polar.get_position().x1)/2 == pytest.approx(right_center)
        assert (bar.get_position().x0+bar.get_position().x1)/2 == pytest.approx(right_center)
        assert polar.get_position().height > .75
        assert polar.get_position().y1 > .95
        assert polar.get_position().width*17.5 == pytest.approx(5.1)
        assert polar.get_position().height*6.5 == pytest.approx(5.1)
        assert bar.get_position().width == pytest.approx(polar.get_position().width)
        for row in (motion[:3], motion[3:]):
            first = row[0].get_yticks()
            assert len(first) == 7 and first[3] == 0.
            np.testing.assert_allclose(first, -first[::-1], atol=1e-15)
            for axis in row:
                np.testing.assert_array_equal(axis.get_yticks(), first)
                np.testing.assert_allclose(axis.get_ylim(), [first[0], first[-1]], atol=1e-15)
                np.testing.assert_array_equal(axis.lines[0].get_xdata(), actual.time_s)
        assert all(axis.get_xlabel() == 'Time (s)' for axis in motion[3:])
        observed.append(True)
        return print_png(canvas, *args, **kwargs)

    monkeypatch.setattr(FigureCanvasAgg, 'print_png', inspect)
    path = tmp_path/'result.png'
    save_motion_plot(actual, path, sea=sea)
    assert observed == [True]
    assert imread(path).shape[:2] == (1170, 3150)


def test_all_native_samples_and_six_motions_with_initial_heading_and_trim():
    original = saved_motion_result()
    times = original.time_s
    expected = np.sin(times[:, None]*np.array([.1, .2, .3, .4, .5, .6]))
    expected *= np.array([.4, .7, .5, 8., 3., 2.])
    heading = Rotation.from_euler('z', .7)
    mean = Rotation.from_euler('y', .12)
    rotations = heading * Rotation.from_euler('xyz', expected[:, 3:], degrees=True) * mean
    positions = np.array([14., -8., 2.]) + expected[:, :3] @ heading.as_matrix().T
    result = replace(original, position_m=positions,
                     orientation_wxyz=rotations.as_quat(scalar_first=True))
    actual_time, actual = motion_histories(result)
    assert actual.shape == (3001, 6)
    np.testing.assert_array_equal(actual_time, times)
    np.testing.assert_allclose(actual, expected, atol=1e-12)
    np.testing.assert_allclose(actual[0], 0., atol=1e-12)
    np.testing.assert_array_equal(result.position_m, positions)


def test_yaw_is_unwrapped_without_dropping_startup():
    original = saved_motion_result([0., 1., 2., 3.])
    yaw = np.array([0., 170., 190., 350.])
    result = replace(original, orientation_wxyz=Rotation.from_euler(
        'z', yaw[:, None], degrees=True).as_quat(scalar_first=True))
    times, histories = motion_histories(result)
    np.testing.assert_array_equal(times, [0., 1., 2., 3.])
    np.testing.assert_allclose(histories[:, 5], yaw, atol=1e-12)


def test_gimbal_lock_warns_and_near_vertical_initial_bow_fails_clearly():
    original = saved_motion_result([0., 1.])
    result = replace(original, orientation_wxyz=Rotation.from_euler(
        'y', np.array([[0.], [90.]]), degrees=True).as_quat(scalar_first=True))
    with pytest.warns(UserWarning, match='ambiguous at gimbal lock'):
        _, values = motion_histories(result)
    assert np.isfinite(values).all()
    vertical = Rotation.from_euler('y', np.pi/2-1e-14).as_quat(scalar_first=True)
    with pytest.raises(ValueError, match='no horizontal heading'):
        motion_histories(replace(original, orientation_wxyz=np.tile(vertical, (2, 1))))


def test_png_dimensions_and_complete_history_publication(tmp_path, monkeypatch):
    from matplotlib.axes import Axes
    from matplotlib.image import imread
    recorded = []
    plot = Axes.plot

    def capture(axis, x, y, *args, **kwargs):
        recorded.append((np.array(x), np.array(y)))
        return plot(axis, x, y, *args, **kwargs)

    monkeypatch.setattr(Axes, 'plot', capture)
    destination = tmp_path/'case'/'result.png'
    destination.parent.mkdir()
    destination.write_bytes(b'old chart')
    result = saved_motion_result()
    assert save_motion_plot(result, destination) == destination
    assert imread(destination).shape[:2] == (1170, 2160)
    assert len(recorded) == 6
    for times, values in recorded:
        np.testing.assert_array_equal(times, result.time_s)
        assert values.shape == (3001,)
    assert list(destination.parent.iterdir()) == [destination]
    assert list(tmp_path.iterdir()) == [destination.parent]


@pytest.mark.parametrize('failure', ['encode', 'publish'])
def test_plot_failure_preserves_old_png_and_removes_staging(tmp_path, monkeypatch, failure):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    import topwater.solve.plot as module
    destination = tmp_path/'case'/'result.png'
    destination.parent.mkdir()
    destination.write_bytes(b'old chart')

    def fail(*args, **kwargs):
        raise OSError('deliberate plot failure')

    if failure == 'encode':
        monkeypatch.setattr(FigureCanvasAgg, 'print_png', fail)
    else:
        monkeypatch.setattr(module.os, 'replace', fail)
    with pytest.raises(OSError, match='deliberate plot failure'):
        save_motion_plot(saved_motion_result([0., 1.]), destination)
    assert destination.read_bytes() == b'old chart'
    assert list(destination.parent.iterdir()) == [destination]
    assert list(tmp_path.iterdir()) == [destination.parent]


def test_plot_rejects_wrong_contract_and_extension(tmp_path):
    with pytest.raises(TypeError, match='requires a Result'):
        motion_histories(object())
    with pytest.raises(ValueError, match='.png'):
        save_motion_plot(saved_motion_result([0., 1.]), tmp_path/'result.pdf')


@pytest.mark.parametrize('extent', [0., 1e-9, 1., 1e9])
def test_seven_nice_ticks_zero_middle_shared_rows_and_limits(tmp_path, monkeypatch, extent):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    import topwater.solve.plot as module
    times = np.array([0., 1., 2.])
    histories = extent*np.array([[0., 0., 0., 0., 0., 0.],
                                 [1., -.4, .7, -8., 3., .5],
                                 [-.2, .8, -.6, 4., -2., 1.]])
    monkeypatch.setattr(module, 'motion_histories', lambda result: (times, histories))
    original = FigureCanvasAgg.print_png
    observed = []

    def inspect(canvas, *args, **kwargs):
        observed.extend((axis.get_yticks().copy(), axis.get_ylim()) for axis in canvas.figure.axes)
        return original(canvas, *args, **kwargs)

    monkeypatch.setattr(FigureCanvasAgg, 'print_png', inspect)
    save_motion_plot(saved_motion_result(times), tmp_path/'case'/'result.png')
    assert len(observed) == 6
    for index, (ticks, limits) in enumerate(observed):
        assert len(ticks) == 7
        assert ticks[3] == 0.
        np.testing.assert_allclose(ticks, -ticks[::-1], rtol=0., atol=0.)
        np.testing.assert_array_equal(limits, ticks[[0, -1]])
        np.testing.assert_array_equal(ticks, observed[(index//3)*3][0])
        np.testing.assert_array_equal(limits, observed[(index//3)*3][1])
        assert np.max(abs(histories[:, index])) <= limits[1]/1.05
        step = ticks[4]
        nice = step/(10.**np.floor(np.log10(step)))
        assert any(np.isclose(nice, value) for value in (1., 2., 2.5, 5., 10.))
