"""Saved six-motion conventions and failure-safe PNG publication."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from topwater.result_plot import motion_histories, save_motion_plot
from topwater.workflow_data import Result


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
    import topwater.result_plot as module
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
