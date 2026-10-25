"""Saved geometry rendering, pose conventions and fail-closed publication."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from topwater.result_render import (
    _boundary_enters_view, _projected_frame, _topology, render_result, sampled_poses,
)
from topwater.workflow_data import Result


def saved_result():
    hull = np.array([[-1., -1., -1.], [1., -1., -1.], [1., 1., -1.], [-1., 1., -1.],
                     [-1., -1., 1.], [1., -1., 1.], [1., 1., 1.], [-1., 1., 1.]])
    hull_faces = np.array([[0, 2, 1], [0, 3, 2], [0, 1, 5], [0, 5, 4],
                           [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6],
                           [3, 0, 4], [3, 4, 7]])
    water = np.array([[-5., -5., 0.], [5., -5., 0.], [5., 5., 0.], [-5., 5., 0.],
                      [-1., -1., 0.], [1., -1., 0.], [1., 1., 0.], [-1., 1., 0.]])
    faces = np.array([[0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
                      [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])
    water_next = np.vstack((water, [0., -5., 0.]))
    faces_next = np.vstack(([[0, 8, 5], [8, 1, 5]], faces[1:]))
    return Result(time_s=np.array([0., .04, .08]), position_m=np.zeros((3, 3)),
                  orientation_wxyz=np.tile([1., 0., 0., 0.], (3, 1)),
                  hull_vertices_m=hull, hull_faces=hull_faces, cg_m=np.zeros(3),
                  surface_time_s=np.array([0., .08]),
                  surface_vertices_m=np.vstack((water, water_next)),
                  surface_faces=np.vstack((faces, faces_next)),
                  surface_vertex_offsets=np.array([0, 8, 17]),
                  surface_face_offsets=np.array([0, 8, 17]))


def test_pose_interpolation_and_heading_relative_hud_preserve_initial_trim():
    original = saved_result()
    heading = Rotation.from_euler('z', .7)
    trim = Rotation.from_euler('y', .05)
    rotations = heading*Rotation.from_euler('z', np.array([0., .1, .2])[:, None])*trim
    result = replace(original, position_m=np.array([[0., 0., 0.], [1., 2., .1], [2., 4., .2]]),
                     orientation_wxyz=rotations.as_quat(scalar_first=True),
                     surface_time_s=np.array([0., .06, .08]),
                     surface_vertices_m=np.vstack((original.surface_vertices_m[:8],
                                                   original.surface_vertices_m[:8],
                                                   original.surface_vertices_m[8:])),
                     surface_faces=np.vstack((original.surface_faces[:8],
                                               original.surface_faces[:8], original.surface_faces[8:])),
                     surface_vertex_offsets=np.array([0, 8, 16, 25]),
                     surface_face_offsets=np.array([0, 8, 16, 25]))
    positions, matrices, angles, ambiguous = sampled_poses(result)
    np.testing.assert_allclose(positions[1], [1.5, 3., .15], atol=1e-14)
    expected = (heading*Rotation.from_euler('z', np.array([0., .15, .2])[:, None])*trim).as_matrix()
    np.testing.assert_allclose(matrices, expected, atol=1e-14)
    np.testing.assert_allclose(angles[:, :2], 0., atol=1e-12)
    np.testing.assert_allclose(angles[:, 2], np.rad2deg([0., .15, .2]), atol=1e-12)
    assert not ambiguous


def test_render_consumes_variable_topology_and_preserves_hull_hole(tmp_path):
    result = saved_result()
    destination = tmp_path/'video.mp4'
    destination.write_bytes(b'previous movie')
    summary = render_result(result, destination)
    assert summary['complete'] and summary['frames'] == 2
    assert summary['topology_changes'] == 1
    assert summary['size'] == (800, 800)
    assert destination.stat().st_size > 1000
    assert list(tmp_path.iterdir()) == [destination]


def test_projection_rejects_uncovered_interior_hole_outside_hull():
    from matplotlib import tri
    result = saved_result()
    water, faces = result.surface_frame(0)
    bad = water.copy()
    bad[4:, 0] += 2.
    loops = _topology(faces, len(bad))
    with pytest.raises(ValueError, match='uncovered hole'):
        _projected_frame(bad, faces, result.hull_vertices_m,
                         np.array([-1.45, -1.45]), np.array([1.45, 1.45]), loops, tri)


def test_projection_rejects_folded_layer_and_inconsistent_winding():
    from matplotlib import tri
    result = saved_result()
    water, faces = result.surface_frame(0)
    bad = faces.copy()
    bad[0] = bad[0, ::-1]
    with pytest.raises(ValueError, match='winding'):
        _topology(bad, len(water))
    folded = water.copy()
    folded[1, :2] = folded[0, :2]
    loops = _topology(faces, len(water))
    with pytest.raises(ValueError, match='coincident|folded'):
        _projected_frame(folded, faces, result.hull_vertices_m,
                         np.array([-1.45, -1.45]), np.array([1.45, 1.45]), loops, tri)


def test_outer_boundary_indent_is_not_hidden_by_corner_coverage():
    xy = np.array([[-5., -5.], [5., -5.], [5., 5.], [0., 0.], [-5., 5.]])
    assert _boundary_enters_view(xy, np.arange(len(xy)), np.array([-1., -1.]),
                                np.array([1., 1.]), 1e-10)


def test_failed_projection_preserves_previous_video(tmp_path, monkeypatch):
    import imageio_ffmpeg
    import topwater.result_render as renderer

    def fake_writer(path, *args, **kwargs):
        from pathlib import Path
        Path(path).write_bytes(b'incomplete new movie')
        while True:
            yield

    def failed_projection(*args, **kwargs):
        raise ValueError('projection rejected')

    monkeypatch.setattr(imageio_ffmpeg, 'write_frames', fake_writer)
    monkeypatch.setattr(renderer, '_projected_frame', failed_projection)
    destination = tmp_path/'video.mp4'
    destination.write_bytes(b'previous movie')
    with pytest.raises(ValueError, match='projection rejected'):
        render_result(saved_result(), destination)
    assert destination.read_bytes() == b'previous movie'
    assert list(tmp_path.iterdir()) == [destination]


def test_failed_decode_preserves_previous_video(tmp_path, monkeypatch):
    import imageio_ffmpeg
    import topwater.result_render as renderer

    def fake_writer(path, *args, **kwargs):
        from pathlib import Path
        Path(path).write_bytes(b'incomplete new movie')
        while True:
            yield

    def failed_decode(*args, **kwargs):
        raise RuntimeError('decode rejected')

    monkeypatch.setattr(imageio_ffmpeg, 'write_frames', fake_writer)
    monkeypatch.setattr(renderer, '_verify_video', failed_decode)
    destination = tmp_path/'video.mp4'
    destination.write_bytes(b'previous movie')
    with pytest.raises(RuntimeError, match='decode rejected'):
        render_result(saved_result(), destination)
    assert destination.read_bytes() == b'previous movie'


def test_nonuniform_surface_cadence_is_not_silently_resampled(tmp_path):
    original = saved_result()
    changed = replace(original, surface_time_s=np.array([0., .03, .08]),
                      surface_vertices_m=np.vstack((original.surface_vertices_m[:8],
                                                    original.surface_vertices_m[:8],
                                                    original.surface_vertices_m[8:])),
                      surface_faces=np.vstack((original.surface_faces[:8],
                                                original.surface_faces[:8], original.surface_faces[8:])),
                      surface_vertex_offsets=np.array([0, 8, 16, 25]),
                      surface_face_offsets=np.array([0, 8, 16, 25]))
    with pytest.raises(ValueError, match='uniformly spaced'):
        render_result(changed, tmp_path/'video.mp4')
