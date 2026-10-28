"""Producer-neutral scene contract; no graphics context or simulation needed."""
import numpy as np
import pytest

from topwater.cli import main
from topwater.visualization import SceneSequence, RenderSettings, convex_envelope, render_sequence
from topwater._visualization_geometry import waterline_segments, mesh_edges
from topwater._visualization_geometry import projected_footprint
from topwater._visualization_2d import header_text


def pyramid():
    points = np.array([[-2., -1., 0.], [2., -1., 0.], [2., 1., 0.], [-2., 1., 0.], [0., 0., 2.]])
    # Deliberately open base and mixed winding: convexity is geometric.
    faces = np.array([[0, 1, 4], [4, 2, 1], [2, 3, 4], [4, 0, 3]])
    return points, faces


def sequence():
    points, faces = pyramid()
    water = np.array([[-20., -20., 0.], [20., -20., 0.], [20., 20., 0.], [-20., 20., 0.]])
    return SceneSequence(np.array([0., .08]), np.array([points, points + [0., 0., .1]]), faces,
                         np.array([water, water]), np.array([[0, 1, 2], [0, 2, 3]]),
                         {'units': 'm', 'coordinate_system': 'right-handed-z-up', 'producer': 'synthetic'})


def test_prepared_scene_roundtrip(tmp_path):
    original = sequence()
    original.save(tmp_path / 'scene.npz')
    loaded = SceneSequence.load(tmp_path / 'scene.npz')
    for name in ('time_s', 'hull_points', 'hull_faces', 'water_points', 'water_faces'):
        np.testing.assert_array_equal(getattr(original, name), getattr(loaded, name))
    assert original.metadata == loaded.metadata


def test_open_convex_mesh_and_scale_independence():
    points, faces = pyramid()
    for scale in (1e-3, 1., 1e3):
        closure, oriented = convex_envelope(points * scale, faces)
        triangles = points[closure]
        normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        assert np.all(np.einsum('ij,ij->i', normals, triangles.mean(axis=1) - points.mean(axis=0)) > 0)
        assert len(oriented) == len(faces)


def test_rejects_nonconvex_boundary_instead_of_replacing_it():
    points, faces = pyramid()
    points = np.vstack((points, [0., 0., .5]))
    with pytest.raises(ValueError, match='convex boundary'):
        convex_envelope(points, np.vstack((faces, [0, 1, 5])))


def test_rejects_wrong_frame_units_and_time():
    data = sequence()
    fields = {name: getattr(data, name) for name in ('time_s', 'hull_points', 'hull_faces', 'water_points', 'water_faces')}
    with pytest.raises(ValueError, match='Normalize'):
        SceneSequence(**fields, metadata={'units': 'ft'})
    fields['time_s'] = np.array([1., 0.])
    with pytest.raises(ValueError, match='increase'):
        SceneSequence(**fields)


def test_rejects_float_connectivity_and_bad_indices():
    data = sequence()
    fields = {name: getattr(data, name) for name in ('time_s', 'hull_points', 'hull_faces', 'water_points', 'water_faces')}
    fields['hull_faces'] = data.hull_faces.astype(float)
    with pytest.raises(ValueError, match='integer'):
        SceneSequence(**fields)
    fields['hull_faces'] = np.array([[0, 1, 999]])
    with pytest.raises(ValueError, match='outside'):
        SceneSequence(**fields)


def test_irregular_cadence_cannot_silently_change_physical_time(tmp_path):
    scene = sequence()
    irregular = SceneSequence(np.array([0., .08, .17]), np.concatenate((scene.hull_points, scene.hull_points[:1])),
        scene.hull_faces, np.concatenate((scene.water_points, scene.water_points[:1])), scene.water_faces)
    with pytest.raises(ValueError, match='uniformly spaced'):
        render_sequence(irregular, tmp_path / 'wrong_speed.mp4')
    assert not (tmp_path / 'wrong_speed.mp4').exists()


def test_cli_help_does_not_need_a_graphics_context(capsys):
    with pytest.raises(SystemExit) as stopped:
        main(['render', '--help'])
    assert stopped.value.code == 0
    help_text = capsys.readouterr().out
    assert 'directory' in help_text and 'docs/workflow.md' in help_text
    for unsupported in ('--start', '--end', '--size', '--output'):
        assert unsupported not in help_text


@pytest.mark.parametrize('arguments', [[], ['render', 'directory', 'extra-output.mp4']])
def test_cli_rejects_missing_command_or_extra_render_path(arguments):
    with pytest.raises(SystemExit) as stopped:
        main(arguments)
    assert stopped.value.code == 2


def test_cli_passes_saved_result_and_fixed_output_file(tmp_path, monkeypatch):
    from topwater import result_render
    from tests.test_workflow_data import result
    saved = result()
    saved.save(tmp_path/'result.npz')
    captured = {}

    def capture(loaded, output, *, settings):
        captured.update(result=loaded, output=output, settings=settings)
        report = dict.fromkeys(('preparation_seconds', 'setup_seconds', 'geometry_seconds',
                                'draw_seconds', 'encoder_start_seconds', 'encoder_send_seconds',
                                'encoder_close_seconds', 'verification_seconds', 'elapsed_s'), 0.)
        report.update(complete=True, frames=len(loaded.surface_time_s), fps=12.5,
                      size=(800, 800), source_end_s=float(loaded.surface_time_s[-1]),
                      video_duration_s=len(loaded.surface_time_s)/12.5)
        return report

    monkeypatch.setattr(result_render, 'render_result', capture)
    assert main(['render', str(tmp_path)]) == 0
    assert captured['output'] == tmp_path/'video.mp4'
    np.testing.assert_array_equal(captured['result'].time_s, saved.time_s)
    assert isinstance(captured['settings'], RenderSettings)
    assert captured['settings'].maximum_seconds == max(240., float(len(saved.surface_time_s)))


def test_output_is_an_mp4_filename_not_a_directory(tmp_path):
    with pytest.raises(ValueError, match='complete MP4 file'):
        render_sequence(sequence(), tmp_path / 'directory')


def test_telemetry_roundtrip_and_validation(tmp_path):
    scene = sequence()
    scene.telemetry = {'roll_deg': np.array([1., 2.]), 'pitch_deg': np.array([.1, .2]), 'heave_m': np.array([.2, .3])}
    scene.save(tmp_path / 'telemetry.npz')
    loaded = SceneSequence.load(tmp_path / 'telemetry.npz')
    for key, values in scene.telemetry.items():
        np.testing.assert_array_equal(loaded.telemetry[key], values)
    fields = {key: getattr(scene, key) for key in ('time_s', 'hull_points', 'hull_faces', 'water_points', 'water_faces')}
    with pytest.raises(ValueError, match='one value per time'):
        SceneSequence(**fields, telemetry={'roll_deg': np.array([1.])})


def test_topdown_header_never_invents_motion_values(tmp_path):
    with pytest.raises(ValueError, match='requires producer telemetry'):
        render_sequence(sequence(), tmp_path / 'header.mp4')


def test_waterline_triangle_intersections_with_coarse_flat_water():
    from scipy.spatial import ConvexHull
    points = np.array([[x, y, z] for x in (-1., 1.) for y in (-1., 1.) for z in (-1., 1.)])
    faces = ConvexHull(points).simplices
    water = np.array([[-2., -2., 0.], [2., -2., 0.], [2., 2., 0.], [-2., 2., 0.]])
    water_faces = np.array([[0, 1, 2], [0, 2, 3]])
    segments = waterline_segments(points, faces, water, water_faces)
    np.testing.assert_allclose(segments[:, :, 2], 0., atol=1e-12)
    np.testing.assert_allclose(abs(segments[:, :, :2]).max(axis=2), 1., atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(segments[:, 1] - segments[:, 0], axis=1).sum(), 8., atol=1e-12)
    assert len(waterline_segments(points, faces, water + [0., 0., 3.], water_faces)) == 0


def test_mesh_rim_is_derived_from_open_boundary():
    _, faces = pyramid()
    all_edges, rim = mesh_edges(faces)
    assert len(all_edges) == 8 and len(rim) == 4
    assert set(rim.ravel()) == {0, 1, 2, 3}


def test_header_columns_stay_fixed_across_signs_and_digit_changes():
    first = header_text(9.92, -9.99, 0., -.1)
    second = header_text(50.00, 10.00, -1., 1.2)
    assert len(first) == len(second)
    for label in ('TIME', 'ROLL', 'PITCH', 'HEAVE'):
        assert first.index(label) == second.index(label)


def test_opaque_footprint_uses_projected_input_geometry():
    points, _ = pyramid()
    for scale in (1., 2.):
        outline = projected_footprint(points * scale + [100., -50., 3.])
        area = abs(np.dot(outline[:, 0], np.roll(outline[:, 1], -1)) -
                   np.dot(outline[:, 1], np.roll(outline[:, 0], -1))) / 2
        assert len(outline) == 4
        np.testing.assert_allclose(area, 8. * scale**2)
