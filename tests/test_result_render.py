"""Saved geometry rendering, pose conventions and fail-closed publication."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from topwater.result_render import (
    _boundary_enters_view, _contact_contours, _mask_path, _projected_frame,
    _submerged_edges_from_contact, _topology, _wetted_display,
    render_result, sampled_poses,
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


def test_render_consumes_variable_topology_and_preserves_hull_hole(tmp_path, capsys):
    result = saved_result()
    destination = tmp_path/'video.mp4'
    destination.write_bytes(b'previous movie')
    summary = render_result(result, destination)
    assert summary['complete'] and summary['frames'] == 2
    assert summary['topology_changes'] == 1
    assert summary['size'] == (800, 800)
    for name in ('validation_seconds', 'mask_seconds', 'edge_seconds', 'contact_seconds'):
        assert summary[name] >= 0.
    progress = capsys.readouterr().out
    assert 'render progress: 1/2 frames' in progress
    assert 'render progress: 2/2 frames' in progress
    assert destination.stat().st_size > 1000
    assert list(tmp_path.iterdir()) == [destination]


def test_full_dry_wireframe_black_outline_and_solid_blue_waterline_layers(tmp_path, monkeypatch):
    from pathlib import Path
    from matplotlib.axes import Axes
    from matplotlib.collections import LineCollection
    from matplotlib.patches import PathPatch
    from topwater._visualization_geometry import mesh_edges
    from topwater.surface_geometry import waterline_segments
    import imageio_ffmpeg
    import topwater.result_render as module
    original = saved_result()
    hull, hull_faces, _, _ = flared_display_case()
    result = replace(original, hull_vertices_m=hull)
    drawn = []
    draw_artist = Axes.draw_artist

    def capture(axes, artist):
        snapshot = dict(label=artist.get_label(), zorder=artist.get_zorder())
        if isinstance(artist, LineCollection):
            snapshot.update(segments=np.asarray(artist.get_segments()).copy(),
                colors=artist.get_colors().copy(), linestyle=artist.get_linestyle(),
                effects=artist.get_path_effects())
        if isinstance(artist, PathPatch):
            snapshot['mask_vertices'] = artist.get_path().vertices.copy()
        drawn.append(snapshot)
        return draw_artist(axes, artist)

    def fake_writer(path, *args, **kwargs):
        Path(path).write_bytes(b'private layer-test stub')
        while True:
            yield

    def unexpected_submerged_clipping(*args, **kwargs):
        pytest.fail('full wireframe rendering must not clip dry edges away')

    monkeypatch.setattr(Axes, 'draw_artist', capture)
    monkeypatch.setattr(imageio_ffmpeg, 'write_frames', fake_writer)
    monkeypatch.setattr(module, '_verify_video', lambda *args: None)
    monkeypatch.setattr(module, '_submerged_edges_from_contact', unexpected_submerged_clipping)
    render_result(result, tmp_path/'layers.mp4')
    assert len(drawn) == 10  # Water, water-domain mask, then all three line layers per frame.
    first = drawn[:5]
    assert [item['zorder'] for item in first] == [0, 1, 2, 3, 4]
    assert [item['label'] for item in first[2:]] == ['boat mesh', 'boat outline', 'physical waterline']
    edges, rim = mesh_edges(hull_faces)
    np.testing.assert_array_equal(first[2]['segments'], hull[edges, :2])
    np.testing.assert_array_equal(first[3]['segments'], hull[rim, :2])
    # Dry rim is at half-width2, while the mask and waterline remain at1.25.
    assert np.max(abs(first[3]['segments'])) == pytest.approx(2.)
    assert np.max(abs(first[1]['mask_vertices'])) == pytest.approx(1.25)
    water, faces = result.surface_frame(0)
    contact = waterline_segments(hull, hull_faces, water, faces)
    np.testing.assert_allclose(first[4]['segments'], contact[:, :, :2], atol=1e-12)
    np.testing.assert_allclose(first[2]['colors'], [[176/255, 176/255, 176/255, 1.]])
    np.testing.assert_array_equal(first[3]['colors'], [[0., 0., 0., 1.]])
    np.testing.assert_array_equal(first[4]['colors'], [[0., 0., 1., 1.]])
    assert all(dashes is None for offset, dashes in first[4]['linestyle'])
    assert not first[4]['effects']  # Solid blue, without a white halo.


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


def tiny_corner_patch():
    # One 10 m triangular patch, with its corner cut into a 10 micrometre face.
    # The tiny face is a well-conditioned right triangle, not a thin sliver.
    points = np.array([[0., 0., 0.], [10., 0., 0.], [0., 10., 0.],
                       [1e-5, 0., 0.], [0., 1e-5, 0.]])
    faces = np.array([[0, 3, 4], [3, 1, 2], [3, 2, 4]])
    return points, faces


def test_projection_retains_tiny_well_conditioned_cut_triangle():
    from matplotlib import tri
    points, faces = tiny_corner_patch()
    loops = _topology(faces, len(points))
    hull = saved_result().hull_vertices_m*.25 + [2., 2., 0.]
    projected = _projected_frame(points, faces, hull,
                                 np.array([1.5, 1.5]), np.array([2.5, 2.5]), loops, tri)
    np.testing.assert_array_equal(projected.triangles, faces)
    assert projected.get_trifinder()(2., 2.) >= 0


@pytest.mark.parametrize('change', ['collinear', 'negative', 'nearvertical', 'folded'])
def test_small_face_projection_still_rejects_unsupported_geometry(change):
    from matplotlib import tri
    points, faces = tiny_corner_patch()
    loops = _topology(faces, len(points))
    if change == 'collinear':
        points[4] = [2e-5, 0., 1.]
    elif change == 'negative':
        faces[0] = faces[0, ::-1]
    elif change == 'nearvertical':
        # A nondegenerate 3D face with numerically collinear projected edges.
        points[4] = [1e-5, 1e-22, 1.]
    else:
        points[4] = [1e-5, 0., 1.]
    hull = saved_result().hull_vertices_m*.25 + [2., 2., 0.]
    with pytest.raises(ValueError, match='nondegenerate upward|coincident|folded'):
        _projected_frame(points, faces, hull,
                         np.array([1.5, 1.5]), np.array([2.5, 2.5]), loops, tri)


def test_projection_rejects_nonfinite_computation_from_finite_points():
    from matplotlib import tri
    points, faces = tiny_corner_patch()
    # All supplied coordinates remain finite; large-face area arithmetic cannot.
    points *= 1e154
    loops = _topology(faces, len(points))
    hull = (saved_result().hull_vertices_m*.25 + [2., 2., 0.])*1e154
    with pytest.raises(ValueError, match='nonfinite area computations'):
        _projected_frame(points, faces, hull,
                         np.array([1.5, 1.5])*1e154, np.array([2.5, 2.5])*1e154, loops, tri)


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


def test_failed_decode_preserves_previous_video(tmp_path, monkeypatch, capsys):
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
    failure = capsys.readouterr().err
    assert 'phase=full video verification' in failure
    assert 'completed=2/2 frames' in failure
    assert 'mask=' in failure and 'encoder_close=' in failure and 'active_phase_elapsed=' in failure


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


def flared_display_case():
    result = saved_result()
    hull = result.hull_vertices_m.copy()
    hull[:4, :2] *= .5
    hull[4:, :2] *= 2.
    water = np.array([[-5., -5., 0.], [5., -5., 0.], [5., 5., 0.], [-5., 5., 0.]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    return hull, result.hull_faces, water, faces


def projected_area(polygons):
    return sum(abs(np.sum(p[:, 0]*np.roll(p[:, 1], -1)
                          -p[:, 1]*np.roll(p[:, 0], -1)))/2 for p in polygons)


def test_wetted_mask_reaches_waterline_without_hiding_dry_flared_hull():
    hull, hull_faces, water, faces = flared_display_case()
    originals = [value.copy() for value in (hull, hull_faces, water, faces)]
    masks, wet_edges, wet_rim, contact = _wetted_display(
        hull, hull_faces, water, faces, _topology(faces, len(water)))
    # At z=0 the half-width is1.25; the dry gunwale half-width is2.
    assert projected_area(masks) == pytest.approx(2.5**2)
    assert max(abs(np.concatenate(masks)).ravel()) == pytest.approx(1.25)
    assert len(wet_edges) and np.max(wet_edges[:, :, 2]) <= 1e-12
    assert not len(wet_rim)
    assert len(contact)
    for current, original in zip((hull, hull_faces, water, faces), originals):
        np.testing.assert_array_equal(current, original)


@pytest.mark.parametrize('posed', [False, True])
def test_sloping_saved_surface_and_moving_hull_share_exact_contact_boundary(posed):
    from scipy.spatial import ConvexHull
    hull, hull_faces, water, faces = flared_display_case()
    if posed:
        hull = hull @ Rotation.from_euler('xyz', [.2, .1, .5]).as_matrix().T + [1., -.5, .1]
    water[:, 2] = .2*water[:, 0]+.1*water[:, 1]
    masks, wet_edges, _, contact = _wetted_display(
        hull, hull_faces, water, faces, _topology(faces, len(water)))
    assert len(masks) and len(contact) and len(wet_edges)
    assert np.all(wet_edges[:, :, 2] <= .2*wet_edges[:, :, 0]+.1*wet_edges[:, :, 1]+1e-12)
    planes = ConvexHull(hull).equations
    boundary = []
    for polygon in masks:
        xyz = np.column_stack((polygon, .2*polygon[:, 0]+.1*polygon[:, 1]))
        np.testing.assert_array_less(xyz @ planes[:, :3].T+planes[:, 3], 1e-10)
        for start, end in zip(xyz, np.roll(xyz, -1, axis=0)):
            if np.any((abs(planes[:, :3] @ start+planes[:, 3]) < 1e-10)
                      & (abs(planes[:, :3] @ end+planes[:, 3]) < 1e-10)):
                if np.linalg.norm(end[:2]-start[:2]) > 1e-12:
                    boundary.append([start[:2], end[:2]])
    boundary = np.asarray(boundary)
    for point in contact.reshape(-1, 3):
        direction = boundary[:, 1]-boundary[:, 0]
        fraction = np.clip(np.sum((point[:2]-boundary[:, 0])*direction, axis=1)
                           /np.sum(direction**2, axis=1), 0., 1.)
        distance = np.linalg.norm(boundary[:, 0]+fraction[:, None]*direction-point[:2], axis=1)
        assert distance.min() < 1e-9


def test_overtopped_opening_completes_display_volume_without_changing_physical_contact(monkeypatch):
    from scipy.spatial import ConvexHull
    from topwater.boat_geometry import physical_hull
    from topwater.surface_geometry import waterline_segments
    import topwater.result_render as module
    hull, hull_faces, water, faces = flared_display_case()
    hull = hull @ Rotation.from_euler('xyz', [.2, .1, .5]).as_matrix().T + [1., -.5, .1]
    water[:, 2] = .2*water[:, 0]+.1*water[:, 1]
    rim_gap = hull[4:, 2]-.2*hull[4:, 0]-.1*hull[4:, 1]
    assert rim_gap.min() < 0. < rim_gap.max()  # Actual partial opening overtopping.
    physical = waterline_segments(hull, hull_faces, water, faces)
    with pytest.raises(ValueError, match='contact is open'):
        _contact_contours(physical, np.linalg.norm(np.ptp(hull, axis=0)))
    called_faces = []

    def record_contact(vertices, triangles, *args, **kwargs):
        called_faces.append(np.array(triangles))
        return waterline_segments(vertices, triangles, *args, **kwargs)

    monkeypatch.setattr(module, 'waterline_segments', record_contact)
    masks, _, wet_rim, returned_contact = _wetted_display(
        hull, hull_faces, water, faces, _topology(faces, len(water)))
    np.testing.assert_array_equal(returned_contact, physical)
    assert len(wet_rim)
    # Independent planar section of the closed convex body: intersect its edges
    # with the known water plane. No invented elevations or triangle masks.
    closed = physical_hull(hull, hull_faces)
    assert len(called_faces) == 2
    np.testing.assert_array_equal(called_faces[0], hull_faces)
    np.testing.assert_array_equal(called_faces[1], closed.faces[closed.artificial_faces])
    edge_ids = np.unique(np.sort(np.concatenate((closed.faces[:, [0, 1]],
        closed.faces[:, [1, 2]], closed.faces[:, [2, 0]])), axis=1), axis=0)
    edges = closed.vertices[edge_ids]
    gap = edges[:, :, 2]-.2*edges[:, :, 0]-.1*edges[:, :, 1]
    crossing = gap[:, 0]*gap[:, 1] < 0.
    fraction = gap[crossing, 0]/(gap[crossing, 0]-gap[crossing, 1])
    section = edges[crossing, 0]+fraction[:, None]*(edges[crossing, 1]-edges[crossing, 0])
    assert projected_area(masks) == pytest.approx(ConvexHull(section[:, :2]).volume, abs=1e-10)


def test_closed_physical_contact_has_no_display_closure_or_extra_intersection(monkeypatch):
    import topwater.boat_geometry as geometry
    import topwater.result_render as module
    hull, hull_faces, water, faces = flared_display_case()
    kernel = module.waterline_segments
    calls = []

    def record_contact(*args, **kwargs):
        calls.append(True)
        return kernel(*args, **kwargs)

    def unexpected_closure(*args, **kwargs):
        pytest.fail('closed physical contact must not compute an opening closure')

    monkeypatch.setattr(module, 'waterline_segments', record_contact)
    monkeypatch.setattr(geometry, 'physical_hull', unexpected_closure)
    masks, _, _, _ = _wetted_display(hull, hull_faces, water, faces, _topology(faces, len(water)))
    assert len(masks) == 1
    assert len(calls) == 1


def test_dry_hull_is_transparent_and_fully_submerged_edges_remain_visible():
    hull, hull_faces, water, faces = flared_display_case()
    loops = _topology(faces, len(water))
    masks, edges, rim, contact = _wetted_display(hull+[0., 0., 3.], hull_faces, water, faces, loops)
    assert not masks and not len(edges) and not len(rim) and not len(contact)
    masks, edges, rim, contact = _wetted_display(hull-[0., 0., 3.], hull_faces, water, faces, loops)
    # Water is above the complete body: no air/water interface inside solid.
    assert not masks and len(edges) and len(rim) and not len(contact)
    assert np.max(edges[:, :, 2]) < 0.


def test_missing_water_hole_cannot_hide_behind_transparent_dry_gunwale():
    hull, hull_faces, _, _ = flared_display_case()
    original = saved_result()
    water, faces = original.surface_frame(0)
    water = water.copy()
    water[4:, :2] *= 1.5  # Inside old gunwale footprint, outside actual wetted volume.
    with pytest.raises(ValueError, match='hole is exposed'):
        _wetted_display(hull, hull_faces, water, faces, _topology(faces, len(water)))


def test_submerged_edge_preserves_multiple_crossings_between_dry_endpoints():
    from matplotlib.tri import Triangulation
    water = np.array([[x, y, z] for x, z in enumerate([-1., 1., -1., 1., -1.])
                      for y in [-1., 1.]])
    faces = np.array([[2*i, 2*i+2, 2*i+3] for i in range(4)]
                     + [[2*i, 2*i+3, 2*i+1] for i in range(4)])
    hull_edge = np.array([[0., 0., 0.], [4., 0., 0.]])
    contact = np.array([[[x, 0., 0.], [x, .5, 0.]] for x in [.5, 1.5, 2.5, 3.5]])
    segments, _ = _submerged_edges_from_contact(hull_edge, np.array([[0, 1]]), contact,
        Triangulation(water[:, 0], water[:, 1], faces), water, faces, 1e-9)
    assert segments.shape == (2, 2, 3)
    np.testing.assert_allclose(segments[:, :, 0], [[.5, 1.5], [2.5, 3.5]], atol=1e-12)
    np.testing.assert_array_equal(segments[:, :, 1:], 0.)


@pytest.mark.parametrize('start,end,expected', [
    ([-2., 0., -1.], [2., 0., -1.], [[-2., -1.], [1., 2.]]),
    ([-6., 0., -1.], [6., 0., -1.], [[-5., -1.], [1., 5.]]),
])
def test_submerged_edge_does_not_bridge_missing_water_between_covered_ends(start, end, expected):
    from matplotlib.tri import Triangulation
    original = saved_result()
    water, faces = original.surface_frame(0)
    loops = _topology(faces, len(water))
    edge = np.array([start, end])
    segments, owners = _submerged_edges_from_contact(edge, np.array([[0, 1]]),
        np.empty((0, 2, 3)), Triangulation(water[:, 0], water[:, 1], faces),
        water, faces, 1e-9, loops)
    assert segments.shape == (2, 2, 3)
    np.testing.assert_allclose(segments[:, :, 0], expected, atol=1e-12)
    np.testing.assert_array_equal(owners, [0, 0])
    # Both displayed portions have real saved coverage. No invented span runs
    # through the central x=-1..1 hole, even with covered physical endpoints.
    assert np.all((segments[:, 1, 0] <= -1.+1e-12) | (segments[:, 0, 0] >= 1.-1e-12))


def test_mask_has_one_contact_subpath_independent_of_water_triangle_count():
    hull, hull_faces, _, _ = flared_display_case()
    xx, yy = np.meshgrid(np.linspace(-5., 5., 21), np.linspace(-5., 5., 21))
    water = np.column_stack((xx.ravel(), yy.ravel(), np.zeros(xx.size)))
    faces = []
    for row in range(20):
        for column in range(20):
            a = row*21+column
            faces.extend([[a, a+1, a+22], [a, a+22, a+21]])
    faces = np.array(faces)
    masks, _, _, contact = _wetted_display(hull, hull_faces, water, faces, _topology(faces, len(water)))
    assert len(contact) > 4
    assert len(masks) == 1  # One actual loop, not one mask for each water tile.
    assert projected_area(masks) == pytest.approx(6.25)


def test_compound_mask_preserves_concavity_and_clockwise_disconnected_fill():
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.collections import PolyCollection
    from matplotlib.figure import Figure
    from matplotlib.patches import PathPatch
    polygons = [np.array([[0., 0.], [2., 0.], [2., 1.], [0., 1.]]),
                np.array([[0., 1.], [1., 1.], [1., 2.], [0., 2.]]),
                np.array([[3., 0.], [3., 1.], [4., 1.], [4., 0.]])]
    outputs = []
    for compound in (False, True):
        figure = Figure(figsize=(4, 3), dpi=40, facecolor='black')
        FigureCanvasAgg(figure)
        axes = figure.add_axes((0, 0, 1, 1))
        axes.set_xlim(-.5, 4.5); axes.set_ylim(-.5, 2.5); axes.set_axis_off()
        if compound:
            axes.add_patch(PathPatch(_mask_path(polygons), facecolor='white',
                                     edgecolor='none', antialiased=False))
        else:
            axes.add_collection(PolyCollection(polygons, facecolor='white',
                                               edgecolor='none', antialiased=False))
        figure.canvas.draw()
        pixels = np.asarray(figure.canvas.buffer_rgba())
        samples = []
        for point in ([.5, .5], [.5, 1.5], [1.5, 1.5], [3.5, .5], [2.5, .5]):
            x, y = axes.transData.transform(point).astype(int)
            samples.append(bool(pixels[pixels.shape[0]-1-y, x, 0]))
        outputs.append(samples)
        figure.clear()
    assert outputs == [[True, True, False, True, False]]*2
    assert len(_mask_path([]).vertices) == 0


def oriented_contact(polygon):
    xyz = np.column_stack((polygon, np.zeros(len(polygon))))
    return np.stack((xyz, np.roll(xyz, -1, axis=0)), axis=1)


def test_contact_contours_keep_concavity_disconnected_components_and_nested_hole():
    outer = np.array([[0., 0.], [0., 4.], [4., 4.], [4., 0.]])
    hole = np.array([[1., 1.], [3., 1.], [3., 3.], [1., 3.]])
    concave = np.array([[5., 0.], [5., 3.], [6., 3.], [6., 1.], [8., 1.], [8., 0.]])
    contact = np.concatenate([oriented_contact(polygon) for polygon in (outer, hole, concave)])
    original = contact.copy()
    contours = _contact_contours(contact[::-1], 8.)
    assert len(contours) == 3
    np.testing.assert_allclose(sorted(projected_area([polygon]) for polygon in contours), [4., 5., 16.])
    np.testing.assert_array_equal(contact, original)
    # Probe the actual Agg nonzero-winding fill; compound Path.contains_points
    # does not implement the renderer's nested winding semantics.
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.patches import PathPatch
    figure = Figure(figsize=(9, 5), dpi=30, facecolor='black')
    FigureCanvasAgg(figure)
    axes = figure.add_axes((0, 0, 1, 1))
    axes.set_xlim(-.5, 8.5); axes.set_ylim(-.5, 4.5); axes.set_axis_off()
    axes.add_patch(PathPatch(_mask_path(contours), facecolor='white', edgecolor='none'))
    figure.canvas.draw()
    pixels = np.asarray(figure.canvas.buffer_rgba())
    for point, expected in (([.5, .5], True), ([2., 2.], False),
                            ([5.5, 2.], True), ([7., 2.], False), ([7., .5], True)):
        x, y = axes.transData.transform(point).astype(int)
        assert bool(pixels[pixels.shape[0]-1-y, x, 0]) == expected
    figure.clear()


@pytest.mark.parametrize('invalid', ['open', 'branched', 'crossed', 'wrong_hole_winding'])
def test_unsupported_contact_contours_fail_clearly(invalid):
    outer = oriented_contact(np.array([[0., 0.], [0., 4.], [4., 4.], [4., 0.]]))
    if invalid == 'open':
        contact = outer[:-1]
    elif invalid == 'branched':
        contact = np.concatenate((outer, np.array([[[0., 0., 0.], [2., 2., 0.]]])))
    elif invalid == 'crossed':
        contact = oriented_contact(np.array([[0., 0.], [4., 4.], [4., 0.], [0., 4.]]))
    else:
        contact = np.concatenate((outer, oriented_contact(
            np.array([[1., 1.], [1., 3.], [3., 3.], [3., 1.]]))))
    with pytest.raises(ValueError, match='waterline|boundaries'):
        _contact_contours(contact, 8.)
