"""Saved geometry rendering, pose conventions and fail-closed publication."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from topwater.render.geometry import _boundary_enters_view
from topwater.render.geometry import _contact_contours
from topwater.render.geometry import _projected_frame
from topwater.render.geometry import _submerged_edges_from_contact
from topwater.render.geometry import _topology
from topwater.render.geometry import _wetted_display
from topwater.render.workflow import render_result
from topwater.render.geometry import sampled_poses
from topwater.workflow_data import Result


@pytest.fixture(autouse=True)
def native_cpu_video(monkeypatch):
    from topwater.render import workflow
    from topwater.render.cpu.backend import DrawingBackend
    monkeypatch.setattr(workflow,'DrawingBackend',lambda size,lower,upper,**kwargs:DrawingBackend(size,lower,upper))


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


def test_public_scene_api_uses_shared_native_cpu_renderer_for_nonrigid_mixed_winding(tmp_path):
    from topwater.render.scene import SceneSequence
    from topwater.render.scene import render_sequence
    hull,faces,water,water_faces=flared_display_case()
    next_hull=hull.copy()
    next_hull[:,0]*=1.1  # Explicit mesh scenes need not be rigid pose records.
    next_hull[:,2]+=.05
    mixed=faces.copy();mixed[::2]=mixed[::2,::-1]
    scene=SceneSequence(np.array([0.,.08]),np.stack((hull,next_hull)),mixed,
        np.stack((water,water)),water_faces,telemetry=dict(roll_deg=np.zeros(2),
            pitch_deg=np.zeros(2),heave_m=np.array([0.,.05])))
    summary=render_sequence(scene,tmp_path/'scene.mp4')
    assert summary['frames']==2 and summary['backend']=='Skia'
    np.testing.assert_array_equal(scene.hull_faces,mixed)


def test_graticule_texture_is_cached_once_and_header_changes(tmp_path, monkeypatch):
    import topwater.render.hud as hud
    import topwater.render.cpu.backend as backend
    cached, headers, frames = [], [], []
    original_graticule, original_header, original_draw = hud.HUD.graticule, hud.HUD.header, backend.DrawingBackend.draw
    def graticule(self):
        rgba=original_graticule(self);cached.append(rgba.copy());return rgba
    def header(self,value):
        headers.append(value);return original_header(self,value)
    def draw(self,*args):
        rgba=original_draw(self,*args);frames.append(rgba.copy());return rgba
    monkeypatch.setattr(hud.HUD,'graticule',graticule)
    monkeypatch.setattr(hud.HUD,'header',header)
    monkeypatch.setattr(backend.DrawingBackend,'draw',draw)
    summary=render_result(saved_result(),tmp_path/'cached.mp4')
    assert len(cached)==1 and cached[0].shape==(800,800,4)
    assert len(headers)==len(frames)==2 and headers[0]!=headers[1]
    assert np.any(cached[0][:,:,3]==0)
    opaque=cached[0][:,:,3]==255
    for frame in frames:
        np.testing.assert_allclose(frame[opaque],cached[0][opaque],atol=1)
    assert summary['backend']=='Skia'
    assert summary['hud']['gauge_metres_per_pixel']==summary['hud']['common_metres_per_pixel']


@pytest.mark.parametrize('change,expected_checks', [('height', 1), ('xy', 2), ('topology', 2)])
def test_projection_cache_requires_identical_xy_and_connectivity(tmp_path, monkeypatch, change, expected_checks):
    from pathlib import Path
    import imageio_ffmpeg
    import topwater.render.geometry as module
    water = np.array([[-5., -5., 0.], [5., -5., 0.], [5., 5., 0.], [-5., 5., 0.]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    next_water, next_faces = water.copy(), faces.copy()
    next_water[:, 2] = .2
    if change == 'xy':
        next_water[0, 0] -= .1
    elif change == 'topology':
        next_faces = np.array([[0, 1, 3], [1, 2, 3]])
    result = replace(saved_result(), surface_vertices_m=np.vstack((water, next_water)),
        surface_faces=np.vstack((faces, next_faces)), surface_vertex_offsets=np.array([0, 4, 8]),
        surface_face_offsets=np.array([0, 2, 4]))
    checks = []
    original = module._projected_frame
    def check(*args):
        checks.append(True)
        return original(*args)
    def writer(path, *args, **kwargs):
        Path(path).write_bytes(b'projection-cache test')
        while True:
            yield
    monkeypatch.setattr(module, '_projected_frame', check)
    monkeypatch.setattr(imageio_ffmpeg, 'write_frames', writer)
    monkeypatch.setattr(module, '_verify_video', lambda *args: None)
    render_result(result, tmp_path/'cached.mp4')
    assert len(checks) == expected_checks


def test_fixed_projection_cache_never_skips_moving_hull_hole_checks(tmp_path):
    original = saved_result()
    water, faces = original.surface_frame(0)
    result = replace(original, position_m=np.array([[0., 0., 0.], [.05, 0., 0.], [.1, 0., 0.]]),
        surface_vertices_m=np.vstack((water, water)), surface_faces=np.vstack((faces, faces)),
        surface_vertex_offsets=np.array([0, 8, 16]), surface_face_offsets=np.array([0, 8, 16]))
    destination = tmp_path/'previous.mp4'
    destination.write_bytes(b'previous movie')
    with pytest.raises(ValueError, match='uncovered hole'):
        render_result(result, destination)
    assert destination.read_bytes() == b'previous movie'


def test_projection_cache_rejects_later_folded_xy(tmp_path):
    water = np.array([[-5., -5., 0.], [5., -5., 0.], [5., 5., 0.], [-5., 5., 0.]])
    next_water = water.copy()
    next_water[1] = [*water[0, :2], .2]  # Valid 3-D face, coincident XY projection.
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    result = replace(saved_result(), surface_vertices_m=np.vstack((water, next_water)),
        surface_faces=np.vstack((faces, faces)), surface_vertex_offsets=np.array([0, 4, 8]),
        surface_face_offsets=np.array([0, 2, 4]))
    destination = tmp_path/'previous.mp4'
    destination.write_bytes(b'previous movie')
    with pytest.raises(ValueError, match='coincident|folded'):
        render_result(result, destination)
    assert destination.read_bytes() == b'previous movie'


def test_origin_camera_failure_preserves_old_video_when_original_zoom_cannot_cover_trajectory(tmp_path):
    original = saved_result()
    translated = replace(original, position_m=np.tile([10., 0., 0.], (3, 1)),
                         surface_vertices_m=original.surface_vertices_m+[10., 0., 0.])
    output = tmp_path/'video.mp4'
    output.write_bytes(b'previous movie')
    with pytest.raises(ValueError, match='origin-centred camera'):
        render_result(translated, output)
    assert output.read_bytes() == b'previous movie'


def test_full_dry_wireframe_black_outline_and_solid_blue_waterline_layers(tmp_path, monkeypatch):
    from topwater.render.cpu.backend import DrawingBackend
    from topwater.render.mesh import mesh_edges
    from topwater.geometry.waterline import waterline_segments
    hull,hull_faces,_,_=flared_display_case()
    result=replace(saved_result(),hull_vertices_m=hull)
    drawn=[];original=DrawingBackend.draw
    def draw(self,water,faces,masks,layers):
        drawn.append((masks,[(name,np.array(segments),color,width,dashed) for name,segments,color,width,dashed in layers]))
        return original(self,water,faces,masks,layers)
    monkeypatch.setattr(DrawingBackend,'draw',draw)
    import topwater.render.geometry as module
    monkeypatch.setattr(module,'_submerged_edges_from_contact',lambda *args:pytest.fail('dry edges must remain visible'))
    render_result(result,tmp_path/'layers.mp4')
    masks,layers=drawn[0]
    assert [x[0] for x in layers]==['reference waterline','boat mesh','boat outline','physical waterline']
    edges,rim=mesh_edges(hull_faces)
    np.testing.assert_array_equal(layers[1][1],hull[edges,:2])
    np.testing.assert_array_equal(layers[2][1],hull[rim,:2])
    assert np.max(abs(layers[2][1]))==pytest.approx(2.)
    assert np.max(abs(np.concatenate(masks)))==pytest.approx(1.25)
    water,faces=result.surface_frame(0)
    np.testing.assert_allclose(layers[3][1],waterline_segments(hull,hull_faces,water,faces)[:,:,:2],atol=1e-12)
    assert layers[0][-1] is True and layers[3][-1] is False
    assert layers[1][2]==[176/255]*3+[1.]
    assert layers[2][2]==[0,0,0,1] and layers[3][2]==[0,0,1,1]


def test_reference_waterline_is_flat_initial_hull_fixed_once_not_first_saved_water(tmp_path, monkeypatch):
    from topwater.render.cpu.backend import DrawingBackend
    import topwater.render.geometry as module
    original=saved_result();hull,_,_,_=flared_display_case()
    water=original.surface_vertices_m.copy();water[:,2]=.2
    result=replace(original,hull_vertices_m=hull,surface_vertices_m=water,
        position_m=np.array([[0.,0.,0.],[.05,0.,.02],[.1,0.,.04]]),
        orientation_wxyz=Rotation.from_euler('x',np.array([[0.],[.025],[.05]])).as_quat(scalar_first=True))
    references,moving,computations=[],[],[]
    original_draw,original_reference=DrawingBackend.draw,module._flat_reference_waterline
    def reference(*args):
        computations.append(True);return original_reference(*args)
    def draw(self,water,faces,masks,layers):
        references.append(np.array(layers[0][1]));moving.append(np.array(layers[3][1]))
        return original_draw(self,water,faces,masks,layers)
    monkeypatch.setattr(module,'_flat_reference_waterline',reference)
    monkeypatch.setattr(DrawingBackend,'draw',draw)
    render_result(result,tmp_path/'reference.mp4')
    assert len(computations)==1 and len(references)==len(moving)==2
    np.testing.assert_array_equal(references[0],references[1])
    assert np.max(abs(references[0]))==pytest.approx(1.25)
    assert np.max(abs(moving[0]))==pytest.approx(1.4)
    assert moving[0].shape!=moving[1].shape or not np.allclose(moving[0],moving[1])


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
    import topwater.render.geometry as renderer

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
    import topwater.render.geometry as renderer

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
    from topwater.boat.geometry import physical_hull
    from topwater.geometry.waterline import waterline_segments
    import topwater.render.geometry as module
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
    import topwater.boat.geometry as geometry
    import topwater.render.geometry as module
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


def mask_samples(polygons,points,lower,upper):
    from topwater.render.cpu.backend import DrawingBackend
    backend=DrawingBackend((270,150),lower,upper)
    try:
        water=backend.vertices(backend.fullscreen,[0,0,1,1])
        pixels=backend.draw(water,np.arange(6).reshape(2,3),polygons,[])
        indices=((np.asarray(points)-lower)/(np.asarray(upper)-lower)*[270,150]).astype(int)
        return [bool(pixels[149-y,x,0]>128) for x,y in indices]
    finally:backend.close()


def test_compound_mask_preserves_concavity_and_clockwise_disconnected_fill():
    polygons=[np.array([[0.,0.],[2.,0.],[2.,1.],[0.,1.]]),
              np.array([[0.,1.],[1.,1.],[1.,2.],[0.,2.]]),
              np.array([[3.,0.],[3.,1.],[4.,1.],[4.,0.]])]
    assert mask_samples(polygons,[[.5,.5],[.5,1.5],[1.5,1.5],[3.5,.5],[2.5,.5]],
        [-.5,-.5],[4.5,2.5])==[True,True,False,True,False]


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
    # Probe the actual WebGPU stencil fill, including a nested reversed hole.
    probes=[[.5,.5],[2.,2.],[5.5,2.],[7.,2.],[7.,.5]]
    assert mask_samples(contours,probes,[-.5,-.5],[8.5,4.5])==[True,False,True,False,True]


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
