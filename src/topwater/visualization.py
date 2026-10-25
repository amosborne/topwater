"""Producer-neutral mesh animation using optional PyVista/VTK rendering.

Inputs are already prepared world geometry, never simulation parameters or loads.
One sequence has fixed triangular connectivity; another hull is another input.
"""
from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
import shutil
import tempfile
import time

import numpy as np
from scipy.spatial import ConvexHull, QhullError, cKDTree


@dataclass
class SceneSequence:
    """World meshes in metres, right handed with z up, at increasing times.

    ``*_points`` have shape (time, vertex, 3); ``*_faces`` are (triangle, 3)
    integer indices. Hull faces may leave a gunwale open, but must lie on the
    convex envelope. The closed envelope is used only to mask displayed water.
    Water faces should be consistently oriented toward the air.
    """
    time_s: np.ndarray
    hull_points: np.ndarray
    hull_faces: np.ndarray
    water_points: np.ndarray
    water_faces: np.ndarray
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.time_s = np.asarray(self.time_s, dtype=float)
        if self.time_s.ndim != 1 or not len(self.time_s) or not np.isfinite(self.time_s).all():
            raise ValueError('time_s must be a nonempty finite vector')
        if np.any(np.diff(self.time_s) <= 0):
            raise ValueError('time_s must increase strictly')
        for name in ('hull', 'water'):
            points = np.asarray(getattr(self, name + '_points'), dtype=float)
            raw_faces = np.asarray(getattr(self, name + '_faces'))
            if points.ndim != 3 or points.shape[0] != len(self.time_s) or points.shape[2] != 3 or not np.isfinite(points).all():
                raise ValueError(f'{name}_points must be finite (time, vertex, 3)')
            if raw_faces.ndim != 2 or raw_faces.shape[1] != 3 or not len(raw_faces) or not np.issubdtype(raw_faces.dtype, np.integer):
                raise ValueError(f'{name}_faces must be integer (triangle, 3)')
            if raw_faces.min() < 0 or raw_faces.max() >= points.shape[1]:
                raise ValueError(f'{name}_faces index outside vertex array')
            setattr(self, name + '_points', points)
            setattr(self, name + '_faces', raw_faces.astype(np.int64, copy=False))
        if self.metadata.get('units', 'm') != 'm' or self.metadata.get('coordinate_system', 'right-handed-z-up') != 'right-handed-z-up':
            raise ValueError('Normalize geometry to metres and right-handed-z-up before rendering')

    def save(self, path):
        """Save prepared geometry once; style changes can reuse it without physics."""
        np.savez(path, time_s=self.time_s, hull_points=self.hull_points,
                 hull_faces=self.hull_faces, water_points=self.water_points,
                 water_faces=self.water_faces, metadata=json.dumps(self.metadata))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as data:
            return cls(**{key: data[key] for key in ('time_s', 'hull_points', 'hull_faces', 'water_points', 'water_faces')},
                       metadata=json.loads(str(data['metadata'])))


def convex_envelope(points, faces):
    """Validate convex boundary faces and return outward oriented closure faces.

    An internal face is rejected rather than silently replacing a nonconvex boat
    by its convex hull. Open boundaries are permitted for the visual hull.
    """
    points = np.asarray(points, dtype=float)
    faces = np.asarray(faces, dtype=int)
    scale = np.linalg.norm(np.ptp(points, axis=0))
    try:
        envelope = ConvexHull(points)
    except QhullError as error:
        raise ValueError('Hull must enclose a three-dimensional convex volume') from error
    if np.any(np.linalg.norm(np.cross(points[faces[:, 1]] - points[faces[:, 0]], points[faces[:, 2]] - points[faces[:, 0]]), axis=1) <= 1e-12 * scale**2):
        raise ValueError('Hull contains degenerate triangles')
    oriented = _outward(points, faces)
    triangle = points[oriented]
    normal = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    # A boundary triangle's outward plane must match a convex-envelope plane.
    # Query by normal, then test its vertices against that support plane. This
    # avoids an O(faces * envelope_planes) dense allocation on refined meshes.
    _, nearest = cKDTree(envelope.equations[:, :3]).query(normal)
    planes = envelope.equations[nearest]
    residual = abs(np.einsum('fvi,fi->fv', triangle, planes[:, :3]) + planes[:, 3, None])
    if np.any(residual.max(axis=1) > 1e-8 * max(scale, 1e-12)):
        raise ValueError('Hull faces must lie on a convex boundary')
    closure = envelope.simplices.copy()
    return _outward(points, closure), oriented


def _outward(points, faces):
    result = faces.copy()
    triangle = points[result]
    normal = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    flip = np.einsum('ij,ij->i', normal, triangle.mean(axis=1) - points.mean(axis=0)) < 0
    result[flip] = result[flip][:, [0, 2, 1]]
    return result


@dataclass
class RenderSettings:
    """Scene styling; none of these settings belong to a particular producer."""
    azimuth_deg: float = -110.
    elevation_deg: float = 24.
    view_angle_deg: float = 25.
    framing: float = 1.45
    target_z: float = 0.
    follow_xy: bool = True
    perspective: bool = True
    hull_color: str = '#dedbd1'
    interior_color: str = '#ae8053'
    hull_opacity: float = 1.
    water_opacity: float = .9
    water_colormap: str = 'Blues_r'
    water_specular: float = .35
    water_specular_power: float = 35.
    shadows: bool = False
    waterline: bool = False
    timestamps: bool = True
    maximum_seconds: float = 240.


def render_sequence(scene, output_path, settings=None):
    """Render prepared meshes with one reused VTK context and actors.

    Render the whole input at its saved uniform cadence to one 800-square MP4.
    Diagnostics are returned, never written as persistent sidecars. VTK filters
    perform clipping and normals; this is visualization, not wetting physics.
    """
    settings = RenderSettings() if settings is None else settings
    if settings.framing <= 1 or not (0 < settings.view_angle_deg < 90) or not (0 < settings.elevation_deg < 90):
        raise ValueError('Framing > 1 and acute camera angles required')
    if len(scene.time_s) < 2 or not np.allclose(np.diff(scene.time_s), scene.time_s[1] - scene.time_s[0], atol=1e-9, rtol=1e-7):
        raise ValueError('Movie input needs at least two uniformly spaced times; resampling belongs to the producer adapter')
    fps = 1 / (scene.time_s[1] - scene.time_s[0])
    output = Path(output_path)
    if output.suffix.lower() != '.mp4':
        raise ValueError('Output must name the complete MP4 file, including .mp4')
    if output.exists():
        raise FileExistsError(f'Output already exists: {output}')
    started = time.monotonic()
    closures = [convex_envelope(points, scene.hull_faces) for points in scene.hull_points]
    validation_seconds = time.monotonic() - started
    try:
        import pyvista as pv
        import imageio.v2 as imageio
        from vtkmodules.vtkRenderingCore import vtkCellPicker
    except ImportError as error:
        raise ImportError('Install topwater[visualization] to render scenes') from error

    def cells(faces):
        return np.column_stack((np.full(len(faces), 3), faces)).ravel()

    hull = pv.PolyData(scene.hull_points[0].copy(), cells(closures[0][1]))
    water = pv.PolyData(scene.water_points[0].copy(), cells(scene.water_faces))
    water['elevation'] = water.points[:, 2]
    closed = pv.PolyData(scene.hull_points[0].copy(), cells(closures[0][0]))
    clipped = water.clip_surface(closed, invert=False).extract_surface(algorithm='dataset_surface')
    clipped = clipped.compute_normals(point_normals=True, cell_normals=False, consistent_normals=True)
    size = (800, 800)
    plotter = pv.Plotter(off_screen=True, window_size=size, lighting='none')
    temporary = tempfile.TemporaryDirectory(prefix='topwater-render-')
    encoded = Path(temporary.name) / 'movie.mp4'
    writer = None
    frame_seconds = []
    coverage = []
    try:
        plotter.set_background('#4384a2')
        plotter.enable_anti_aliasing('fxaa')
        plotter.enable_depth_peeling(number_of_peels=8)
        hull_actor = plotter.add_mesh(hull, color=settings.hull_color, opacity=settings.hull_opacity,
            show_edges=False, smooth_shading=False, ambient=.35, diffuse=.65, specular=.12,
            backface_params={'color': settings.interior_color, 'ambient': .4, 'diffuse': .6})
        limits = (float(scene.water_points[:, :, 2].min()), float(scene.water_points[:, :, 2].max()))
        if limits[1] - limits[0] < 1e-8:
            limits = (limits[0] - .01, limits[1] + .01)
        water_actor = plotter.add_mesh(clipped, scalars='elevation', cmap=settings.water_colormap,
            clim=limits, opacity=settings.water_opacity, show_edges=False, show_scalar_bar=False,
            smooth_shading=True, ambient=.35, diffuse=.65, specular=settings.water_specular,
            specular_power=settings.water_specular_power)
        line_actor = None
        if settings.waterline:
            edges = clipped.extract_feature_edges(boundary_edges=True, feature_edges=False,
                manifold_edges=False, non_manifold_edges=False)
            line_actor = plotter.add_mesh(edges, color='#163e50', line_width=2, lighting=False)
        centers = (scene.hull_points.min(axis=1) + scene.hull_points.max(axis=1)) / 2
        centers[:, 2] = settings.target_z
        if not settings.follow_xy:
            centers[:] = centers[0]
        radius = float(np.linalg.norm(scene.hull_points - centers[:, None], axis=2).max())
        span = 2 * radius * settings.framing
        distance = span / (2 * np.tan(np.deg2rad(settings.view_angle_deg / 2)))
        elevation, azimuth = np.deg2rad([settings.elevation_deg, settings.azimuth_deg])
        toward = np.array([np.cos(elevation) * np.cos(azimuth), np.cos(elevation) * np.sin(azimuth), np.sin(elevation)])
        center0 = centers[0].copy()
        sunlight = pv.Light(position=center0 + span * np.array([-.6, -.5, 1.5]), focal_point=center0,
                            color='white', intensity=.7, light_type='scene light')
        plotter.add_light(sunlight)
        plotter.add_light(pv.Light(light_type='headlight', intensity=.35))
        if settings.shadows:
            plotter.enable_shadows()
        plotter.camera.view_angle = settings.view_angle_deg
        plotter.camera.parallel_projection = not settings.perspective
        plotter.camera.parallel_scale = span / 2
        plotter.camera.up = (0, 0, 1)
        picker = vtkCellPicker()
        picker.PickFromListOn()
        picker.AddPickList(water_actor)
        picker.SetTolerance(1e-6)
        plotter.camera_position = (center0 + toward * distance, center0, (0, 0, 1))
        plotter.show(auto_close=False, interactive=False)
        capabilities = plotter.render_window.ReportCapabilities()
        writer = imageio.get_writer(encoded, fps=fps, codec='libx264',
            quality=8, macro_block_size=None, ffmpeg_params=['-threads', '1'])
        setup_seconds = time.monotonic() - started - validation_seconds
        for index, stamp in enumerate(scene.time_s):
            before = time.monotonic()
            if before - started > settings.maximum_seconds:
                raise TimeoutError('Render exceeded maximum_seconds')
            hull.points = scene.hull_points[index]
            hull.faces = cells(closures[index][1])
            closed.points = scene.hull_points[index]
            closed.faces = cells(closures[index][0])
            water.points = scene.water_points[index]
            water['elevation'] = water.points[:, 2]
            clipped = water.clip_surface(closed, invert=False).extract_surface(algorithm='dataset_surface')
            clipped = clipped.compute_normals(point_normals=True, cell_normals=False, consistent_normals=True)
            water_actor.mapper.dataset = clipped
            if line_actor is not None:
                line_actor.mapper.dataset = clipped.extract_feature_edges(boundary_edges=True,
                    feature_edges=False, manifold_edges=False, non_manifold_edges=False)
            center = centers[index]
            plotter.camera_position = (center + toward * distance, center, (0, 0, 1))
            plotter.reset_camera_clipping_range()
            sunlight.position = center + span * np.array([-.6, -.5, 1.5])
            sunlight.focal_point = center
            if settings.timestamps:
                plotter.add_text(f't = {stamp:.2f} s', position='upper_left', font_size=12,
                                 color='white', shadow=True, name='time')
            plotter.render()
            width, height = size
            corners = [(1, 1), (width - 2, 1), (1, height - 2), (width - 2, height - 2)]
            hits = [bool(picker.Pick(x, y, 0, plotter.renderer)) for x, y in corners]
            if not all(hits):
                raise ValueError('Water mesh does not cover the camera corners; extend the prepared water surface')
            coverage.append(hits)
            pixels = plotter.screenshot(return_img=True)
            if writer is not None:
                writer.append_data(pixels)
            frame_seconds.append(time.monotonic() - before)
        if writer is not None:
            writer.close()
            writer = None
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(encoded), str(output))
        report = dict(complete=True, output=str(output), frames=len(scene.time_s), fps=float(fps), size=size,
            source_start_s=float(scene.time_s[0]), source_end_s=float(scene.time_s[-1]),
            video_duration_s=float(len(scene.time_s) / fps),
            elapsed_s=time.monotonic() - started, validation_seconds=validation_seconds,
            context_setup_seconds=setup_seconds, frame_seconds=frame_seconds,
            all_camera_corners_hit_water=bool(np.all(coverage)), settings=asdict(settings),
            metadata=scene.metadata, pyvista_version=pv.__version__, vtk_version=pv.vtk_version_info,
            backend_capabilities=capabilities,
            limitations='Convex display-volume water masking; not physical wetting/flooding. Fixed mesh connectivity per sequence.')
        return report
    finally:
        try:
            if writer is not None:
                writer.close()
        finally:
            try:
                plotter.close()
            finally:
                temporary.cleanup()
