"""Minimal top-down scalar field and mesh annotations using standard Matplotlib."""
from dataclasses import asdict
from pathlib import Path
import shutil
import tempfile
import time

import numpy as np

from ._visualization_geometry import mesh_edges, waterline_segments, projected_footprint


def header_text(stamp, roll, pitch, heave):
    """Fixed-width columns for a monospace motion header."""
    return (f'TIME {stamp:7.2f} s   ROLL {roll:+7.2f}°   '
            f'PITCH {pitch:+7.2f}°   HEAVE {heave:+7.2f} m')


def render_topdown(scene, output_path, settings):
    from .visualization import convex_envelope

    if len(scene.time_s) < 2 or not np.allclose(np.diff(scene.time_s), scene.time_s[1] - scene.time_s[0], atol=1e-9, rtol=1e-7):
        raise ValueError('Movie input needs at least two uniformly spaced times; resampling belongs to the producer adapter')
    output = Path(output_path)
    if output.suffix.lower() != '.mp4':
        raise ValueError('Output must name the complete MP4 file, including .mp4')
    if output.exists():
        raise FileExistsError(f'Output already exists: {output}')
    required = ('roll_deg', 'pitch_deg', 'heave_m')
    if any(name not in scene.telemetry for name in required):
        raise ValueError('Top-down header requires producer telemetry: roll_deg, pitch_deg, heave_m')
    if settings.framing <= 1:
        raise ValueError('Framing must exceed one')
    started = time.monotonic()
    for points in scene.hull_points:
        convex_envelope(points, scene.hull_faces)
    validation_seconds = time.monotonic() - started
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import colors, patheffects, pyplot as plt, tri
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Polygon
    import imageio_ffmpeg

    fps = float(1 / (scene.time_s[1] - scene.time_s[0]))
    center = ((scene.hull_points.min(axis=1) + scene.hull_points.max(axis=1)) / 2).mean(axis=0)[:2]
    half = float(abs(scene.hull_points[:, :, :2] - center).max()) * settings.framing
    lower, upper = center - half, center + half
    # Cull whole triangles only when their all-time XY bounds miss the fixed view.
    point_lo = scene.water_points[:, :, :2].min(axis=0)
    point_hi = scene.water_points[:, :, :2].max(axis=0)
    face_lo, face_hi = point_lo[scene.water_faces].min(axis=1), point_hi[scene.water_faces].max(axis=1)
    mask = np.all((face_hi >= lower) & (face_lo <= upper), axis=1)
    faces = scene.water_faces[mask]
    if not len(faces):
        raise ValueError('No water triangles cover the requested view')
    nodes = np.unique(faces)
    local_faces = np.searchsorted(nodes, faces)
    edges, rim = mesh_edges(scene.hull_faces)
    magnitude = max(float(abs(scene.water_points[:, :, 2]).max()), 1e-12)
    normalization = colors.TwoSlopeNorm(vmin=-magnitude, vcenter=0., vmax=magnitude)
    figure = plt.figure(figsize=(8, 8), dpi=100)
    axes = figure.add_axes((0, 0, 1, 1), aspect='equal')
    axes.set_xlim(lower[0], upper[0]); axes.set_ylim(lower[1], upper[1]); axes.set_axis_off()
    points = scene.water_points[0, nodes]
    triangulation = tri.Triangulation(points[:, 0].copy(), points[:, 1].copy(), local_faces)
    water = axes.tripcolor(triangulation, points[:, 2], shading='gouraud', cmap='coolwarm', norm=normalization, edgecolors='none', linewidth=0, zorder=0)
    interior_artist = Polygon(projected_footprint(scene.hull_points[0]), closed=True,
                              facecolor='#f3f1e8', edgecolor='none', alpha=1., zorder=1)
    axes.add_patch(interior_artist)
    mesh_artist = LineCollection([], colors='#b0b0b0', linewidths=.7, zorder=2)
    rim_artist = LineCollection([], colors='black', linewidths=1.5, zorder=3)
    contact_artist = LineCollection([], colors='#003a70', linewidths=1.3, zorder=4)
    contact_artist.set_path_effects([patheffects.Stroke(linewidth=2.6, foreground='white'), patheffects.Normal()])
    for artist in (mesh_artist, rim_artist, contact_artist):
        axes.add_collection(artist)
    header = figure.text(.5, .985, '', ha='center', va='top', fontsize=14.5, fontfamily='DejaVu Sans Mono', color='black',
                         bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .85, 'pad': 3})
    frame_times, geometry_times, draw_times, pipe_times, segment_counts = [], [], [], [], []
    corners = np.array([[lower[0], lower[1]], [lower[0], upper[1]], [upper[0], lower[1]], [upper[0], upper[1]]])
    for artist in (water, interior_artist, mesh_artist, rim_artist, contact_artist, header):
        artist.set_animated(True)
    figure.canvas.draw()
    background = figure.canvas.copy_from_bbox(figure.bbox)
    setup_seconds = time.monotonic() - started - validation_seconds
    try:
        with tempfile.TemporaryDirectory(prefix='topwater-render-') as temporary:
            encoded = Path(temporary) / 'movie.mp4'
            writer = imageio_ffmpeg.write_frames(str(encoded), (800, 800), fps=fps, codec='libx264',
                pix_fmt_in='rgba', pix_fmt_out='yuv420p', bitrate='3000k', output_params=['-threads', '1'])
            encoder_start = time.monotonic()
            writer.send(None)
            encoder_startup_seconds = time.monotonic() - encoder_start
            try:
                    for index, stamp in enumerate(scene.time_s):
                        before = time.monotonic()
                        if before - started > settings.maximum_seconds:
                            raise TimeoutError('Render exceeded maximum_seconds')
                        points = scene.water_points[index, nodes]
                        coverage = tri.Triangulation(points[:, 0], points[:, 1], local_faces)
                        if np.any(coverage.get_trifinder()(corners[:, 0], corners[:, 1]) < 0):
                            raise ValueError('Water mesh does not cover the fixed view corners')
                        triangulation.x[:] = points[:, 0]; triangulation.y[:] = points[:, 1]
                        water.set_array(points[:, 2])
                        hull = scene.hull_points[index]
                        interior_artist.set_xy(projected_footprint(hull))
                        mesh_artist.set_segments(hull[edges, :2])
                        rim_artist.set_segments(hull[rim, :2])
                        contact = waterline_segments(hull, scene.hull_faces, scene.water_points[index], faces)
                        contact_artist.set_segments(contact[:, :, :2])
                        segment_counts.append(len(contact))
                        header.set_text(header_text(stamp, scene.telemetry['roll_deg'][index],
                            scene.telemetry['pitch_deg'][index], scene.telemetry['heave_m'][index]))
                        geometry_times.append(time.monotonic() - before)
                        drawing = time.monotonic()
                        figure.canvas.restore_region(background)
                        for artist in (water, interior_artist, mesh_artist, rim_artist, contact_artist):
                            axes.draw_artist(artist)
                        figure.draw_artist(header)
                        figure.canvas.blit(figure.bbox)
                        draw_times.append(time.monotonic() - drawing)
                        sending = time.monotonic()
                        writer.send(np.asarray(figure.canvas.buffer_rgba()))
                        pipe_times.append(time.monotonic() - sending)
                        frame_times.append(time.monotonic() - before)
            finally:
                closing = time.monotonic()
                writer.close()
                encoder_close_seconds = time.monotonic() - closing
            output.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(encoded), str(output))
        return dict(complete=True, output=str(output), frames=len(scene.time_s), fps=fps, size=(800, 800),
            source_start_s=float(scene.time_s[0]), source_end_s=float(scene.time_s[-1]), video_duration_s=len(scene.time_s) / fps,
            elapsed_s=time.monotonic() - started, validation_seconds=validation_seconds, context_setup_seconds=setup_seconds,
            frame_seconds=frame_times, geometry_seconds=geometry_times, agg_draw_seconds=draw_times,
            encoder_send_seconds=pipe_times, encoder_startup_seconds=encoder_startup_seconds,
            encoder_close_seconds=encoder_close_seconds,
            encoder_timing_note='send is pipe transfer/backpressure, not isolated encoder CPU time; encoding overlaps plotting',
            visible_water_faces=len(faces), total_water_faces=len(scene.water_faces), contact_segments=segment_counts,
            all_view_corners_have_water=True, height_normalization_m=(-magnitude, 0., magnitude),
            settings=asdict(settings), metadata=scene.metadata, matplotlib_version=matplotlib.__version__,
            backend_capabilities='Matplotlib Agg CPU blit; public imageio-ffmpeg frame stream, no second full-figure draw',
            limitations='Top-down projected hull annotations and actual triangle/triangle contact; no flooding interpretation or height amplification.')
    finally:
        plt.close(figure)
