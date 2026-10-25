"""Stream saved result geometry through the shared Matplotlib presentation.

This module interpolates rigid poses, not fluid fields. Surface frames are used
exactly as saved, including their connectivity and holes. No solver or sea input
is read and no per-frame hull array is materialized for the whole sequence.
"""
from pathlib import Path
import os
import shutil
import sys
import tempfile
import time
import warnings

import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation, Slerp
from scipy.sparse import coo_array
from scipy.sparse.csgraph import connected_components

from ._visualization_2d import header_text
from ._visualization_geometry import mesh_edges, projected_footprint, waterline_segments
from .visualization import RenderSettings
from .workflow_data import Result


def sampled_poses(result):
    """Interpolate translation and shortest-path orientation at saved water times.

    HUD rotations are XYZ Euler angles of departure from the initial attitude,
    expressed in the initial horizontal heading frame. Initial trim is therefore
    not reported as wave-driven pitch. Euler labels can be nonunique at gimbal
    lock; the quaternions used for geometry remain authoritative.
    """
    times = result.surface_time_s
    positions = np.column_stack([np.interp(times, result.time_s, result.position_m[:, j])
                                 for j in range(3)])
    rotations = Slerp(result.time_s, Rotation.from_quat(result.orientation_wxyz,
                                                       scalar_first=True))(times)
    initial = Rotation.from_quat(result.orientation_wxyz[0], scalar_first=True)
    bow = initial.as_matrix()[:, 0]
    if np.linalg.norm(bow[:2]) <= 1e-12:
        raise ValueError('initial bow has no horizontal heading for the top-down HUD')
    heading = Rotation.from_euler('z', np.arctan2(bow[1], bow[0]))
    relative = heading.inv() * rotations * initial.inv() * heading
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', UserWarning)
        angles = relative.as_euler('xyz', degrees=True)
    return positions, rotations.as_matrix(), angles, bool(caught)


def _topology(faces, vertex_count):
    """Validate an oriented connected planar-sheet topology and return its loops."""
    directed = np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    owners = np.tile(np.arange(len(faces)), 3)
    ordered = np.sort(directed, axis=1)
    sort = np.lexsort((ordered[:, 1], ordered[:, 0]))
    edges = ordered[sort]
    starts = np.r_[0, 1 + np.flatnonzero(np.any(edges[1:] != edges[:-1], axis=1))]
    counts = np.diff(np.r_[starts, len(edges)])
    if np.any(counts > 2):
        raise ValueError('projected water mesh has non-manifold edges')
    interior = starts[counts == 2]
    if len(interior):
        original = directed[sort]
        if np.any(np.all(original[interior] == original[interior+1], axis=1)):
            raise ValueError('projected water mesh has inconsistent face winding')
        pair = np.column_stack((owners[sort[interior]], owners[sort[interior+1]]))
        graph = coo_array((np.ones(len(pair)), (pair[:, 0], pair[:, 1])),
                          shape=(len(faces), len(faces))).tocsr()
        if connected_components(graph, directed=False, return_labels=False) != 1:
            raise ValueError('top-down water mesh must be one edge-connected sheet')
    elif len(faces) != 1:
        raise ValueError('top-down water mesh must be one edge-connected sheet')
    boundary = directed[sort[starts[counts == 1]]]
    if not len(boundary):
        raise ValueError('top-down water sheet needs an outer boundary')
    if len(np.unique(boundary[:, 0])) != len(boundary) or \
            set(boundary[:, 0]) != set(boundary[:, 1]):
        raise ValueError('water boundaries must be disjoint simple loops')
    following = {int(a): int(b) for a, b in boundary}
    remaining, loops = set(following), []
    while remaining:
        start = min(remaining)
        loop, vertex = [], start
        while vertex in remaining:
            remaining.remove(vertex)
            loop.append(vertex)
            vertex = following[vertex]
        if vertex != start or len(loop) < 3:
            raise ValueError('water boundaries must be disjoint simple loops')
        loops.append(np.asarray(loop, dtype=np.int64))
    if vertex_count-len(starts)+len(faces) != 2-len(loops):
        raise ValueError('top-down water mesh is not a planar sheet with holes')
    # Edge manifoldness alone does not exclude two fans pinched at a vertex.
    link_edges = [[] for _ in range(vertex_count)]
    boundary_vertices = set(boundary.ravel())
    for a, b, c in faces:
        link_edges[a].append((b, c))
        link_edges[b].append((c, a))
        link_edges[c].append((a, b))
    for vertex, links in enumerate(link_edges):
        neighbors = {}
        for a, b in links:
            neighbors.setdefault(a, []).append(b)
            neighbors.setdefault(b, []).append(a)
        degrees = np.asarray([len(items) for items in neighbors.values()])
        ends = 2 if vertex in boundary_vertices else 0
        if np.any(degrees > 2) or np.count_nonzero(degrees == 1) != ends:
            raise ValueError('water mesh has a non-manifold vertex fan')
        pending, reached = [next(iter(neighbors))], set()
        while pending:
            current = pending.pop()
            if current not in reached:
                reached.add(current)
                pending.extend(neighbors[current])
        if len(reached) != len(neighbors):
            raise ValueError('water mesh has disconnected vertex fans')
    return loops


def _check_boundary_crossings(xy, loops, tolerance):
    """Reject crossing/touching boundary edges, with an X-sorted broad phase."""
    a = np.concatenate([xy[loop] for loop in loops])
    b = np.concatenate([xy[np.roll(loop, -1)] for loop in loops])
    ids_a = np.concatenate(loops)
    ids_b = np.concatenate([np.roll(loop, -1) for loop in loops])
    for loop in loops:
        u = xy[loop]-xy[np.roll(loop, 1)]
        v = xy[np.roll(loop, -1)]-xy[loop]
        turn = u[:, 0]*v[:, 1]-u[:, 1]*v[:, 0]
        if np.any((abs(turn) <= tolerance) & (np.einsum('ij,ij->i', u, v) < -tolerance)):
            raise ValueError('projected water boundary doubles back on itself')
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    order = np.argsort(lo[:, 0], kind='stable')
    for offset, i in enumerate(order):
        candidates = order[offset+1:]
        end = np.searchsorted(lo[candidates, 0], hi[i, 0]+tolerance, side='right')
        candidates = candidates[:end]
        if not len(candidates):
            continue
        candidates = candidates[(hi[candidates, 1] >= lo[i, 1]-tolerance)
                                & (lo[candidates, 1] <= hi[i, 1]+tolerance)]
        candidates = candidates[(ids_a[candidates] != ids_a[i]) & (ids_b[candidates] != ids_a[i])
                                & (ids_a[candidates] != ids_b[i]) & (ids_b[candidates] != ids_b[i])]
        if not len(candidates):
            continue
        p, q = a[candidates], b[candidates]
        u = b[i]-a[i]
        v = q-p
        cross = lambda x, y: x[..., 0]*y[..., 1]-x[..., 1]*y[..., 0]
        ca, cb = cross(u, p-a[i]), cross(u, q-a[i])
        cc, cd = cross(v, a[i]-p), cross(v, b[i]-p)
        if np.any((np.minimum(ca, cb) <= tolerance) & (np.maximum(ca, cb) >= -tolerance)
                  & (np.minimum(cc, cd) <= tolerance) & (np.maximum(cc, cd) >= -tolerance)):
            raise ValueError('projected water boundaries intersect or touch')


def _boundary_enters_view(xy, loop, lower, upper, tolerance):
    """Whether an outer edge crosses the interior of the rectangular view."""
    a, b = xy[loop], xy[np.roll(loop, -1)]
    delta = b-a
    first, last = np.zeros(len(a)), np.ones(len(a))
    valid = np.ones(len(a), dtype=bool)
    for axis in range(2):
        left, right = lower[axis]+tolerance, upper[axis]-tolerance
        parallel = abs(delta[:, axis]) <= tolerance
        valid &= ~parallel | ((a[:, axis] > left) & (a[:, axis] < right))
        lo = np.divide(left-a[:, axis], delta[:, axis], out=np.full(len(a), -np.inf), where=~parallel)
        hi = np.divide(right-a[:, axis], delta[:, axis], out=np.full(len(a), np.inf), where=~parallel)
        first = np.maximum(first, np.minimum(lo, hi))
        last = np.minimum(last, np.maximum(lo, hi))
    return bool(np.any(valid & (last > first)))


def _projected_frame(points, faces, hull, lower, upper, loops, tri_module):
    """Check a single-valued upward water sheet, including interior view holes."""
    xy = points[:, :2]
    span = max(float(np.ptp(xy, axis=0).max()), 1e-12)
    normalized = (xy-xy[0])/span
    eps = 1e-11
    if len(np.unique(xy, axis=0)) != len(xy):
        raise ValueError('top-down water projection has coincident vertices or folded layers')
    t = normalized[faces]
    u, v = t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]
    if np.any(u[:, 0]*v[:, 1]-u[:, 1]*v[:, 0] <= eps):
        raise ValueError('top-down water projection must have nondegenerate upward triangles')
    _check_boundary_crossings(normalized, loops, eps)
    signed = np.asarray([np.sum(normalized[loop, 0]*normalized[np.roll(loop, -1), 1]
                                 -normalized[loop, 1]*normalized[np.roll(loop, -1), 0])/2
                         for loop in loops])
    outer = np.flatnonzero(signed > eps)
    if len(outer) != 1 or np.any(abs(signed) <= eps):
        raise ValueError('top-down water sheet needs one outer loop and nondegenerate holes')
    envelope = ConvexHull(hull[:, :2])
    for index, loop in enumerate(loops):
        if index != outer[0] and np.any(xy[loop] @ envelope.equations[:, :2].T
                                        + envelope.equations[:, 2] > span*eps):
            raise ValueError('water mesh contains an uncovered hole outside the hull footprint')
    if _boundary_enters_view(xy, loops[outer[0]], lower, upper, span*eps):
        raise ValueError('water outer boundary enters the fixed view')
    triangulation = tri_module.Triangulation(xy[:, 0].copy(), xy[:, 1].copy(), faces)
    corners = np.array([[lower[0], lower[1]], [lower[0], upper[1]],
                        [upper[0], lower[1]], [upper[0], upper[1]]])
    try:
        found = triangulation.get_trifinder()(corners[:, 0], corners[:, 1])
    except RuntimeError as error:
        raise ValueError('invalid or overlapping top-down water triangulation') from error
    if np.any(found < 0):
        raise ValueError('water mesh does not cover the fixed view corners')
    return triangulation


def _verify_video(path, expected_frames, expected_fps):
    import imageio_ffmpeg
    reader = imageio_ffmpeg.read_frames(str(path), pix_fmt='rgb24')
    try:
        metadata = next(reader)
        # FFmpeg's human-readable header rounds fps to two decimal places.
        if tuple(metadata['size']) != (800, 800) or abs(metadata['fps']-expected_fps) > .006:
            raise RuntimeError('encoded video dimensions or cadence differ from the saved surface')
        count = sum(1 for _ in reader)
        if count != expected_frames:
            raise RuntimeError(f'encoded video has {count} frames; expected {expected_frames}')
    finally:
        reader.close()


def _publish_video(encoded, output):
    """Replace only a complete movie; cross-volume staging stays outside example."""
    output.parent.mkdir(parents=True, exist_ok=True)
    if os.stat(encoded).st_dev == os.stat(output.parent).st_dev:
        os.replace(encoded, output)
        return
    with tempfile.TemporaryDirectory(prefix='.topwater-publish-', dir=output.parent.parent) as staging:
        ready = Path(staging)/'video.mp4'
        shutil.copyfile(encoded, ready)
        with ready.open('rb') as stream:
            os.fsync(stream.fileno())
        os.replace(ready, output)


def render_result(result, output_path, settings=None):
    """Render all saved surface frames to a validated, atomically replaced MP4."""
    if not isinstance(result, Result):
        raise TypeError('render_result requires a validated workflow_data.Result')
    settings = RenderSettings() if settings is None else settings
    if not np.isfinite(settings.framing) or settings.framing <= 1:
        raise ValueError('framing must be finite and exceed one')
    if not np.isfinite(settings.maximum_seconds) or settings.maximum_seconds <= 0:
        raise ValueError('maximum_seconds must be positive and finite')
    stamps = result.surface_time_s
    cadence = stamps[1]-stamps[0]
    if not np.allclose(np.diff(stamps), cadence, atol=1e-9, rtol=1e-7):
        raise ValueError('movie needs uniformly spaced saved surface times; no fluid resampling is performed')
    output = Path(output_path)
    if output.suffix.lower() != '.mp4':
        raise ValueError('output must name a complete .mp4 file')
    started = time.monotonic()
    positions, rotations, angles, ambiguous_euler = sampled_poses(result)
    if ambiguous_euler:
        print('WARNING render: Euler HUD angles are nonunique at gimbal lock; '
              'geometry uses the saved quaternions', file=sys.stderr)
    local_hull = result.hull_vertices_m-result.cg_m
    low, high = np.full(2, np.inf), np.full(2, -np.inf)
    for position, rotation in zip(positions, rotations):
        xy = (position+local_hull @ rotation.T)[:, :2]
        low, high = np.minimum(low, xy.min(axis=0)), np.maximum(high, xy.max(axis=0))
    center = (low+high)/2
    half = float((high-low).max()/2)*settings.framing
    if half <= 0:
        raise ValueError('projected hull has no finite view extent')
    lower, upper = center-half, center+half
    magnitude = max(float(abs(result.surface_vertices_m[:, 2]).max()), 1e-12)
    edges, rim = mesh_edges(result.hull_faces)
    preparation_seconds = time.monotonic()-started

    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import colors, patheffects, pyplot as plt, tri
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Polygon
    import imageio_ffmpeg

    figure = plt.figure(figsize=(8, 8), dpi=100)
    axes = figure.add_axes((0, 0, 1, 1), aspect='equal')
    axes.set_xlim(lower[0], upper[0]); axes.set_ylim(lower[1], upper[1]); axes.set_axis_off()
    normalization = colors.TwoSlopeNorm(vmin=-magnitude, vcenter=0., vmax=magnitude)
    hull = positions[0]+local_hull @ rotations[0].T
    interior = Polygon(projected_footprint(hull), closed=True, facecolor='#f3f1e8',
                       edgecolor='none', alpha=1., zorder=1)
    axes.add_patch(interior)
    mesh_artist = LineCollection([], colors='#b0b0b0', linewidths=.7, zorder=2)
    rim_artist = LineCollection([], colors='black', linewidths=1.5, zorder=3)
    contact_artist = LineCollection([], colors='#003a70', linewidths=1.3, zorder=4)
    contact_artist.set_path_effects([patheffects.Stroke(linewidth=2.6, foreground='white'), patheffects.Normal()])
    for artist in (mesh_artist, rim_artist, contact_artist):
        axes.add_collection(artist)
    header = figure.text(.5, .985, '', ha='center', va='top', fontsize=14.5,
                         fontfamily='DejaVu Sans Mono', color='black',
                         bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .85, 'pad': 3})
    for artist in (interior, mesh_artist, rim_artist, contact_artist, header):
        artist.set_animated(True)
    figure.canvas.draw()
    background = figure.canvas.copy_from_bbox(figure.bbox)
    setup_seconds = time.monotonic()-started-preparation_seconds
    geometry_seconds = draw_seconds = pipe_seconds = 0.
    topology_changes = 0
    water = previous_faces = previous_count = loops = drawn_triangulation = None
    fps = float(1/cadence)
    try:
        with tempfile.TemporaryDirectory(prefix='topwater-render-') as temporary:
            encoded = Path(temporary)/'video.mp4'
            writer = imageio_ffmpeg.write_frames(str(encoded), (800, 800), fps=fps,
                codec='libx264', pix_fmt_in='rgba', pix_fmt_out='yuv420p', bitrate='3000k',
                input_params=['-r', format(fps, '.17g')],
                ffmpeg_timeout=min(30., settings.maximum_seconds),
                output_params=['-threads', '1'])
            encoder_start = time.monotonic()
            writer.send(None)
            encoder_start_seconds = time.monotonic()-encoder_start
            try:
                for index, stamp in enumerate(stamps):
                    before = time.monotonic()
                    if before-started > settings.maximum_seconds:
                        raise TimeoutError('render exceeded maximum_seconds')
                    points, faces = result.surface_frame(index)
                    hull = positions[index]+local_hull @ rotations[index].T
                    changed = previous_count != len(points) or previous_faces is None \
                        or not np.array_equal(previous_faces, faces)
                    if changed:
                        loops = _topology(faces, len(points))
                    projected = _projected_frame(points, faces, hull, lower, upper, loops, tri)
                    if changed:
                        if water is not None:
                            water.remove()
                            topology_changes += 1
                        water = axes.tripcolor(projected, points[:, 2], shading='gouraud',
                            cmap='coolwarm', norm=normalization, edgecolors='none', linewidth=0, zorder=0)
                        water.set_animated(True)
                        drawn_triangulation = projected
                        previous_faces, previous_count = faces, len(points)
                    else:
                        drawn_triangulation.x[:] = points[:, 0]
                        drawn_triangulation.y[:] = points[:, 1]
                        water.set_array(points[:, 2])
                    interior.set_xy(projected_footprint(hull))
                    mesh_artist.set_segments(hull[edges, :2])
                    rim_artist.set_segments(hull[rim, :2])
                    contact = waterline_segments(hull, result.hull_faces, points, faces)
                    contact_artist.set_segments(contact[:, :, :2])
                    header.set_text(header_text(stamp, angles[index, 0], angles[index, 1],
                                                positions[index, 2]-positions[0, 2]))
                    geometry_seconds += time.monotonic()-before
                    drawing = time.monotonic()
                    figure.canvas.restore_region(background)
                    for artist in (water, interior, mesh_artist, rim_artist, contact_artist):
                        axes.draw_artist(artist)
                    figure.draw_artist(header)
                    figure.canvas.blit(figure.bbox)
                    draw_seconds += time.monotonic()-drawing
                    sending = time.monotonic()
                    writer.send(np.asarray(figure.canvas.buffer_rgba()))
                    pipe_seconds += time.monotonic()-sending
            finally:
                closing = time.monotonic()
                writer.close()
                encoder_close_seconds = time.monotonic()-closing
            verifying = time.monotonic()
            _verify_video(encoded, len(stamps), fps)
            verification_seconds = time.monotonic()-verifying
            if time.monotonic()-started > settings.maximum_seconds:
                raise TimeoutError('render exceeded maximum_seconds')
            _publish_video(encoded, output)
        return dict(complete=True, output=str(output), frames=len(stamps), fps=fps,
                    size=(800, 800), source_end_s=float(stamps[-1]),
                    video_duration_s=len(stamps)/fps, preparation_seconds=preparation_seconds,
                    setup_seconds=setup_seconds, geometry_seconds=geometry_seconds,
                    draw_seconds=draw_seconds, encoder_send_seconds=pipe_seconds,
                    encoder_start_seconds=encoder_start_seconds, encoder_close_seconds=encoder_close_seconds,
                    verification_seconds=verification_seconds, elapsed_s=time.monotonic()-started,
                    topology_changes=topology_changes, ambiguous_euler_labels=ambiguous_euler,
                    matplotlib_version=matplotlib.__version__, backend='Matplotlib Agg CPU blit')
    finally:
        plt.close(figure)
