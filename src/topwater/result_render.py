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
from scipy.spatial import ConvexHull, cKDTree
from scipy.spatial.transform import Rotation, Slerp
from scipy.sparse import coo_array
from scipy.sparse.csgraph import connected_components

from ._visualization_2d import header_text
from ._visualization_geometry import mesh_edges, waterline_segments
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
    with np.errstate(over='ignore', invalid='ignore'):
        span = max(float(np.ptp(xy, axis=0).max()), 1e-12)
    if not np.isfinite(span):
        raise ValueError('top-down water projection has a nonfinite computed span')
    normalized = (xy-xy[0])/span
    eps = 1e-11
    if len(np.unique(xy, axis=0)) != len(xy):
        raise ValueError('top-down water projection has coincident vertices or folded layers')
    # Projection degeneracy is local shape, not area relative to the full patch.
    # Form edges before global normalization so small valid cut faces survive.
    t = xy[faces]
    with np.errstate(over='ignore', invalid='ignore'):
        u, v = t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]
        determinant = u[:, 0]*v[:, 1]-u[:, 1]*v[:, 0]
        threshold = 64*np.finfo(float).eps*np.linalg.norm(u, axis=1)*np.linalg.norm(v, axis=1)
    if not np.isfinite(determinant).all() or not np.isfinite(threshold).all():
        raise ValueError('top-down water projection has nonfinite area computations')
    if np.any(determinant <= threshold):
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


def _contact_contours(contact, scale):
    """Join the existing oriented physical contact; never intersect new meshes."""
    from matplotlib.path import Path as ArtistPath
    if not len(contact):
        return []
    endpoints = contact.reshape(-1, 3)
    tolerance = scale*1e-9
    parent = np.arange(len(endpoints))

    def representative(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for a, b in sorted(cKDTree(endpoints).query_pairs(tolerance)):
        a, b = representative(a), representative(b)
        parent[max(a, b)] = min(a, b)
    nodes = np.array([representative(i) for i in range(len(endpoints))])
    directed = set()
    for a, b in nodes.reshape(-1, 2):
        if a != b:
            if (b, a) in directed:
                raise ValueError('waterline contact has inconsistent orientation')
            directed.add((int(a), int(b)))
    following, incoming = {}, {}
    for a, b in directed:
        if a in following or b in incoming:
            raise ValueError('waterline contact is branched; transparent hull needs simple closed contours')
        following[a], incoming[b] = b, a
    if not following or set(following) != set(incoming):
        raise ValueError('waterline contact is open or collapsed; transparent hull needs closed contours')
    remaining, loops = set(following), []
    while remaining:
        start = min(remaining)
        loop, node = [], start
        while node in remaining:
            remaining.remove(node)
            loop.append(node)
            node = following[node]
        if node != start or len(loop) < 3:
            raise ValueError('waterline contact does not form simple closed contours')
        loops.append(np.array(loop))
    xy = (endpoints[:, :2]-endpoints[:, :2].mean(axis=0))/scale
    _check_boundary_crossings(xy, loops, 1e-11)
    paths = [ArtistPath(np.vstack((xy[loop], xy[loop[0]])), closed=True) for loop in loops]
    for index, loop in enumerate(loops):
        polygon = xy[loop]
        area = np.sum(polygon[:, 0]*np.roll(polygon[:, 1], -1)
                      -polygon[:, 1]*np.roll(polygon[:, 0], -1))
        depth = sum(path.contains_point(polygon[0]) for j, path in enumerate(paths) if j != index)
        # cross(outward hull normal, upward water normal) puts solid on right.
        if not np.isfinite(area) or area == 0. or (area < 0.) != (depth % 2 == 0):
            raise ValueError('waterline projection has ambiguous area or inconsistent nested winding')
    return [endpoints[loop, :2] for loop in loops]


def _submerged_edges_from_contact(hull, edges, contact, projected, points, faces, tolerance,
                                  coverage_loops=()):
    """Split at existing contact/coverage boundaries; classify P1 midpoints."""
    endpoints = contact.reshape(-1, 3)
    boundary = np.concatenate([np.stack((points[loop, :2], points[np.roll(loop, -1), :2]), axis=1)
                               for loop in coverage_loops]) if len(coverage_loops) else np.empty((0, 2, 2))
    if len(boundary):
        selected = np.all(boundary.max(axis=1) >= hull[:, :2].min(axis=0)-tolerance, axis=1) \
            & np.all(boundary.min(axis=1) <= hull[:, :2].max(axis=0)+tolerance, axis=1)
        boundary = boundary[selected]
    segments, owners = [], []
    for owner, (start, end) in enumerate(hull[edges]):
        delta = end-start
        length = np.linalg.norm(delta)
        parameter = (endpoints-start) @ delta/length**2
        distance = np.linalg.norm(endpoints-start-parameter[:, None]*delta, axis=1)
        selected = (distance <= tolerance) & (parameter >= 0.) & (parameter <= 1.)
        coverage_cuts = []
        horizontal_length = np.linalg.norm(delta[:2])
        if len(boundary) and horizontal_length > tolerance:
            direction = boundary[:, 1]-boundary[:, 0]
            offset = boundary[:, 0]-start[:2]
            boundary_length = np.linalg.norm(direction, axis=1)
            cross = lambda a, b: a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]
            determinant = cross(delta[:2], direction)
            parallel = abs(determinant) <= 64*np.finfo(float).eps*horizontal_length*boundary_length
            at = np.divide(cross(offset, direction), determinant,
                           out=np.zeros(len(boundary)), where=~parallel)
            along = np.divide(cross(offset, delta[:2]), determinant,
                              out=np.zeros(len(boundary)), where=~parallel)
            crossing = ~parallel & (at >= -tolerance/horizontal_length) \
                & (at <= 1.+tolerance/horizontal_length) \
                & (along >= -tolerance/boundary_length) & (along <= 1.+tolerance/boundary_length)
            coverage_cuts.extend(np.clip(at[crossing], 0., 1.))
            collinear = parallel & (abs(cross(offset, delta[:2])) <= tolerance*horizontal_length)
            parameters = ((boundary[collinear]-start[:2]) @ delta[:2]/horizontal_length**2).ravel()
            coverage_cuts.extend(parameters[(parameters >= 0.) & (parameters <= 1.)])
        cuts = np.sort(np.r_[0., parameter[selected], coverage_cuts, 1.])
        cuts = cuts[np.r_[True, np.diff(cuts) > tolerance/length]]
        if cuts[-1] < 1.:
            cuts[-1] = 1.
        for left, right in zip(cuts[:-1], cuts[1:]):
            segments.append([start+left*delta, start+right*delta])
            owners.append(owner)
    segments = np.asarray(segments).reshape(-1, 2, 3)
    owners = np.asarray(owners, dtype=int)
    middle = segments.mean(axis=1)
    found = projected.get_trifinder()(middle[:, 0], middle[:, 1])
    covered = found >= 0
    triangles = points[faces[found[covered]]]
    a, b = triangles[:, 1, :2]-triangles[:, 0, :2], triangles[:, 2, :2]-triangles[:, 0, :2]
    offset = middle[covered, :2]-triangles[:, 0, :2]
    determinant = a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]
    u = (offset[:, 0]*b[:, 1]-offset[:, 1]*b[:, 0])/determinant
    v = (a[:, 0]*offset[:, 1]-a[:, 1]*offset[:, 0])/determinant
    elevation = (1-u-v)*triangles[:, 0, 2]+u*triangles[:, 1, 2]+v*triangles[:, 2, 2]
    wet = np.zeros(len(segments), dtype=bool)
    wet[covered] = middle[covered, 2] <= elevation+tolerance
    return segments[wet], owners[wet]


def _wetted_display(hull, hull_faces, points, faces, loops, timings=None, projected=None,
                    compute_submerged_edges=True):
    """Reuse contact for the mask and optional submerged-only diagnostic edges.

    Rendering shows the full wireframe, so it skips submerged edge computation.
    """
    scale = max(float(np.linalg.norm(np.ptp(hull, axis=0))), 1e-12)
    tolerance = scale*1e-9
    started = time.monotonic()
    try:
        contact = waterline_segments(hull, hull_faces, points, faces)
    finally:
        if timings is not None:
            timings['contact_seconds'] += time.monotonic()-started
    started = time.monotonic()
    # Legacy omitted-water holes must still be hidden by solid, not dry hull.
    holes = []
    for loop in loops:
        xy = points[loop, :2]
        if np.sum(xy[:, 0]*np.roll(xy[:, 1], -1)-xy[:, 1]*np.roll(xy[:, 0], -1)) < 0.:
            holes.append(loop)
    if holes:
        planes = ConvexHull(hull).equations
        for loop in holes:
            if np.any(points[loop] @ planes[:, :3].T+planes[:, 3] > tolerance):
                raise ValueError('water hole is exposed by transparent dry hull; saved water is missing')
    try:
        contours = _contact_contours(contact, scale)
    except ValueError as error:
        if not str(error).startswith('waterline contact is open or collapsed'):
            raise
        # A surface crossing the physical opening leaves genuine physical
        # contact open. Complete ONLY the display-volume boundary with the
        # existing internal closure; it is never returned as physical contact.
        # Dry closed-contact frames take no closure/intersection fallback.
        from .boat_geometry import physical_hull
        closed = physical_hull(hull, hull_faces)
        closure_trace = waterline_segments(closed.vertices, closed.faces[closed.artificial_faces],
                                           points, faces)
        contours = _contact_contours(np.concatenate((contact, closure_trace)), scale)
    if timings is not None:
        timings['mask_seconds'] += time.monotonic()-started
    if not compute_submerged_edges:
        empty = np.empty((0, 2, 3))
        return contours, empty, empty, contact
    started = time.monotonic()
    if projected is None:
        from matplotlib.tri import Triangulation
        projected = Triangulation(points[:, 0], points[:, 1], faces)
    edges, rim = mesh_edges(hull_faces)
    submerged, owners = _submerged_edges_from_contact(
        hull, edges, contact, projected, points, faces, tolerance, loops)
    rim_set = {tuple(edge) for edge in rim}
    on_rim = np.array([tuple(edge) in rim_set for edge in edges])
    submerged_rim = submerged[on_rim[owners]]
    if timings is not None:
        timings['edge_seconds'] += time.monotonic()-started
    return contours, submerged, submerged_rim, contact


def _flat_reference_waterline(hull, hull_faces, lower, upper):
    """Initial posed physical hull intersected with flat world z=0, once."""
    water = np.array([[lower[0], lower[1], 0.], [upper[0], lower[1], 0.],
                      [upper[0], upper[1], 0.], [lower[0], upper[1], 0.]])
    return waterline_segments(hull, hull_faces, water, np.array([[0, 1, 2], [0, 2, 3]]))


def _mask_path(polygons):
    """One standard compound fill preserving exact contact-loop winding."""
    from matplotlib.path import Path as ArtistPath
    groups = {}
    for polygon in polygons:
        groups.setdefault(len(polygon), []).append(polygon)
    paths = [ArtistPath.make_compound_path_from_polys(np.asarray(group))
             for group in groups.values()]
    return ArtistPath.make_compound_path(*paths)


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
    reference = _flat_reference_waterline(positions[0]+local_hull @ rotations[0].T,
                                          result.hull_faces, lower, upper)
    preparation_seconds = time.monotonic()-started

    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import colors, pyplot as plt, tri
    from matplotlib.collections import LineCollection
    from matplotlib.patches import PathPatch
    import imageio_ffmpeg

    figure = plt.figure(figsize=(8, 8), dpi=100)
    axes = figure.add_axes((0, 0, 1, 1), aspect='equal')
    axes.set_xlim(lower[0], upper[0]); axes.set_ylim(lower[1], upper[1]); axes.set_axis_off()
    normalization = colors.TwoSlopeNorm(vmin=-magnitude, vcenter=0., vmax=magnitude)
    interior = PathPatch(_mask_path([]), facecolor='#f3f1e8', edgecolor='none',
                         antialiased=False, alpha=1., zorder=1)
    axes.add_patch(interior)
    reference_artist = LineCollection(reference[:, :, :2], colors='#777777', linewidths=1.,
                                      linestyles='dashed', zorder=1.5, label='reference waterline')
    mesh_artist = LineCollection([], colors='#b0b0b0', linewidths=.7, zorder=2,
                                 label='boat mesh')
    rim_artist = LineCollection([], colors='black', linewidths=1.5, zorder=3,
                                label='boat outline')
    contact_artist = LineCollection([], colors='blue', linewidths=1.3, linestyles='solid',
                                    zorder=4, label='physical waterline')
    for artist in (reference_artist, mesh_artist, rim_artist, contact_artist):
        axes.add_collection(artist)
    header = figure.text(.5, .985, '', ha='center', va='top', fontsize=14.5,
                         fontfamily='DejaVu Sans Mono', color='black',
                         bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .85, 'pad': 3})
    for artist in (interior, reference_artist, mesh_artist, rim_artist, contact_artist, header):
        artist.set_animated(True)
    figure.canvas.draw()
    background = figure.canvas.copy_from_bbox(figure.bbox)
    setup_seconds = time.monotonic()-started-preparation_seconds
    geometry_seconds = draw_seconds = pipe_seconds = 0.
    validation_seconds = encoder_start_seconds = encoder_close_seconds = verification_seconds = 0.
    mask_path_seconds = 0.
    display_timings = dict(mask_seconds=0., edge_seconds=0., contact_seconds=0.)
    completed_frames, phase = 0, 'encoder startup'
    phase_started = time.monotonic()
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
                    phase, phase_started = 'frame budget check', before
                    if before-started > settings.maximum_seconds:
                        raise TimeoutError('render exceeded maximum_seconds')
                    phase = 'geometry validation'
                    phase_started = time.monotonic()
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
                    validation_seconds += time.monotonic()-before
                    phase = 'wetted display geometry'
                    phase_started = time.monotonic()
                    masks, _, _, contact = _wetted_display(
                        hull, result.hull_faces, points, faces, loops, display_timings, projected,
                        compute_submerged_edges=False)
                    path_started = time.monotonic()
                    interior.set_path(_mask_path(masks))
                    mask_path_seconds += time.monotonic()-path_started
                    mesh_artist.set_segments(hull[edges, :2])
                    rim_artist.set_segments(hull[rim, :2])
                    contact_artist.set_segments(contact[:, :, :2])
                    header.set_text(header_text(stamp, angles[index, 0], angles[index, 1],
                                                positions[index, 2]-positions[0, 2]))
                    geometry_seconds += time.monotonic()-before
                    drawing = time.monotonic()
                    phase = 'Agg drawing'
                    phase_started = time.monotonic()
                    figure.canvas.restore_region(background)
                    for artist in (water, interior, reference_artist, mesh_artist, rim_artist, contact_artist):
                        axes.draw_artist(artist)
                    figure.draw_artist(header)
                    figure.canvas.blit(figure.bbox)
                    draw_seconds += time.monotonic()-drawing
                    sending = time.monotonic()
                    phase = 'encoder send'
                    phase_started = time.monotonic()
                    writer.send(np.asarray(figure.canvas.buffer_rgba()))
                    pipe_seconds += time.monotonic()-sending
                    completed_frames += 1
                    if completed_frames == 1 or completed_frames % 50 == 0 or completed_frames == len(stamps):
                        print(f'render progress: {completed_frames}/{len(stamps)} frames; '
                              f'elapsed={time.monotonic()-started:.3f}s '
                              f'validation={validation_seconds:.3f}s '
                              f'mask={display_timings["mask_seconds"]:.3f}s '
                              f'path={mask_path_seconds:.3f}s '
                              f'edges={display_timings["edge_seconds"]:.3f}s '
                              f'contact={display_timings["contact_seconds"]:.3f}s '
                              f'draw={draw_seconds:.3f}s send={pipe_seconds:.3f}s', flush=True)
            finally:
                phase = 'encoder finalization' if completed_frames == len(stamps) else phase
                closing = time.monotonic()
                if completed_frames == len(stamps):
                    phase_started = closing
                try:
                    writer.close()
                finally:
                    encoder_close_seconds = time.monotonic()-closing
            verifying = time.monotonic()
            phase = 'full video verification'
            phase_started = verifying
            try:
                _verify_video(encoded, len(stamps), fps)
            finally:
                verification_seconds = time.monotonic()-verifying
            if time.monotonic()-started > settings.maximum_seconds:
                raise TimeoutError('render exceeded maximum_seconds')
            phase = 'video publication'
            phase_started = time.monotonic()
            _publish_video(encoded, output)
        return dict(complete=True, output=str(output), frames=len(stamps), fps=fps,
                    size=(800, 800), source_end_s=float(stamps[-1]),
                    video_duration_s=len(stamps)/fps, preparation_seconds=preparation_seconds,
                    setup_seconds=setup_seconds, geometry_seconds=geometry_seconds,
                    draw_seconds=draw_seconds, encoder_send_seconds=pipe_seconds,
                    encoder_start_seconds=encoder_start_seconds, encoder_close_seconds=encoder_close_seconds,
                    verification_seconds=verification_seconds, elapsed_s=time.monotonic()-started,
                    topology_changes=topology_changes, ambiguous_euler_labels=ambiguous_euler,
                    validation_seconds=validation_seconds, mask_path_seconds=mask_path_seconds, **display_timings,
                    matplotlib_version=matplotlib.__version__, backend='Matplotlib Agg CPU blit')
    except BaseException:
        print(f'render failed: phase={phase}; completed={completed_frames}/{len(stamps)} frames; '
              f'elapsed={time.monotonic()-started:.3f}s '
              f'active_phase_elapsed={time.monotonic()-phase_started:.3f}s '
              f'preparation={preparation_seconds:.3f}s '
              f'setup={setup_seconds:.3f}s validation={validation_seconds:.3f}s '
              f'mask={display_timings["mask_seconds"]:.3f}s '
              f'path={mask_path_seconds:.3f}s '
              f'edges={display_timings["edge_seconds"]:.3f}s '
              f'contact={display_timings["contact_seconds"]:.3f}s '
              f'draw={draw_seconds:.3f}s send={pipe_seconds:.3f}s '
              f'encoder_start={encoder_start_seconds:.3f}s encoder_close={encoder_close_seconds:.3f}s '
              f'verification={verification_seconds:.3f}s', file=sys.stderr, flush=True)
        raise
    finally:
        plt.close(figure)
