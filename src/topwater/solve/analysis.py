"""Optional completed-field diagnostics; never part of normal reconstruction."""
import json
import sys
import traceback
from time import perf_counter
import numpy as np
from scipy.integrate import trapezoid
from scipy.spatial import ConvexHull
from .validation import _finite, _grid
from .progress import _StageProgress
from .evidence import SurfaceExportError, _save_frame_failure


def _clip_area_polygon(points, normal, offset, tolerance):
    """Convex planar polygon intersection with normal.point+offset<=0."""
    if len(points) < 3:
        return np.empty((0, 3))
    # Keep the former individual dot products (and rounding), evaluating each
    # once. Most prism/plane clips do not actually cross the polygon.
    distances = [float(point @ normal+offset) for point in points]
    if min(distances) > tolerance:
        return np.empty((0, 3))
    if max(distances) <= tolerance:
        edges = points-np.concatenate((points[-1:], points[:-1]))
        # This sufficient condition implies every original norm-based cleanup
        # test keeps its vertex. Borderline/small-edge cases use the unchanged
        # cleanup below; no tolerance or indistinguishable-point rule changes.
        if np.all(np.max(abs(edges), axis=1) > tolerance):
            return points
    result = []
    previous = points[-1]
    dp = distances[-1]
    for point, dc in zip(points, distances):
        before, after = dp <= tolerance, dc <= tolerance
        if before != after:
            denominator = dp-dc
            if abs(denominator) > np.finfo(float).tiny:
                result.append(previous+(dp/denominator)*(point-previous))
        if after:
            result.append(point)
        previous, dp = point, dc
    if len(result) < 3:
        return np.empty((0, 3))
    clean = [result[0]]
    for point in result[1:]:
        if np.linalg.norm(point-clean[-1]) > tolerance:
            clean.append(point)
    if len(clean) > 1 and np.linalg.norm(clean[0]-clean[-1]) <= tolerance:
        clean.pop()
    return np.asarray(clean) if len(clean) >= 3 else np.empty((0, 3))


def _polygon_area(points):
    if len(points) < 3:
        return 0.
    relative = points-points[0]
    following = np.concatenate((relative[1:], relative[:1]))
    # Same products/subtractions, row order and C(n,3) sum as np.cross/roll,
    # avoiding generic axis dispatch for these tiny integration polygons.
    cross = np.column_stack((relative[:, 1]*following[:, 2]-relative[:, 2]*following[:, 1],
                             relative[:, 2]*following[:, 0]-relative[:, 0]*following[:, 2],
                             relative[:, 0]*following[:, 1]-relative[:, 1]*following[:, 0]))
    return float(np.linalg.norm(cross.sum(axis=0))/2)


class WaterGraphPartition:
    """Fixed XY triangle prisms for exact represented-P1 hull wet-set area.

    The triangulation is a single-valued upward graph. Shared prism edges use
    deterministic ownership when a vertical hull facet lies exactly on one,
    preventing double area counting without moving the physical hull.
    """
    def __init__(self, xy, faces):
        self.xy = _finite(xy, (None, 2), 'water xy')
        raw = np.asarray(faces)
        if raw.ndim != 2 or raw.shape[1] != 3 or raw.dtype.kind not in 'iu':
            raise ValueError('Water graph connectivity must be integer triangles')
        if not len(raw) or raw.min() < 0 or raw.max() >= len(xy):
            raise ValueError('Water graph connectivity is empty or invalid')
        self.faces = raw.astype(np.int64, copy=False)
        triangle = self.xy[self.faces]
        edge = np.roll(triangle, -1, axis=1)-triangle
        cross = edge[:, 0, 0]*edge[:, 1, 1]-edge[:, 0, 1]*edge[:, 1, 0]
        if np.any(cross <= 0):
            raise ValueError('Water graph triangles must have positive XY area and upward winding')
        normals = np.stack([edge[:, :, 1], -edge[:, :, 0]], axis=-1)
        normals /= np.linalg.norm(normals, axis=-1)[:, :, None]
        self.normals = normals
        self.offsets = -np.sum(normals*triangle, axis=-1)
        self.lower, self.upper = triangle.min(axis=1), triangle.max(axis=1)
        edges = np.stack((self.faces, self.faces[:, [1, 2, 0]]), axis=2).reshape(-1, 2)
        unique, inverse = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True)
        owners = np.full(len(unique), len(self.faces), dtype=np.int64)
        np.minimum.at(owners, inverse, np.repeat(np.arange(len(self.faces)), 3))
        self.owner = owners[inverse].reshape(-1, 3)
        self.tolerance = 1e-10*max(1., float(np.ptp(self.xy, axis=0).max()))

    def wet_set_area(self, current_hull, physical_faces, height_m, reference_normal,
                     reference_offset, initial_wetted_area):
        """Material-associated wet-set symmetric difference / initial wet area.

        ``reference_normal.point+reference_offset<=0`` marks the material that
        was wet at the initial mean pose, expressed in the current posed axes.
        All areas are actual 3D physical-hull area, not projected area or point
        classification. The only surface approximation is the supplied P1 graph.
        """
        hull = _finite(current_hull, (None, 3), 'current_hull')
        height = _finite(height_m, (len(self.xy),), 'height_m')
        normal = _finite(reference_normal, (3,), 'reference_normal')
        reference_offset = float(reference_offset)
        initial_area = float(initial_wetted_area)
        if not np.isfinite(initial_area) or initial_area <= 0:
            raise ValueError('Initial physical wetted area must be positive')
        water = np.c_[self.xy, height][self.faces]
        water_normal = np.cross(water[:, 1]-water[:, 0], water[:, 2]-water[:, 0])
        water_normal /= np.linalg.norm(water_normal, axis=1)[:, None]
        water_offset = -np.sum(water_normal*water[:, 0], axis=1)
        current_area = mean_area = common_area = 0.
        eps = self.tolerance
        for triangle in hull[np.asarray(physical_faces, dtype=np.int64)]:
            lo, hi = triangle[:, :2].min(axis=0), triangle[:, :2].max(axis=0)
            candidate = np.flatnonzero(np.all(self.upper >= lo-eps, axis=1)
                                       & np.all(self.lower <= hi+eps, axis=1))
            for i in candidate:
                polygon = triangle.copy()
                for j in range(3):
                    n = np.r_[self.normals[i, j], 0.]
                    d = self.offsets[i, j]
                    if len(polygon) and np.all(abs(polygon @ n+d) <= eps) and self.owner[i, j] != i:
                        polygon = np.empty((0, 3)); break
                    polygon = _clip_area_polygon(polygon, n, d, eps)
                    if len(polygon) < 3:
                        break
                polygon_area = _polygon_area(polygon)
                if len(polygon) < 3 or polygon_area <= eps**2:
                    continue
                wet = _clip_area_polygon(polygon, water_normal[i], water_offset[i], eps)
                mean = _clip_area_polygon(polygon, normal, reference_offset, eps)
                common = _clip_area_polygon(wet, normal, reference_offset, eps)
                wet_area = polygon_area if wet is polygon else _polygon_area(wet)
                mean_piece_area = polygon_area if mean is polygon else _polygon_area(mean)
                common_piece_area = wet_area if common is wet else _polygon_area(common)
                current_area += wet_area
                mean_area += mean_piece_area
                common_area += common_piece_area
        if abs(mean_area-initial_area) > 1e-7*initial_area:
            raise ValueError('Wet-set area partition does not reproduce initial wetted area; graph coverage or prism ownership is unresolved')
        difference = max(0., current_area+mean_area-2*common_area)
        return difference/initial_area, current_area


def _grid_p1_values(axes, values, points):
    """Values on the exact fixed-diagonal P1 graph, including new cut vertices."""
    x, y = axes
    query = _finite(points, (None, 2), 'P1 points')
    ix = np.clip(np.searchsorted(x, query[:, 0], side='right')-1, 0, len(x)-2)
    iy = np.clip(np.searchsorted(y, query[:, 1], side='right')-1, 0, len(y)-2)
    u = (query[:, 0]-x[ix])/(x[ix+1]-x[ix])
    v = (query[:, 1]-y[iy])/(y[iy+1]-y[iy])
    if np.any(u < -1e-9) or np.any(u > 1+1e-9) or np.any(v < -1e-9) or np.any(v > 1+1e-9):
        raise ValueError('P1 query outside exported water domain')
    z = np.asarray(values).reshape(len(y), len(x))
    lower = v <= u
    answer = np.empty(len(query))
    answer[lower] = ((1-u[lower])*z[iy[lower], ix[lower]]
                     +(u[lower]-v[lower])*z[iy[lower], ix[lower]+1]
                     +v[lower]*z[iy[lower]+1, ix[lower]+1])
    upper = ~lower
    answer[upper] = ((1-v[upper])*z[iy[upper], ix[upper]]
                     +u[upper]*z[iy[upper]+1, ix[upper]+1]
                     +(v[upper]-u[upper])*z[iy[upper]+1, ix[upper]])
    return answer


def _waterline_distance_xy(points, segments):
    if not len(segments):
        raise ValueError('Total surface has no resolved physical waterline; near-band metric is undefined')
    a, b = segments[:, 0, :2], segments[:, 1, :2]
    edge = b-a
    length_squared = np.sum(edge*edge, axis=1)
    # Vertical waterline segments have a point projection; their horizontal
    # distance is still defined, without division by zero or discarding them.
    result = np.empty(len(points))
    for start in range(0, len(points), 1024):
        query = points[start:start+1024]
        relative = query[:, None]-a[None]
        fraction = np.clip(np.sum(relative*edge[None], axis=2)/np.where(length_squared > 0, length_squared, 1.), 0, 1)
        result[start:start+len(query)] = np.linalg.norm(
            query[:, None]-(a[None]+fraction[:, :, None]*edge[None]), axis=2).min(axis=1)
    return result


def _waterline_band_xy(points, segments, radius):
    """Same exact distance mask, with a conservative global box broad phase."""
    if not len(segments):
        raise ValueError('Total surface has no resolved physical waterline; near-band metric is undefined')
    # Candidate expansion only: allow for floating closest-point arithmetic
    # at translated/oblique coordinates. The final distance threshold below
    # is unchanged; this guard cannot admit a point to the reported band.
    guard = 64*np.finfo(float).eps*max(1., float(np.max(abs(segments[:, :, :2]))), radius)
    lower = np.nextafter(segments[:, :, :2].min(axis=(0, 1))-radius-guard, -np.inf)
    upper = np.nextafter(segments[:, :, :2].max(axis=(0, 1))+radius+guard, np.inf)
    candidate = np.all((points >= lower) & (points <= upper), axis=1)
    selected = np.zeros(len(points), dtype=bool)
    selected[candidate] = _waterline_distance_xy(points[candidate], segments) <= radius
    return selected


def _split_metric_polygon(polygon, normal, offset):
    """Temporary P1 integration rows (x,y,total,incident), never an output mesh.

    Interpolation of the complete row carries the SAME source-P1 incident
    value to a cut point. No analytic incident re-evaluation or vertex welding
    is performed. A plane boundary has zero projected measure.
    """
    inside, outside = [], []
    has_positive = False
    previous = polygon[-1]
    dp = float(previous[:3] @ normal+offset)
    for current in polygon:
        dc = float(current[:3] @ normal+offset)
        has_positive |= dc > 0
        if (dp < 0 < dc) or (dc < 0 < dp):
            cut = previous+(dp/(dp-dc))*(current-previous)
            inside.append(cut); outside.append(cut)
        if dc <= 0:
            inside.append(current)
        if dc >= 0:
            outside.append(current)
        previous, dp = current, dc
    return (np.asarray(inside).reshape(-1, 4),
            np.asarray(outside).reshape(-1, 4) if has_positive else np.empty((0, 4)))


def _metric_triangle_rows(triangle):
    """Exact projected P1 moments on algebraic integration triangles only."""
    a, b = triangle[:, 1, :2]-triangle[:, 0, :2], triangle[:, 2, :2]-triangle[:, 0, :2]
    area = abs(a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0])/2
    q = triangle[:, :, 2]-triangle[:, :, 3]
    inc = triangle[:, :, 3]
    first = area*q.mean(axis=1)
    square = area/6*(np.sum(q*q, axis=1)+q[:, 0]*q[:, 1]+q[:, 1]*q[:, 2]+q[:, 2]*q[:, 0])
    inc_square = area/6*(np.sum(inc*inc, axis=1)+inc[:, 0]*inc[:, 1]+inc[:, 1]*inc[:, 2]+inc[:, 2]*inc[:, 0])
    peak_index = np.argmax(abs(q), axis=1)
    maximum = abs(q[np.arange(len(q)), peak_index])
    peak = triangle[np.arange(len(q)), peak_index, :3]
    return np.column_stack((area, first, square, maximum, inc_square, peak))


def _delta_moments(vertices, faces, incident_at_vertices, hull, beam, segments, *, startup_bounds=None,
                   startup_vertical_offset_m=0.):
    """Exact exterior-support P1 moments, without constructing an output mesh.

    A convex solid's exterior is NOT convex. For each boundary source triangle,
    collect the outside of the next plane and carry only the remaining inside
    polygon forward. The collected pieces are disjoint except for boundaries.
    Bulk exterior classification requires a COMMON separating plane, never
    merely three vertices individually outside different hull planes.

    The near band is a documented discrete approximation: select ORIGINAL
    source triangles by centroid distance to the independent total waterline,
    then integrate only their exterior portions. It is not an exact distance
    region, nor the former remeshed-triangle-centroid convention.
    """
    triangle = vertices[faces]
    rows = np.concatenate((triangle, np.asarray(incident_at_vertices)[faces][..., None]), axis=2)
    planes = ConvexHull(hull).equations
    distances = np.einsum('tvi,pi->tvp', triangle, planes[:, :3])+planes[:, 3]
    exterior = np.any(np.all(distances > 0, axis=1), axis=1)
    interior = np.all(distances <= 0, axis=(1, 2))
    near = _waterline_band_xy(triangle[:, :, :2].mean(axis=1), segments, beam)
    output = {name: np.zeros(8) for name in ('full', 'one_beam')}
    regional = np.zeros(len(faces), dtype=bool)
    if startup_bounds is not None:
        lo, hi = np.asarray(startup_bounds)
        regional = (np.all(triangle[:, :, :2].max(axis=1) >= lo, axis=1)
                    & np.all(triangle[:, :, :2].min(axis=1) <= hi, axis=1))
    startup_area, startup_peak = 0., 0.

    def accumulate(parts, near_mask, regional_mask):
        nonlocal startup_area, startup_peak
        stats = _metric_triangle_rows(parts)
        positive = stats[:, 0] > 0
        for name, mask in [('full', positive), ('one_beam', positive & near_mask)]:
            selected = stats[mask]
            if not len(selected):
                continue
            was_empty = output[name][0] == 0
            output[name][[0, 1, 2, 4]] += selected[:, [0, 1, 2, 4]].sum(axis=0)
            winner = selected[np.argmax(selected[:, 3])]
            if winner[3] > output[name][3] or was_empty:
                output[name][3] = winner[3]
                output[name][5:8] = winner[5:8]
        local = positive & regional_mask
        if np.any(local):
            startup_area += float(stats[local, 0].sum())
            startup_peak = max(startup_peak, float(np.max(abs(parts[local, :, 2]+startup_vertical_offset_m))))

    if np.any(exterior):
        accumulate(rows[exterior], near[exterior], regional[exterior])
    for index in np.flatnonzero(~(exterior | interior)):
        remaining = rows[index]
        for plane in planes:
            if len(remaining) < 3:
                break
            remaining, outside = _split_metric_polygon(remaining, plane[:3], plane[3])
            if len(outside) < 3:
                continue
            # Fan moments are algebra only: no connectivity is registered,
            # published, validated as a mesh, or fed back to the field. Zero
            # area fan pieces contribute exactly zero and need no ear repair.
            fan = np.stack((np.repeat(outside[:1], len(outside)-2, axis=0),
                            outside[1:-1], outside[2:]), axis=1)
            accumulate(fan, near[index], regional[index])
    if startup_bounds is not None and startup_area <= 0:
        raise ValueError('Initial total surface lacks positive exterior boat-region support')
    return output, (startup_peak if startup_bounds is not None else None)


def _analyze_completed_field_impl(completed, recovery_directory, cache_path, *, timings=None):
    """Geometry and metrics only; all physical field samples are immutable."""
    from topwater.geometry.waterline import validate_opening_clearance, waterline_segments

    started = perf_counter()
    stage_seconds = {} if timings is None else timings

    def measured(name, function, *args, **kwargs):
        start = perf_counter()
        try:
            return function(*args, **kwargs)
        finally:
            stage_seconds[name] = stage_seconds.get(name, 0.)+perf_counter()-start

    surface_time = _grid(completed['surface_time_s'], 'cached surface_time_s', zero=True)
    positions = _finite(completed['surface_position_m'], (len(surface_time), 3), 'cached surface positions')
    poses = _finite(completed['surface_rotation'], (len(surface_time), 3, 3), 'cached surface rotations')
    full_vertices = _finite(completed['body_vertices_m'], (None, 3), 'cached body vertices')
    local_hull = full_vertices-_finite(completed['body_cg_m'], (3,), 'cached body CG')
    full_faces = np.asarray(completed['body_faces'], dtype=np.int64)
    rim = np.asarray(completed['opening_rim'], dtype=np.int64)
    caps = np.asarray(completed['artificial_faces'], dtype=np.int64)
    physical_faces = np.asarray(completed['physical_faces'], dtype=np.int64)
    xy = _finite(completed['xy_native_m'], (None, 2), 'cached XY')
    faces = np.asarray(completed['water_faces'], dtype=np.int64)
    axes = [np.asarray(completed['axis_x_m']), np.asarray(completed['axis_y_m'])]
    total = _finite(completed['total_elevation_m'], (len(surface_time), len(xy)), 'cached total field')
    incident = _finite(completed['incident_elevation_m'], total.shape, 'cached incident field')
    embedding = _finite(completed['world_from_native'], (3, 3), 'cached embedding')
    initial_area = float(completed['initial_wetted_area_m2'])
    beam = float(completed['beam_m'])
    calm_limit = float(completed['calm_eta_limit_m'])
    world_vertical_offset = float(completed.get('world_vertical_offset_m', 0.))
    packed_vertices, packed_faces = [], []
    vertex_offsets, face_offsets = [0], [0]
    moments = {name: [] for name in ('full', 'one_beam')}
    local_initial_peak = None
    progress = _StageProgress('surface geometry and diagnostics', len(surface_time), 'frames')
    for i, t in enumerate(surface_time):
        world_hull = local_hull @ poses[i].T+positions[i]
        native_hull = world_hull @ embedding
        native_points = np.c_[xy, total[i]]
        world_points = native_points @ embedding.T
        evidence = dict(frame_index=np.array(i), time_s=np.array(t),
            water_vertices_world_m=world_points, water_vertices_native_m=native_points,
            water_faces=faces, hull_vertices_world_m=world_hull, hull_vertices_native_m=native_hull,
            hull_faces=full_faces, physical_faces=physical_faces, opening_rim=rim, artificial_faces=caps,
            incident_elevation_m=incident[i], axis_x_m=axes[0], axis_y_m=axes[1],
            initial_wetted_area_m2=np.array(initial_area))
        stage = 'opening-clearance'
        try:
            measured(stage, validate_opening_clearance, world_points, faces, world_hull, full_faces,
                     opening_rim=rim, artificial_faces=caps)
            stage = 'independent-waterline'
            waterline_world = measured(stage, waterline_segments, world_hull, full_faces,
                                      world_points, faces, artificial_faces=caps)
            evidence['waterline_segments_world_m'] = waterline_world
            stage = 'waterline-delta-metrics'
            # Only integration bookkeeping uses temporary exterior polygons.
            # The original samples/connectivity remain the published geometry.
            startup_bounds = ([native_hull[:, :2].min(axis=0), native_hull[:, :2].max(axis=0)]
                              if i == 0 else None)
            frame_moments, initial_peak = measured(stage, _delta_moments, native_points, faces,
                incident[i], native_hull, beam, waterline_world @ embedding, startup_bounds=startup_bounds,
                startup_vertical_offset_m=world_vertical_offset)
            for name in moments:
                moments[name].append(frame_moments[name])
            stage = 'startup-total-elevation'
            if i == 0:
                local_initial_peak = initial_peak
                if local_initial_peak > calm_limit+1e-12:
                    raise ValueError(f'Initial represented total surface exceeds local calm eta limit: {local_initial_peak:g} > {calm_limit:g} m')
            stage = 'pack-frame'
            packed_vertices.append(world_points); packed_faces.append(faces)
            vertex_offsets.append(vertex_offsets[-1]+len(world_points))
            face_offsets.append(face_offsets[-1]+len(faces))
            progress.update(i+1)
        except Exception as error:
            print(f'Surface {stage} failed at frame {i}, time {t:g}s: {type(error).__name__}: {error}', file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            _save_frame_failure(recovery_directory, evidence, error, stage, i, t, cache_path)
            raise SurfaceExportError(f'Surface {stage} failed at frame {i}, time {t:g}s: {error}; '
                f'completed fields and failure geometry retained at {recovery_directory}', recovery_directory,
                stage=stage, frame_index=i, time_s=t) from error
    summary = {}
    for name, frames in moments.items():
        value = np.asarray(frames)
        area_time = float(trapezoid(value[:, 0], x=surface_time))
        if area_time <= 0:
            raise ValueError(f'Surface delta metric has no positive {name} area-time support')
        incident_rms = float(np.sqrt(max(0., trapezoid(value[:, 4], x=surface_time)/area_time)))
        rms = float(np.sqrt(max(0., trapezoid(value[:, 2], x=surface_time)/area_time)))
        resolved_peak = int(np.argmax(value[:, 3]))
        denominator_floor = 64*np.finfo(float).eps*max(1., float(np.max(abs(incident))))
        summary[name] = dict(area_time_m2_s=area_time,
            mean_delta_m=float(trapezoid(value[:, 1], x=surface_time)/area_time),
            rms_delta_m=rms, incident_rms_m=incident_rms,
            relative_rms=(rms/incident_rms if incident_rms > denominator_floor else None),
            relative_rms_status=('defined' if incident_rms > denominator_floor else 'degenerate incident RMS'),
            maximum_abs_delta_m=float(value[resolved_peak, 3]),
            resolved_peak_time_s=float(surface_time[resolved_peak]),
            resolved_peak_world_xyz_m=(value[resolved_peak, 5:8] @ embedding.T).tolist())
    print('Surface total-minus-incident metrics (projected area/time; one-beam band uses '
          'ORIGINAL source-triangle-centroid horizontal distance to independent total-surface/physical-hull '
          'waterline; exact exterior portions only, not an exact distance-band region): '+json.dumps(summary), file=sys.stderr)
    print(f'Surface startup: represented total local |eta|<={local_initial_peak:g}m on positive-area '
          'exterior portions of source faces whose XY boxes overlap the initial hull XY bounds '
          '(conservative support may extend beyond that rectangle). '
          'The separate full-hull incident eta/velocity bounds remain in force; '
          'no total disturbed-fluid velocity bound has been established.', file=sys.stderr)
    packed = (np.concatenate(packed_vertices), np.concatenate(packed_faces),
              np.asarray(vertex_offsets, dtype=np.int64), np.asarray(face_offsets, dtype=np.int64))
    stage_seconds['total-postprocess'] = perf_counter()-started
    print('Surface postprocess wall times_s: '+json.dumps(stage_seconds), file=sys.stderr)
    progress.finish()
    return packed
