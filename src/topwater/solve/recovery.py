"""Completed-field recovery evidence and solver-independent replay."""
from pathlib import Path
import sys
import traceback
import numpy as np
from time import perf_counter
from .validation import _finite, _grid
from .evidence import SurfaceExportError, _save_frame_failure


def replay_surface_cache(cache_path, *, diagnostic_directory, output_path=None, analyze=False):
    """Pack frozen completed fields and motion; optional historical analysis only.

    No Capytaine, MATLAB, body integration or wave synthesis is invoked. The
    optional output is the ordinary final Result format. Keep the cache and
    geometry source identities fixed when using the replay.
    """
    path = Path(cache_path)
    with np.load(path, allow_pickle=False) as archive:
        data = {name: archive[name].copy() for name in archive.files}
    if int(data['cache_version'].item()) != 1:
        raise ValueError('Unsupported private completed-surface cache')
    directory = Path(diagnostic_directory)
    directory.mkdir(parents=True, exist_ok=False)
    print(f'Surface geometry replay: frozen field source={str(data["producer_source_sha256"].item())}; '
          f'cache={path.resolve()}', file=sys.stderr)
    packed = _export_completed_field(data, directory, path, analyze=analyze)
    if output_path is not None:
        from topwater.workflow_data import Result
        # The internal closed hull can contain cap-only vertices. Public Result
        # stores physical panels only, so compact that snapshot without changing
        # the closed geometry/rim/cap identities in the recovery cache.
        physical = np.asarray(data['physical_faces'], dtype=np.int64)
        used, remapped = np.unique(physical, return_inverse=True)
        physical_vertices = data['body_vertices_m'][used]
        physical_faces = remapped.reshape(physical.shape)
        world_shift = np.array([0., 0., float(data.get('world_vertical_offset_m', 0.))])
        result = Result(time_s=data['motion_time_s'], position_m=data['motion_position_m']+world_shift,
            orientation_wxyz=data['motion_orientation_wxyz'], hull_vertices_m=physical_vertices,
            hull_faces=physical_faces, cg_m=data['body_cg_m'], surface_time_s=data['surface_time_s'],
            surface_vertices_m=packed[0]+world_shift, surface_faces=packed[1],
            surface_vertex_offsets=packed[2], surface_face_offsets=packed[3],
            surface_interval_s=data.get('surface_interval_s'), view_bounds_xy_m=data.get('view_bounds_xy_m'),
            video_size_px=data.get('video_size_px'))
        result.save(output_path, overwrite=False)
    return packed


def _export_completed_field(completed, recovery_directory, cache_path, *, timings=None, analyze=False):
    if analyze:
        from .analysis import _analyze_completed_field_impl
    try:
        return (_analyze_completed_field_impl if analyze else _export_completed_field_impl)(
            completed, recovery_directory, cache_path, timings=timings)
    except SurfaceExportError:
        raise
    except (KeyboardInterrupt, SystemExit) as error:
        error.recovery_directory = Path(recovery_directory)
        print(f'Surface export interrupted; completed fields retained at {cache_path}', file=sys.stderr)
        raise
    except Exception as error:
        stage = 'completed-field-validation-or-summary'
        # No frame/time is established if completed-data validation itself
        # failed. Reading the malformed axis here would mask that first error.
        stamp = None
        print(f'Surface {stage} failed: {type(error).__name__}: {error}', file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        _save_frame_failure(recovery_directory, {}, error, stage, -1, stamp, cache_path)
        raise SurfaceExportError(f'Surface {stage} failed: {error}; completed fields retained at {cache_path}',
            recovery_directory, stage=stage, frame_index=-1, time_s=stamp) from error


def _export_completed_field_impl(completed, recovery_directory, cache_path, *, timings=None):
    """Validate and pack the unchanged grid; moving contact belongs to render."""
    started = perf_counter()
    stamp = _grid(completed['surface_time_s'], 'cached surface_time_s', zero=True)
    _finite(completed['surface_position_m'], (len(stamp), 3), 'cached surface positions')
    _finite(completed['surface_rotation'], (len(stamp), 3, 3), 'cached surface rotations')
    xy = _finite(completed['xy_native_m'], (None, 2), 'cached XY')
    total = _finite(completed['total_elevation_m'], (len(stamp), len(xy)), 'cached total field')
    _finite(completed['incident_elevation_m'], total.shape, 'cached incident field')
    embedding = _finite(completed['world_from_native'], (3, 3), 'cached embedding')
    raw = np.asarray(completed['water_faces'])
    if raw.ndim != 2 or raw.shape[1] != 3 or raw.dtype.kind not in 'iu':
        raise ValueError('Water graph connectivity must be integer triangles')
    if not len(raw) or raw.min() < 0 or raw.max() >= len(xy):
        raise ValueError('Water graph connectivity is empty or invalid')
    faces = raw.astype(np.int64, copy=False)
    triangle = xy[faces]
    a, b = triangle[:, 1]-triangle[:, 0], triangle[:, 2]-triangle[:, 0]
    if np.any(a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0] <= 0):
        raise ValueError('Water graph triangles must have positive XY area and upward winding')
    # Keep the exact per-frame transform/order used by the former assembly.
    vertices = np.concatenate([np.c_[xy, height] @ embedding.T for height in total])
    packed = (vertices, np.tile(faces, (len(stamp), 1)),
              np.arange(len(stamp)+1, dtype=np.int64)*len(xy),
              np.arange(len(stamp)+1, dtype=np.int64)*len(faces))
    elapsed = perf_counter()-started
    if timings is not None:
        timings['total-postprocess'] = elapsed
    print(f'Surface packing wall time_s: {elapsed:.6g}', file=sys.stderr)
    return packed
