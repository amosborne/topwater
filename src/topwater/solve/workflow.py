"""One preparation, motion, reconstruction and resource-lifetime workflow."""
from contextlib import contextmanager
from dataclasses import replace
import json
from pathlib import Path
import shutil
import sys
import tempfile
import numpy as np
from topwater.workflow_data import Result
from .preparation import prepare_case, waves_from_sea
from .field_backend import select_field_backend
from .sampling import surface_sample_times
from .settings import DEFAULTS, DEFAULT_DWELL_S, VIDEO_SIZE_PX
from .types import ScratchCleanupError
from . import surface as producer


@contextmanager
def scratch_directory(solver):
    """Retain recovery data and preserve the original exception on cleanup failure."""
    directory = Path(tempfile.mkdtemp(prefix=f'topwater-{solver}-'))
    label = 'WEC-Sim' if solver == 'wecsim' else 'Native'
    retain = False
    primary = None
    try:
        yield directory
    except BaseException as error:
        primary = error
        if getattr(error, 'retain_workspace', False) or getattr(error, 'recovery_directory', None) is not None:
            retain = True
            print(f'{label}: temporary recovery directory retained: {directory}', file=sys.stderr)
            error.add_note(f'{label} recovery workspace retained: {directory}')
        raise
    finally:
        if not retain:
            try:
                shutil.rmtree(directory)
            except Exception as error:
                message = f'{label}: temporary cleanup failed; recovery directory retained: {directory}: {error}'
                print(message, file=sys.stderr)
                if primary is None:
                    raise ScratchCleanupError(message) from error
                primary.add_note(message)


def solve(boat, sea, dwell_s=DEFAULT_DWELL_S, *, rest_s=10., ramp_s=5.,
          solver='native', influence=None, device='auto', options=None):
    options = DEFAULTS if options is None else options
    if not isinstance(device, str) or not device.strip():
        raise ValueError('device must be auto, cpu, gpu, or a GPU name')
    if solver == 'native':
        from .native import backend as motion_backend
    elif solver == 'wecsim':
        if device != 'auto':
            raise ValueError('--device requires --solver=native')
        from .wecsim import backend as motion_backend
    else:
        raise ValueError('solver must be native or wecsim')
    runtime = motion_backend.preflight(options, influence)
    setup = prepare_case(boat.to_body(), waves_from_sea(sea), dwell_s,
        rest_s=rest_s, ramp_s=ramp_s, time_step_s=options.time_step_s,
        envelope_tolerance_m=options.packet_tolerance_m)
    setup.train.certify_incident_startup(setup)
    with scratch_directory(solver) as workspace:
        field_backend = None
        try:
            prepared = motion_backend.prepare(setup, workspace, options, runtime)
            hydro = prepared.hydrodynamics
            motion = motion_backend.integrate(setup, prepared, workspace, options, influence, runtime)
            field_backend = select_field_backend(hydro, 'cpu' if solver == 'wecsim' else device,
                cpu_workers=options.cpu_workers)
            surface_time = surface_sample_times(motion.time_s, every=options.render_every)
            offset = np.array([0., 0., setup.world_vertical_offset_m])
            computational = replace(motion, position_m=motion.position_m-offset)
            marker = workspace/'surface-world-offset.json'
            marker.write_text(json.dumps(dict(world_vertical_offset_m=float(offset[2]),
                stored_surface_frame='computational_mean_zero',
                recovery='Translate recovered hull/body/surface world z by this offset once.'))+'\n', encoding='utf-8')
            try:
                surface = producer.reconstruct_total_surface(hydro, setup, computational,
                    surface_time, field_backend=field_backend, surface_options=options)
                vertices, faces, vo, fo = surface
            except BaseException as error:
                recovery = getattr(error, 'recovery_directory', None)
                if recovery is not None:
                    try:
                        shutil.copy2(marker, Path(recovery)/marker.name)
                    except OSError as copy_error:
                        error.add_note(f'World-offset marker remains at {marker}: {copy_error}')
                raise
            vertices = np.asarray(vertices).copy()
            vertices[:, 2] += offset[2]
            return Result(time_s=motion.time_s, position_m=motion.position_m,
                orientation_wxyz=motion.orientation_wxyz, hull_vertices_m=boat.vertices_m,
                hull_faces=boat.faces, cg_m=boat.cg_m, surface_time_s=surface_time,
                surface_vertices_m=vertices, surface_faces=faces,
                surface_vertex_offsets=vo, surface_face_offsets=fo,
                surface_interval_s=options.surface_interval_s,
                view_bounds_xy_m=getattr(surface, 'view_bounds_xy_m', None), video_size_px=VIDEO_SIZE_PX)
        finally:
            if field_backend is not None:
                primary = sys.exception()
                try:
                    field_backend.close()
                except Exception as error:
                    if primary is None:
                        raise
                    primary.add_note(f'Field cleanup also failed: {error}')
