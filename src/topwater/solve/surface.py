"""Shared water-field reconstruction from prepared hydrodynamics and saved motion.

Harmonics use Re(A exp(-i omega t)) in fixed mean-world axes about the mean CG.
Current reconstruction preserves the complete native velocity FFT. Explicit
completed-field analysis remains separate from publication and integration.
"""
from dataclasses import dataclass
from .fourier import surface_radiation_frequencies, offline_radiation_response
from .types import MODES
from .validation import _finite, _grid
from .progress import _StageProgress, _summarize_capytaine_warnings
from topwater.solve.settings import DEFAULTS, VIDEO_SIZE_PX, SURFACE_PIXELS_PER_SAMPLE
from functools import wraps
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import ConvexHull, QhullError


@dataclass(frozen=True)
class SurfaceOutput:
    vertices_m: np.ndarray
    faces: np.ndarray
    vertex_offsets: np.ndarray
    face_offsets: np.ndarray
    view_bounds_xy_m: np.ndarray

    def __iter__(self):
        return iter((self.vertices_m, self.faces, self.vertex_offsets, self.face_offsets))

    def __getitem__(self, index):
        return tuple(self)[index]

    def __len__(self):
        return 4


def surface_spacing_m(view_bounds_xy_m):
    """Maximum world spacing for ten pixels between samples in the saved view."""
    lower, upper = np.asarray(view_bounds_xy_m)
    return float(np.min((upper-lower)/VIDEO_SIZE_PX)*SURFACE_PIXELS_PER_SAMPLE)


from .evidence import _save_completed_cache
from .recovery import _export_completed_field


def _single_threaded(function):
    """Set the existing one-thread policy once, not once per scalar BEM solve."""
    @wraps(function)
    def call(*args, **kwargs):
        from threadpoolctl import threadpool_limits
        with threadpool_limits(1):
            return function(*args, **kwargs)
    return call


def exterior_field_points(hydro, xy, *, allow_boundary=False):
    """Admit only points strictly outside the mean-waterline footprint.

    Capytaine's linear field is exterior to the fixed mean hull.  Moving-hull
    clipping cannot make an interior field evaluation physical.  The waterline
    convex hull is a conservative footprint: concave waterlines may reject
    extra points, but an interior point is never admitted by that approximation.
    Explicit boundary-trace mode admits the mean-waterline boundary within
    roundoff tolerance, but still rejects strictly interior Green queries.
    """
    points = _finite(xy, (None, 2), 'xy')
    vertices = _finite(hydro.body.mesh.vertices, (None, 3), 'mean hull vertices')
    scale = max(1., float(np.max(abs(vertices))), float(np.ptp(vertices, axis=0).max()))
    tolerance = 64*np.finfo(float).eps*scale
    waterline = vertices[abs(vertices[:, 2]) <= tolerance, :2]
    if len(waterline) < 3:
        raise ValueError('Exterior field evaluation requires a resolved mean waterline')
    try:
        footprint = ConvexHull(waterline)
    except QhullError as error:
        raise ValueError('Mean waterline must enclose a two-dimensional footprint') from error
    distances = points @ footprint.equations[:, :2].T + footprint.equations[:, 2]
    if allow_boundary:
        if np.any(np.all(distances < -tolerance, axis=1)):
            raise ValueError('Strictly interior mean-hull Green-function evaluation is prohibited')
    elif np.any(np.all(distances <= tolerance, axis=1)):
        raise ValueError('Field points must lie strictly outside the mean hull footprint; interior/on-waterline Green-function evaluation is prohibited')
    return points


def exterior_trace_points(hydro, xy, offset_m, *, allow_boundary=False):
    """Map a mean-footprint interior strip to an explicit exterior trace.

    Only disturbance fields may use these mapped points. The incident field
    must still be evaluated at the original XY. This is a bounded geometric
    continuation for export, not an exact nonlinear contact condition.
    """
    points = _finite(xy, (None, 2), 'xy')
    offset = float(offset_m)
    if not np.isfinite(offset) or offset < 0 or (offset == 0 and not allow_boundary):
        raise ValueError('Exterior trace offset must be positive and finite unless explicit boundary mode permits zero')
    vertices = _finite(hydro.body.mesh.vertices, (None, 3), 'mean hull vertices')
    scale = max(1., float(np.max(abs(vertices))))
    rim = vertices[abs(vertices[:, 2]) <= 64*np.finfo(float).eps*scale, :2]
    if len(rim) < 3:
        raise ValueError('Exterior trace requires a resolved mean waterline')
    footprint = ConvexHull(rim)
    polygon = rim[footprint.vertices]
    a, b = polygon, np.roll(polygon, -1, axis=0)
    edge = b-a
    normal = np.c_[edge[:, 1], -edge[:, 0]]
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    mapped = points.copy()
    continued = np.zeros(len(points), dtype=bool)
    for i, point in enumerate(points):
        fraction = np.clip(np.sum((point-a)*edge, axis=1)/np.sum(edge*edge, axis=1), 0, 1)
        projected = a+fraction[:, None]*edge
        distance = np.linalg.norm(projected-point, axis=1)
        j = int(np.argmin(distance))
        inside = np.all(point @ footprint.equations[:, :2].T + footprint.equations[:, 2] <= 0)
        if inside or distance[j] < offset:
            mapped[i] = projected[j]+offset*normal[j]
            continued[i] = True
    exterior_field_points(hydro, mapped, allow_boundary=allow_boundary)
    return mapped, continued


def incident_spectrum(train, xy):
    """Exact already-weighted train incident spectrum at physical target XY."""
    points = _finite(xy, (None, 2), 'xy')
    omega = _grid(train.omega, 'train.omega', zero=True)
    C = _finite(train.C, (None, 2, len(omega)), 'train.C', complex_value=True)
    theta = _finite(train.theta, (len(C),), 'train.theta')
    k = _finite(train.k, (len(omega),), 'train.k')
    output = np.zeros((len(omega), len(points)), dtype=complex)
    for j in range(len(C)):
        s = points @ np.array([np.cos(theta[j]), np.sin(theta[j])])
        output += C[j, 0, :, None]*np.exp(1j*k[:, None]*s)
        output += C[j, 1, :, None]*np.exp(-1j*k[:, None]*s)
    return output


@_summarize_capytaine_warnings
@_single_threaded
def combined_diffraction(hydro, train, xy, *, return_forces=False, checkpoint=None, return_incident=True, allow_boundary=False, field_backend):
    """Exact sideband combined incident/scattered spectra on the frozen hull.

    ``train.C`` already includes spectral quadrature weights.  Its two sign
    branches are two positive-frequency problems with opposite propagation
    headings; no negative-frequency outgoing Green function is evaluated.
    Returns complex (omega, point) elevation spectra and optionally the direct
    full FK+diffraction spectral wrench for comparison with prepared loads.
    This does not change those prepared body loads.
    """
    from topwater.solve.diffraction import packet_diffraction
    checkpoint = (lambda: None) if checkpoint is None else checkpoint
    xy = exterior_field_points(hydro, xy, allow_boundary=allow_boundary)
    omega = _grid(train.omega, 'train.omega', zero=True)
    if omega[0] != 0 or np.any(train.C[:, :, 0] != 0):
        raise ValueError('A zero-DC signed-kappa train is required')
    incident = incident_spectrum(train, xy) if return_incident else np.zeros((len(omega), len(xy)), complex)
    scattered = np.zeros_like(incident)
    force = np.zeros((len(omega), 6), complex)
    response = packet_diffraction(hydro, train)
    progress = _StageProgress('diffraction transfers', len(omega)-1, 'positive sidebands')
    for iw in range(1, len(omega)):
        checkpoint()
        result = response.solve(iw)
        if result is not None:
            scattered[iw] = field_backend.free_surface_elevation(xy, result)
        if return_forces:
            force[iw] = response.force(iw)
        progress.update(iw)
    progress.finish()
    return (incident, scattered, force) if return_forces else (incident, scattered)


def synthesize_spectrum(omega, coefficients, time_s):
    """Synthesize already quadrature-weighted positive-frequency coefficients."""
    omega = _grid(omega, 'omega', zero=True)
    values = _finite(coefficients, (len(omega), None), 'coefficients', complex_value=True)
    times = _finite(time_s, (None,), 'time_s')
    return (np.exp(-1j*times[:, None]*omega) @ values).real




def _regular_grid(bounds, spacing):
    lower, upper = np.asarray(bounds[0]), np.asarray(bounds[1])
    axes = [np.linspace(a, b, int(np.ceil((b-a)/spacing))+1) for a, b in zip(lower, upper)]
    x, y = np.meshgrid(*axes, indexing='xy')
    xy = np.c_[x.ravel(), y.ravel()]
    ids = np.arange(x.size).reshape(x.shape)
    a, b, c, d = ids[:-1, :-1], ids[:-1, 1:], ids[1:, :-1], ids[1:, 1:]
    faces = np.concatenate([np.stack([a, b, d], axis=-1).reshape(-1, 3),
                            np.stack([a, d, c], axis=-1).reshape(-1, 3)])
    return xy, faces, axes


def reconstruct_total_surface(hydro, setup, motion, surface_time_s, *, radiation_settings=None, field_backend=None, surface_options=None):
    """Run the shared reconstruction pipeline with owned or supplied resources."""
    if field_backend is not None:
        return _reconstruct_total_surface(hydro, setup, motion, surface_time_s,
            radiation_settings=radiation_settings, field_backend=field_backend, surface_options=surface_options)
    from topwater.solve.field_backend import select_field_backend
    backend = select_field_backend(hydro, 'cpu', cpu_workers=None if surface_options is None else surface_options.cpu_workers)
    try:
        return _reconstruct_total_surface(hydro, setup, motion, surface_time_s,
            radiation_settings=radiation_settings, field_backend=backend, surface_options=surface_options)
    finally:
        primary = sys.exception()
        try:
            backend.close()
        except Exception as error:
            if primary is None:
                raise
            primary.add_note(f'Field cleanup also failed: {error}')


def _reconstruct_total_surface(hydro, setup, motion, surface_time_s, *, radiation_settings, field_backend, surface_options):
    """Return packed original total-field grids; hull interiors are invalid water.

    Numerical policy: fixed mean-hull linear diffraction/radiation, positive
    sideband frequencies, intact finite-band offline complex response over the complete
    trajectory, exact incident at targets, explicit exterior-trace continuation
    of disturbances near/inside the mean footprint. This is an export of a
    selected linear reference approximation, not nonlinear contact dynamics.
    No body pose, velocity, load or hydrodynamic coefficient is modified.
    """
    from scipy.spatial.transform import Rotation, Slerp
    from threadpoolctl import threadpool_limits
    from topwater.solve.settings import DEFAULTS, VIDEO_SIZE_PX, SURFACE_PIXELS_PER_SAMPLE
    from topwater.solve.settings import VIDEO_SIZE_PX
    options = DEFAULTS if surface_options is None else surface_options

    stamp = _grid(motion.time_s, 'motion.time_s', zero=True)
    surface_time = _grid(surface_time_s, 'surface_time_s', zero=True)
    if stamp[0] != 0 or surface_time[0] != 0 or stamp[-1] != surface_time[-1]:
        raise ValueError('Body and surface histories must share exact zero/start/end times')
    if not np.allclose(np.diff(stamp), stamp[1], rtol=0, atol=1e-12):
        raise ValueError('Radiation integration requires the uniform native time grid')
    if hydro.depth != np.inf:
        raise ValueError('The current incoming-train total-surface policy requires infinite depth')
    velocity = _finite(motion.velocity_native, (len(stamp), 6), 'native radiation-input velocity')
    acceleration = _finite(motion.acceleration_native, velocity.shape, 'native radiation-input acceleration')
    if radiation_settings is None:
        radiation_settings = options.radiation_settings(field_backend.radiation_chunk_points)
    position = _finite(motion.position_m, (len(stamp), 3), 'world position')
    rotation = Rotation.from_quat(motion.orientation_wxyz, scalar_first=True)
    matrices = rotation.as_matrix()
    poses = Slerp(stamp, rotation)(surface_time).as_matrix()
    positions = np.column_stack([np.interp(surface_time, stamp, position[:, j]) for j in range(3)])
    body = setup.body
    local_hull = body.mesh.vertices-body.cg
    hull_world = np.einsum('tij,vj->tvi', matrices, local_hull)+position[:, None, :]
    embedding = _finite(setup.world_from_native, (3, 3), 'world_from_native')
    beam = float(np.ptp(body.mesh.vertices[:, 1]))
    if beam <= 0:
        raise ValueError('A resolved boat beam is required for surface export')
    physical_ids = np.unique(setup.physical_faces)
    physical = hull_world[:, physical_ids, :2]
    view_lo, view_hi = physical.min(axis=(0, 1)), physical.max(axis=(0, 1))
    view_center = (view_lo+view_hi)/2
    view_half = float(np.max(view_hi-view_lo)/2+beam*options.surface_padding_beams)
    view_bounds = np.array([view_center-view_half, view_center+view_half])
    # Cover the solver's world view even when the native axes are oblique.
    corners = np.array([[view_center[0]+sx*view_half, view_center[1]+sy*view_half, 0.]
                        for sx in (-1, 1) for sy in (-1, 1)]) @ embedding
    spacing = surface_spacing_m(view_bounds)
    lo, hi = corners[:, :2].min(axis=0)-spacing, corners[:, :2].max(axis=0)+spacing
    center = (lo+hi)/2
    half = float(np.max(hi-lo)/2)
    xy, faces, axes = _regular_grid([center-half, center+half], spacing)
    mapped, continued = exterior_trace_points(hydro, xy, 0., allow_boundary=True)
    field_xy, inverse = np.unique(mapped, axis=0, return_inverse=True)
    print(f'Surface: evaluating total linear field at {len(field_xy)} exterior/boundary points; '
          f'{int(continued.sum())} target points use the exact mean-waterline disturbance trace.', file=sys.stderr)
    print('Surface numerical approximation: intact complex offline response; '
          f'temporal policy={"native unpadded" if radiation_settings.period_s is None else "explicit historical padded"}; '
          'all velocity FFT bins. Raw startup is preserved. '
          'Whole-domain numerical accuracy and remaining-tail uncertainty require separate qualification.', file=sys.stderr)
    if not isinstance(radiation_settings.chunk_points, int) or radiation_settings.chunk_points < 1:
        raise ValueError('Radiation spatial chunk size must be a positive integer')
    omega = surface_radiation_frequencies(stamp, radiation_settings)
    chunk_points = field_backend.prepare_radiation(stamp, surface_time, omega, radiation_settings, len(field_xy))
    radiation_evidence = []
    with threadpool_limits(1):
        progress = _StageProgress('offline radiation reconstruction', len(field_xy), 'field points')
        radiation = np.empty((len(surface_time), len(field_xy)))
        for start in range(0, len(field_xy), chunk_points):
            end = min(start+chunk_points, len(field_xy))
            # radiation_results retains the source bank across spatial chunks;
            # only native field evaluation repeats, never the body/source solve.
            transfer = field_backend.radiation_transfer(hydro, field_xy[start:end], omega[1:])
            values = field_backend.prepend_zero(transfer)
            response, evidence = offline_radiation_response(omega, values, stamp, velocity, acceleration,
                output_time_s=surface_time, settings=radiation_settings, workspace=field_backend.fourier)
            radiation[:, start:end] = response
            evidence.update(field_point_start=start, field_point_end=end)
            radiation_evidence.append(evidence)
            print('SURFACE_NUMERICS '+json.dumps(evidence, allow_nan=False, separators=(',', ':')), file=sys.stderr)
            del transfer, values
            progress.update(end)
        progress.finish()
        field_backend.finish_radiation()
        _, scattered_spectrum = combined_diffraction(
            hydro, setup.train, field_xy, return_forces=False, return_incident=False,
            allow_boundary=True, field_backend=field_backend)
        progress = _StageProgress('field synthesis', 1, 'completed field')
        scattered = synthesize_spectrum(setup.train.omega, scattered_spectrum, surface_time)
        del scattered_spectrum
        incident_coefficients = field_backend.incident_spectrum(setup.train, xy)
        incident = synthesize_spectrum(setup.train.omega, incident_coefficients, surface_time)
        total = incident+scattered[:, inverse]+radiation[:, inverse]
    progress.finish()
    # Preserve the expensive completed field before validation/packing. The caller's private OS-temp context removes it on success and
    # retains it on SurfaceExportError for a geometry-only replay.
    completed = dict(cache_version=np.array(1, dtype=np.int64),
        producer_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        surface_time_s=surface_time, surface_position_m=positions, surface_rotation=poses,
        motion_time_s=stamp, motion_position_m=position,
        motion_orientation_wxyz=np.asarray(motion.orientation_wxyz), motion_velocity_native=velocity,
        body_vertices_m=body.mesh.vertices, body_cg_m=body.cg, body_faces=body.mesh.faces,
        opening_rim=body.mesh.opening_rim, artificial_faces=body.mesh.artificial_faces,
        physical_faces=setup.physical_faces, xy_native_m=xy, water_faces=faces,
        axis_x_m=axes[0], axis_y_m=axes[1], total_elevation_m=total, incident_elevation_m=incident,
        radiation_numerics_json=np.array(json.dumps(radiation_evidence, allow_nan=False)),
        world_from_native=embedding, mean_rotation=Rotation.from_quat(
            setup.equilibrium.orientation, scalar_first=True).as_matrix(),
        mean_cg_m=setup.mean_cg_m, initial_wetted_area_m2=np.array(setup.hydrostatic.wetted_area),
        beam_m=np.array(beam),
        world_vertical_offset_m=np.array(float(getattr(setup, 'world_vertical_offset_m', 0.))),
        calm_eta_limit_m=np.array(float(getattr(setup.train, 'envelope_tolerance_m',
            .001*float(np.sum(setup.train.carriers.amplitudes))))))
    if hasattr(motion, 'acceleration_native'):
        completed['motion_acceleration_native'] = np.asarray(motion.acceleration_native)
    if hasattr(motion, 'q'):
        completed['motion_q_native'] = np.asarray(motion.q)
    completed['view_bounds_xy_m'] = view_bounds
    completed['video_size_px'] = np.asarray(VIDEO_SIZE_PX)
    if surface_options is not None:
        completed['surface_interval_s'] = np.array(options.surface_interval_s)
    directory = Path(hydro.directory).parent/'surface-recovery'
    cache = _save_completed_cache(directory, completed)
    return SurfaceOutput(*_export_completed_field(completed, directory, cache), view_bounds)
