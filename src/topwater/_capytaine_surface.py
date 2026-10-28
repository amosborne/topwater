"""Native Capytaine coefficients and reviewed building blocks for surface export.

Capytaine is imported only when a producer is explicitly called.  All harmonic
quantities use Re(A exp(-i omega t)), fixed mean-world axes and the mean CG.
Coefficient production does not qualify transient field reconstruction.
The full-complex offline approximation requires explicit saved-case numerical
qualification; historical causal-kernel routines remain for preserved studies.
"""
from dataclasses import dataclass
from functools import wraps
import hashlib
from importlib.metadata import version
import json
import logging
from pathlib import Path
import sys
from time import perf_counter
import traceback

import numpy as np
from scipy.integrate import trapezoid
from scipy.fft import rfft, irfft, next_fast_len
from scipy.spatial import ConvexHull, QhullError


MODES = ('Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw')


class _StageProgress:
    """Wall-time progress at existing work boundaries; never a physical clock."""
    def __init__(self, name, total, unit):
        self.name, self.total, self.unit = name, int(total), unit
        self.started = self.last_report = perf_counter()
        print(f'STAGE {name} started total={self.total} unit={unit}', file=sys.stderr)

    def update(self, completed):
        now = perf_counter()
        if now-self.last_report >= 5.:
            print(f'PROGRESS {self.name} completed={completed}/{self.total} '
                  f'unit={self.unit} elapsed_s={now-self.started:.1f}', file=sys.stderr)
            self.last_report = now

    def finish(self):
        print(f'STAGE {self.name} complete completed={self.total}/{self.total} '
              f'unit={self.unit} elapsed_s={perf_counter()-self.started:.1f}', file=sys.stderr)


class SurfaceExportError(RuntimeError):
    """Primary geometry/metric failure with a private replayable field cache."""
    def __init__(self, message, recovery_directory, *, stage, frame_index, time_s):
        super().__init__(message)
        self.recovery_directory = Path(recovery_directory)
        self.stage = stage
        self.frame_index = int(frame_index)
        self.time_s = None if time_s is None else float(time_s)


def _save_completed_cache(directory, data):
    """Private disposable intermediate, never placed in the example folder."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    path = directory/'completed-field.npz'
    with path.open('xb') as stream:
        np.savez(stream, **data)
    return path


def _save_frame_failure(directory, evidence, error, stage, index, stamp, cache_path):
    """Best-effort secondary diagnostics must never replace the primary error."""
    directory = Path(directory)
    report = dict(stage=stage, frame_index=int(index), time_s=(None if stamp is None else float(stamp)),
                  error_type=type(error).__name__, error_message=str(error),
                  traceback=traceback.format_exc(),
                  completed_cache=str(cache_path))
    try:
        directory.mkdir(parents=True, exist_ok=True)
        with (directory/'failure.json').open('x', encoding='utf-8') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
        with (directory/'failure-frame.npz').open('xb') as stream:
            np.savez(stream, **evidence)
    except Exception as diagnostic_error:
        print(f'Surface failure evidence could not be fully saved: {diagnostic_error}; '
              'the original failure remains authoritative.', file=sys.stderr)


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
        from .workflow_data import Result
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
            surface_vertex_offsets=packed[2], surface_face_offsets=packed[3])
        result.save(output_path, overwrite=False)
    return packed


class _CapytaineWarningFilter(logging.Filter):
    """Keep the first physical warning of each known family, then summarize."""
    def __init__(self):
        super().__init__()
        self.counts = {}

    def filter(self, record):
        text = record.getMessage()
        family = next((prefix for prefix in ('Mesh resolution', 'Irregular frequencies', 'Water depth')
                       if text.startswith(prefix)), None)
        if family is None:
            return True
        self.counts[family] = self.counts.get(family, 0)+1
        return self.counts[family] == 1


def _summarize_capytaine_warnings(function):
    @wraps(function)
    def call(*args, **kwargs):
        logger = logging.getLogger('capytaine.bem.problems_checks')
        summary = _CapytaineWarningFilter()
        logger.addFilter(summary)
        try:
            return function(*args, **kwargs)
        finally:
            logger.removeFilter(summary)
            for family, count in summary.counts.items():
                if count > 1:
                    print(f'Capytaine: {family} warning applies to {count} solves in '
                          f'{function.__name__}; first warning retained. These limitations remain unqualified.',
                          file=sys.stderr)
    return call


def _single_threaded(function):
    """Set the existing one-thread policy once, not once per scalar BEM solve."""
    @wraps(function)
    def call(*args, **kwargs):
        from threadpoolctl import threadpool_limits
        with threadpool_limits(1):
            return function(*args, **kwargs)
    return call


def _finite(value, shape, name, *, complex_value=False):
    result = np.asarray(value, dtype=complex if complex_value else float)
    if result.ndim != len(shape) or any(n is not None and n != m for n, m in zip(shape, result.shape)):
        raise ValueError(f'{name} must have shape {shape}')
    if not np.isfinite(result).all():
        raise ValueError(f'{name} must be finite')
    return result


def _grid(values, name, *, zero=False):
    values = _finite(values, (None,), name)
    if len(values) < 2 or np.any(np.diff(values) <= 0) or values[0] < 0 or (not zero and values[0] == 0):
        raise ValueError(f'{name} must increase strictly from a {"nonnegative" if zero else "positive"} value')
    return values


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


@dataclass
class Hydrodynamics:
    directory: Path
    netcdf_path: Path
    native_mat_path: Path
    radiation_omega: np.ndarray
    excitation_omega: np.ndarray
    excitation_headings: np.ndarray
    added_mass: np.ndarray
    radiation_damping: np.ndarray
    Ainf: np.ndarray
    excitation: np.ndarray
    M: np.ndarray
    K: np.ndarray
    cg: np.ndarray
    cb: np.ndarray
    volume: float
    body: object
    solver: object
    radiation_results: dict
    diffraction_results: dict
    density: float
    gravity: float
    depth: float

    @property
    def omega(self):
        return self.excitation_omega

    @property
    def headings(self):
        return self.excitation_headings


def _native_excitation_payload(omega, headings, froude_krylov, diffraction, total, stiffness):
    """Dense native coefficients for the standard directional WEC-Sim H5.

    Arrays are dimensional (N or N*m per metre elevation), ordered
    (omega, heading, Surge/Sway/Heave/Roll/Pitch/Yaw), and use exp(-i*omega*t).
    Incident phase is referenced to native (0,0,0) and absolute time zero;
    rotational modes are about the mean CG. The consumer performs the single
    WEC-Sim phase/normalization conversion. No wave history or force is built.

    The raw table starts at 0.1 rad/s. Its separate constant-pressure limit
    preserves the declared quadratic bridge below that frequency; it is not a
    fabricated zero-frequency BEM solve. Radiation's wider grid is unchanged.
    """
    frequencies = _grid(omega, 'excitation_omega_rad_s')
    angles = _grid(headings, 'excitation_headings_rad', zero=True)
    if abs(frequencies[0]-.1) > 1e-12:
        raise ValueError('The declared excitation bridge requires a 0.1 rad/s first coefficient')
    if angles[0] != 0 or abs(angles[-1]-2*np.pi) > 1e-12:
        raise ValueError('Native excitation headings must include the 0/2pi periodic seam')
    shape = (len(frequencies), len(angles), 6)
    fk = _finite(froude_krylov, shape, 'excitation_froude_krylov_native', complex_value=True)
    sc = _finite(diffraction, shape, 'excitation_diffraction_native', complex_value=True)
    excitation = _finite(total, shape, 'excitation_native', complex_value=True)
    for name, values in (('FK', fk), ('diffraction', sc), ('total', excitation)):
        if not np.array_equal(values[:, 0], values[:, -1]):
            raise ValueError(f'Native {name} excitation periodic seam does not close')
    if not np.array_equal(excitation, fk+sc):
        raise ValueError('Native total excitation must equal its genuine FK plus diffraction coefficients')
    K = _finite(stiffness, (6, 6), 'hydrostatic_stiffness')
    return dict(excitation_omega_rad_s=frequencies,
                excitation_headings_rad=angles,
                excitation_froude_krylov_native=fk,
                excitation_diffraction_native=sc,
                excitation_native=excitation,
                excitation_zero_limit_native=K[:, 2])


@_summarize_capytaine_warnings
def build_hydrodynamics(wet_triangles, cg, mass_matrix, workspace, *,
                        hydrostatic_stiffness, cb, volume,
                        density=1025., gravity=9.81, depth=np.inf,
                        radiation_omega=None, excitation_omega=None, headings=None,
                        checkpoint=None, green_cache_dir=None):
    """Produce the original coefficient policies, retaining every source solve.

    ``wet_triangles`` are the frozen physical mean wetted surface, not a hull
    regenerated by this function.  ``K``, ``cb`` and ``volume`` come from the
    exact hydrostatic setup.  No mesh refinement or fitted coefficient repair
    occurs.  The official import netCDF retains the original three headings;
    native.mat separately carries genuine dense FK/diffraction coefficients
    with a periodic full-heading seam for native WEC-Sim excitation preparation.

    This function performs numerical work. Its files belong in per-run scratch
    space; do not share a mutable workspace between concurrent producers.
    """
    if version('capytaine') != '3.0.0':
        raise ValueError('The reviewed coefficient producer requires Capytaine 3.0.0')
    import capytaine as cpt
    import xarray as xr
    from capytaine.green_functions.delhommeau import Delhommeau
    from capytaine.io.legacy import export_hydrostatics
    from capytaine.bem.airy_waves import froude_krylov_force
    from scipy.io import savemat
    from threadpoolctl import threadpool_limits
    checkpoint = (lambda: None) if checkpoint is None else checkpoint

    triangles = _finite(wet_triangles, (None, 3, 3), 'wet_triangles')
    center = _finite(cg, (3,), 'cg')
    M = _finite(mass_matrix, (6, 6), 'mass_matrix')
    K = _finite(hydrostatic_stiffness, (6, 6), 'hydrostatic_stiffness')
    cb = _finite(cb, (3,), 'cb')
    if not len(triangles) or not np.isfinite([density, gravity, volume]).all() or min(density, gravity, volume) <= 0:
        raise ValueError('Positive density, gravity, displaced volume and wet triangles are required')
    if not (depth == np.inf or np.isfinite(depth) and depth > 0):
        raise ValueError('depth must be positive or infinite')
    rad_w = _grid(np.arange(1, 121) / 10 if radiation_omega is None else radiation_omega, 'radiation_omega')
    exc_w = _grid(np.arange(1, 35) / 10 if excitation_omega is None else excitation_omega, 'excitation_omega')
    angles = _grid(np.linspace(0, 2*np.pi, 25) if headings is None else headings, 'headings', zero=True)
    if angles[0] != 0 or abs(angles[-1] - 2*np.pi) > 1e-12:
        raise ValueError('headings must include the periodic 0/2pi seam')
    if not all(np.any(np.isclose(angles[:-1], a, rtol=0, atol=1e-12)) for a in (0., np.pi/4, np.pi/2)):
        raise ValueError('The original import headings 0/45/90 degrees must be represented')
    directory = Path(workspace)
    directory.mkdir(parents=True, exist_ok=False)
    # Match the established Capytaine triangle convention without changing any
    # facet, winding or collocation/centroid quadrature.
    starts = np.argmax(np.linalg.norm(triangles-triangles.mean(axis=1)[:, None], axis=2), axis=1)
    triangles = np.array([np.roll(t, -int(i), axis=0) for t, i in zip(triangles, starts)])
    vertices, ids = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
    mesh = cpt.Mesh(vertices=vertices, faces=ids.reshape(-1, 3), auto_clean=False, auto_check=False)
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=center),
                            center_of_mass=center, mass=float(M[0, 0]), name='boat')
    labels = dict(influenced_dof=list(MODES), radiating_dof=list(MODES))
    body.inertia_matrix = xr.DataArray(M, dims=list(labels), coords=labels)
    body.hydrostatic_stiffness = xr.DataArray(K, dims=list(labels), coords=labels)
    # A caller-owned job-temp cache can be shared by producer instances.
    # None preserves the original Green-function settings and no global cache.
    cache = None if green_cache_dir is None else str(Path(green_cache_dir))
    if cache is not None:
        Path(cache).mkdir(parents=True, exist_ok=True)
    tabulation = _StageProgress('Green-function tabulation', 1, 'table')
    solver = cpt.BEMSolver(green_function=Delhommeau(tabulation_cache_dir=cache))
    tabulation.finish()
    radiation, diffraction, imported = {}, {}, []
    A = np.empty((len(rad_w), 6, 6)); B = np.empty_like(A)
    Ainf = np.empty((6, 6))
    H = np.empty((len(exc_w), len(angles), 6), dtype=complex)
    FK = np.empty_like(H)
    SC = np.empty_like(H)
    common = dict(body=body, rho=float(density), g=float(gravity), water_depth=depth)
    import_angles = (0., np.pi/4, np.pi/2)
    solve_count = 0
    solve_total = (6+6*len(rad_w)+len(exc_w)*(len(angles)-1)
                   +3*len(np.setdiff1d(rad_w, exc_w)))
    progress = _StageProgress('coefficient preparation', solve_total, 'BEM solves')
    with threadpool_limits(1):
        for j, mode in enumerate(MODES):
            checkpoint()
            result = solver.solve(cpt.RadiationProblem(**common, omega=np.inf, radiating_dof=mode), keep_details=True)
            Ainf[:, j] = [result.added_mass[m] for m in MODES]
            radiation[(np.inf, mode)] = result
            solve_count += 1
            progress.update(solve_count)
        all_w = np.unique(np.r_[rad_w, exc_w])
        for w in all_w:
            checkpoint()
            rad_index = np.flatnonzero(rad_w == w)
            exc_index = np.flatnonzero(exc_w == w)
            if len(rad_index):
                for j, mode in enumerate(MODES):
                    checkpoint()
                    result = solver.solve(cpt.RadiationProblem(**common, omega=float(w), radiating_dof=mode), keep_details=True)
                    radiation[(float(w), mode)] = result
                    imported.append(result)
                    A[rad_index[0], :, j] = [result.added_mass[m] for m in MODES]
                    B[rad_index[0], :, j] = [result.radiation_damping[m] for m in MODES]
                    solve_count += 1
                    progress.update(solve_count)
            query_angles = angles[:-1] if len(exc_index) else np.asarray(import_angles)
            for a in query_angles:
                checkpoint()
                result = solver.solve(cpt.DiffractionProblem(**common, omega=float(w), wave_direction=float(a)), keep_details=True)
                diffraction[(float(w), float(a))] = result
                if len(rad_index) and any(abs(a-b) < 1e-12 for b in import_angles):
                    imported.append(result)
                if len(exc_index):
                    j = int(np.flatnonzero(angles[:-1] == a)[0])
                    fk = froude_krylov_force(result.problem)
                    H[exc_index[0], j] = [result.forces[m] + fk[m] for m in MODES]
                    FK[exc_index[0], j] = [fk[m] for m in MODES]
                    SC[exc_index[0], j] = [result.forces[m] for m in MODES]
                solve_count += 1
                progress.update(solve_count)
            if len(exc_index):
                H[exc_index[0], -1] = H[exc_index[0], 0]
                FK[exc_index[0], -1] = FK[exc_index[0], 0]
                SC[exc_index[0], -1] = SC[exc_index[0], 0]
        data = cpt.assemble_dataset(imported, hydrostatics=True).sel(influenced_dof=list(MODES), radiating_dof=list(MODES))
        data['hydrostatic_stiffness'].data = K
        data['center_of_buoyancy'].data = cb
        data['disp_mass'].data = np.asarray(density*volume)
        if 'body_name' not in data:
            data['body_name'] = data['body'].copy()
        with xr.set_options(netcdf_engine_order=['scipy']):
            cpt.export_dataset(directory/'canonical.nc', data)
        export_hydrostatics(str(directory), body)
    for name, value in [('A', A), ('B', B), ('Ainf', Ainf), ('excitation', H)]:
        if not np.isfinite(value).all():
            raise ValueError(f'Incomplete native {name}')
    excitation_payload = _native_excitation_payload(exc_w, angles, FK, SC, H, K)
    savemat(directory/'native.mat', dict(Ainf=Ainf, M=M, K=K, cg=center, cb=cb,
             volume=float(volume), rho=float(density), g=float(gravity), **excitation_payload))
    progress.finish()
    return Hydrodynamics(directory, directory/'canonical.nc', directory/'native.mat', rad_w, exc_w,
                         angles, A, B, Ainf, H, M, K, center, cb, float(volume), body,
                         solver, radiation, diffraction, float(density), float(gravity), depth)


@_summarize_capytaine_warnings
@_single_threaded
def combined_diffraction(hydro, train, xy, *, return_forces=False, checkpoint=None, return_incident=True, allow_boundary=False):
    """Exact sideband combined incident/scattered spectra on the frozen hull.

    ``train.C`` already includes spectral quadrature weights.  Its two sign
    branches are two positive-frequency problems with opposite propagation
    headings; no negative-frequency outgoing Green function is evaluated.
    Returns complex (omega, point) elevation spectra and optionally the direct
    full FK+diffraction spectral wrench for comparison with prepared loads.
    This does not change those prepared body loads.
    """
    import capytaine as cpt
    from capytaine.bem.problems_and_results import LinearPotentialFlowProblem
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation, froude_krylov_force
    checkpoint = (lambda: None) if checkpoint is None else checkpoint
    xy = exterior_field_points(hydro, xy, allow_boundary=allow_boundary)
    omega = _grid(train.omega, 'train.omega', zero=True)
    C = _finite(train.C, (None, 2, len(omega)), 'train.C', complex_value=True)
    theta = _finite(train.theta, (len(C),), 'train.theta')
    incident = np.zeros((len(omega), len(xy)), dtype=complex)
    scattered = np.zeros_like(incident)
    force = np.zeros((len(omega), 6), dtype=complex)
    if omega[0] != 0 or np.any(C[:, :, 0] != 0):
        raise ValueError('A zero-DC signed-kappa train is required')
    common = dict(body=hydro.body, rho=hydro.density, g=hydro.gravity, water_depth=hydro.depth)
    progress = _StageProgress('diffraction transfers', len(omega)-1, 'positive sidebands')
    for iw in range(1, len(omega)):
        checkpoint()
        w = float(omega[iw])
        rhs = np.zeros(hydro.body.mesh_including_lid.nb_faces, dtype=complex)
        fk = np.zeros(6, dtype=complex)
        for j in range(len(C)):
            for branch in range(2):
                coefficient = C[j, branch, iw]
                if coefficient == 0:
                    continue
                problem = cpt.DiffractionProblem(**common, omega=w,
                    wave_direction=float((theta[j]+branch*np.pi) % (2*np.pi)))
                rhs += coefficient * problem.boundary_condition
                if return_incident:
                    incident[iw] += coefficient * airy_waves_free_surface_elevation(xy, problem)
                if return_forces:
                    values = froude_krylov_force(problem)
                    fk += coefficient*np.array([values[m] for m in MODES])
        # The custom boundary condition is the exact spectral linear combination
        # of all incident carriers, not a unit plane wave at the nominal angle.
        problem = LinearPotentialFlowProblem(**common, omega=w, boundary_condition=rhs)
        result = hydro.solver.solve(problem, keep_details=True)
        scattered[iw] = hydro.solver.compute_free_surface_elevation(xy, result)
        if return_forces:
            force[iw] = fk + np.array([result.forces[m] for m in MODES])
        progress.update(iw)
    progress.finish()
    return (incident, scattered, force) if return_forces else (incident, scattered)


def synthesize_spectrum(omega, coefficients, time_s):
    """Synthesize already quadrature-weighted positive-frequency coefficients."""
    omega = _grid(omega, 'omega', zero=True)
    values = _finite(coefficients, (len(omega), None), 'coefficients', complex_value=True)
    times = _finite(time_s, (None,), 'time_s')
    return (np.exp(-1j*times[:, None]*omega) @ values).real


@_summarize_capytaine_warnings
@_single_threaded
def radiation_velocity_transfer(hydro, xy, omega=None, *, checkpoint=None, shared_green=False, allow_boundary=False):
    """Return H_eta,v=i H_eta,q/omega, shape (omega, point, six modes)."""
    import capytaine as cpt
    checkpoint = (lambda: None) if checkpoint is None else checkpoint
    xy = exterior_field_points(hydro, xy, allow_boundary=allow_boundary)
    frequencies = hydro.radiation_omega if omega is None else _grid(omega, 'omega')
    values = np.empty((len(frequencies), len(xy), 6), dtype=complex)
    if shared_green:
        # DefaultMatrixEngine.build_S_matrix merges (copies) this immutable
        # source mesh on every call, discarding all cached geometric properties.
        # Use its exact public Green-function evaluator on one merged mesh.
        # No quadrature, field points or source coefficients change.
        source_mesh = hydro.body.mesh_including_lid.merged()
        field_points = np.c_[xy, np.zeros(len(xy))]
    progress = _StageProgress('radiation transfers', len(frequencies), 'frequencies')
    for iw, w in enumerate(frequencies):
        results = []
        for j, mode in enumerate(MODES):
            checkpoint()
            key = (float(w), mode)
            if key not in hydro.radiation_results:
                hydro.radiation_results[key] = hydro.solver.solve(cpt.RadiationProblem(
                    body=hydro.body, rho=hydro.density, g=hydro.gravity,
                    water_depth=hydro.depth, omega=float(w), radiating_dof=mode), keep_details=True)
            results.append(hydro.radiation_results[key])
            if not shared_green:
                eta_q = hydro.solver.compute_free_surface_elevation(xy, hydro.radiation_results[key])
                values[iw, :, j] = 1j*eta_q/float(w)
        if shared_green:
            # This is precisely the public field evaluator's S @ sources,
            # batching the six source columns. No new BEM equation or solver.
            representative = results[0]
            S, _ = hydro.solver.engine.green_function.evaluate(field_points,
                source_mesh, free_surface=representative.free_surface,
                water_depth=representative.water_depth, wavenumber=representative.encounter_wavenumber)
            if not np.isfinite(S).all():
                raise ValueError('Nonfinite Capytaine single-layer field matrix')
            sources = np.column_stack([r.sources for r in results])
            # i/omega * (i*omega/g) * potential = -potential/g.
            values[iw] = -(S @ sources)/hydro.gravity
        progress.update(iw+1)
    progress.finish()
    return values


def radiation_acceleration_direct(hydro, xy):
    """Exact discrete infinite-frequency acceleration-to-elevation term.

    Capytaine represents the endpoint symbolically. Its public elevation result
    divided by omega² is finite; eta=E_inf*omega²*q=-E_inf*acceleration.
    This identifies the discrete endpoint term, not its spatial accuracy.
    """
    xy = exterior_field_points(hydro, xy)
    direct = np.empty((len(xy), 6))
    for j, mode in enumerate(MODES):
        result = hydro.radiation_results[(np.inf, mode)]
        elevation = hydro.solver.compute_free_surface_elevation(xy, result)
        endpoint = np.asarray(elevation/(result.omega*result.omega), dtype=complex)
        if np.max(abs(endpoint.imag)) > 1e-12*max(1., float(np.max(abs(endpoint.real)))):
            raise ValueError('Infinite-frequency surface acceleration term is not real')
        direct[:, j] = -endpoint.real
    return direct


@dataclass(frozen=True)
class OfflineRadiationSettings:
    """Explicit finite-band OFFLINE policy, requiring case-specific qualification.

    Defaults use exactly the saved native samples, without padding or resampling.
    An explicit period_s selects the historical padded policy; its optional
    future polynomial is numerical continuation, not predicted body motion.
    """
    period_s: float | None = None
    fft_dt_s: float | None = None
    future_s: float = 0.
    maximum_omega_rad_s: float = 48.
    omega_step_rad_s: float = .025
    chunk_points: int = 128


def offline_radiation_response(omega, velocity_transfer, time_s, velocity, acceleration,
                               *, output_time_s=None, settings=None):
    """Full complex H_eta,v applied to the complete saved native velocity.

    H uses exp(-i*w*t); NumPy's inverse uses exp(+i*w*t), hence conj(H).
    By default the transform contains only the saved native velocity samples.
    Explicit period_s settings retain the historical padded/Hermite policy.
    Frequencies above the declared band are zero. These choices are logged,
    NOT a universal causality claim
    or an automatic numerical-error certificate. A*a is already part of full H:
    never add it a second time or force its first sample to zero.
    """
    settings = OfflineRadiationSettings() if settings is None else settings
    stamp = _grid(time_s, 'motion time_s', zero=True)
    w = _grid(omega, 'transfer omega', zero=True)
    H = _finite(velocity_transfer, (len(w), None, 6), 'velocity_transfer', complex_value=True)
    v = _finite(velocity, (len(stamp), 6), 'native velocity')
    a = _finite(acceleration, v.shape, 'native acceleration')
    output = stamp if output_time_s is None else _grid(output_time_s, 'output time_s', zero=True)
    native_dt = stamp[1]-stamp[0]
    native = settings.period_s is None
    dt = (native_dt if native else native_dt/2) if settings.fft_dt_s is None else float(settings.fft_dt_s)
    future, band = settings.future_s, settings.maximum_omega_rad_s
    if (not np.isfinite([dt, future, band]).all() or min(dt, band) <= 0 or future < 0
            or stamp[0] != 0 or w[0] != 0 or np.any(v[0] != 0)
            or not np.allclose(np.diff(stamp), native_dt, rtol=0, atol=1e-12)):
        raise ValueError('Offline radiation requires finite settings, uniform zero-start history and exact initial velocity rest')
    if native:
        if abs(dt-native_dt) > 1e-12 or future != 0:
            raise ValueError('Native FFT requires saved timestep and no future continuation')
        n, zero = len(stamp), 0
        period = n*dt  # DFT period includes one sample interval beyond the last knot.
        selected = np.rint(output/dt).astype(int)
        if (band > np.pi/dt or w[-1] < band or output[0] < 0 or output[-1] > stamp[-1]):
            raise ValueError('Offline FFT band must cover history without frequency extrapolation')
        if not np.allclose(selected*dt, output, rtol=0, atol=1e-10):
            raise ValueError('Output clocks must be exactly represented by the native FFT dt')
        values = v  # No interpolation, oversampling, padding or invented future.
    else:
        period = float(settings.period_s)
        if not np.isfinite(period) or period <= 0:
            raise ValueError('Explicit padded FFT period must be positive and finite')
        n = int(round(period/dt))
        if (n % 4 or abs(n*dt-period) > 1e-10 or band > np.pi/dt
                or w[-1] < band or output[0] < 0 or output[-1] > stamp[-1]):
            raise ValueError('Offline FFT period/grid/band must cover the requested history without frequency extrapolation')
        knot_indices = np.rint(stamp/dt).astype(int)
        selected = np.rint(output/dt).astype(int)
        steps = int(round(future/dt))
        if (not np.allclose(knot_indices*dt, stamp, rtol=0, atol=1e-10)
                or not np.allclose(selected*dt, output, rtol=0, atol=1e-10)
                or steps < 0 or abs(steps*dt-future) > 1e-10):
            raise ValueError('Saved motion, output clocks and future duration must be exactly represented by FFT dt')
        zero = n//4
        last = int(knot_indices[-1])
        if zero+last+steps >= n:
            raise ValueError('FFT period does not fit the saved history and declared smooth future')
        values = np.zeros((n, 6))
        inside = np.arange(last+1)*dt
        for mode in range(6):
            values[zero:zero+last+1, mode] = np.interp(inside, stamp, v[:, mode])
        values[zero+knot_indices] = v  # Preserve every native sample exactly.
        if not np.array_equal(values[zero+knot_indices], v):
            raise ValueError('Offline interpolation changed saved velocity samples')
        u = np.arange(1, steps+1)/steps if steps else np.empty(0)
        values[zero+last+1:zero+last+steps+1] = (
            (2*u**3-3*u**2+1)[:, None]*v[-1]
            +future*(u**3-2*u**2+u)[:, None]*a[-1])
    # Accumulate modes into one output spectrum: no FFT-frequency x point x DOF
    # tensor. Spatial chunking by the caller bounds this remaining allocation.
    f = 2*np.pi*np.fft.rfftfreq(n, dt)
    keep = f <= band
    vf = rfft(values, axis=0)
    spectrum = np.zeros((len(f), H.shape[1]), dtype=complex)
    for point in range(H.shape[1]):
        for mode in range(6):
            response = (np.interp(f[keep], w, H[:, point, mode].real)
                        +1j*np.interp(f[keep], w, H[:, point, mode].imag))
            spectrum[keep, point] += np.conj(response)*vf[keep, mode]
    field = irfft(spectrum, n=n, axis=0)
    pre = field[:zero]
    result = field[zero+selected].copy()
    evidence = dict(method='intact full-complex offline FFT(v)',
        period_s=float(period), fft_dt_s=float(dt), n_fft=n, fft_time_origin_s=float(-zero*dt),
        maximum_omega_rad_s=float(band), hard_zero_above_band=True,
        transfer_interpolation='linear real and imaginary, supplied zero-frequency endpoint',
        saved_motion_interpolation=('none; exact native samples' if native else 'linear with every original native velocity knot exact'),
        prehistory=('none; unpadded periodic transform' if native else 'zero native velocity'), future_duration_s=float(future),
        future_policy=('none' if native else 'C1 Hermite from saved final velocity/acceleration to zero velocity/acceleration'),
        pre_response_max_abs_m=(None if native else np.max(np.abs(pre), axis=0).tolist()),
        temporal_policy=('native unpadded' if native else 'explicit historical padded'),
        raw_initial_radiation_m=field[zero].tolist(), initial_sample_modified=False,
        direct_term_policy='intact H, no separately added acceleration; raw initial value preserved',
        universal_causality_claim=False, numerical_accuracy_qualified=False,
        unresolved_response_above_band=True)
    return result, evidence


@dataclass
class SurfaceKernel:
    lag_s: np.ndarray
    values: np.ndarray
    direct_velocity: np.ndarray
    accepted: bool = False
    direct_acceleration: np.ndarray | None = None


# Availability gate only. The historical padded 60 s qualification does not
# qualify the current native/unpadded policy or exact mean-waterline trace.
_SURFACE_POLICY_ACCEPTED = True


def require_surface_readiness():
    """Check full-complex reconstruction capability, not numerical accuracy."""
    if not _SURFACE_POLICY_ACCEPTED:
        raise ValueError('WEC-Sim full-complex offline surface reconstruction is unavailable')


def cosine_kernel(omega, velocity_transfer, lag_s, *, direct_velocity, direct_acceleration=None):
    """Candidate regular causal kernel, never an automatic acceptance gate.

    An explicitly identified instantaneous velocity term is subtracted before
    cosine inversion.  Other pole/derivative terms must be removed analytically
    by the caller; this function must not be used to hide them with finite-band
    integration.  Supplied frequency endpoints/tails remain caller-reviewed.
    The returned candidate cannot publish a total field until accepted.
    """
    omega = _grid(omega, 'omega', zero=True)
    transfer = _finite(velocity_transfer, (len(omega), None, 6), 'velocity_transfer', complex_value=True)
    direct = _finite(direct_velocity, transfer.shape[1:], 'direct_velocity')
    acceleration = (np.zeros_like(direct) if direct_acceleration is None else
                    _finite(direct_acceleration, transfer.shape[1:], 'direct_acceleration'))
    lag = _grid(lag_s, 'lag_s', zero=True)
    if lag[0] != 0:
        raise ValueError('lag_s must start at zero')
    regular = (transfer-direct[None]+1j*omega[:, None, None]*acceleration[None]).real
    delta = np.diff(omega)
    weights = np.r_[delta[0]/2, (delta[:-1]+delta[1:])/2, delta[-1]/2]
    weighted = (weights[:, None, None]*regular).reshape(len(omega), -1)
    values = np.empty((len(lag), *transfer.shape[1:]))
    for start in range(0, len(lag), 64):
        cosine = np.cos(lag[start:start+64, None]*omega)
        values[start:start+64] = (2/np.pi*(cosine @ weighted)).reshape(-1, *transfer.shape[1:])
    return SurfaceKernel(lag, values, direct, direct_acceleration=acceleration)


def recover_transfer(kernel, omega):
    """Full complex response of a finite candidate kernel at held-out omega."""
    omega = _finite(omega, (None,), 'omega')
    output = np.empty((len(omega), *kernel.values.shape[1:]), dtype=complex)
    for i, w in enumerate(omega):
        output[i] = kernel.direct_velocity + trapezoid(
            kernel.values*np.exp(1j*w*kernel.lag_s)[:, None, None], x=kernel.lag_s, axis=0)
        if kernel.direct_acceleration is not None:
            output[i] -= 1j*w*kernel.direct_acceleration
    return output


def convolve_radiation(kernel, time_s, velocity, acceleration=None):
    """Zero-prehistory radiation from actual fixed-mean-CG modal velocities.

    These are the six native radiation-input velocities, not Euler angle rates.
    Lag coverage must extend through the whole requested trajectory, regardless
    of the different memory cutoff used by WEC-Sim's body-force operator.
    """
    if not kernel.accepted:
        raise ValueError('Transient surface kernel has not passed endpoint/tail and full-complex-transfer review')
    times = _grid(time_s, 'time_s', zero=True)
    v = _finite(velocity, (len(times), 6), 'velocity')
    dt = float(times[1]-times[0])
    if times[0] != 0 or not np.allclose(np.diff(times), dt, rtol=0, atol=1e-12):
        raise ValueError('Radiation history requires a uniform grid starting at zero')
    if np.any(v[0] != 0):
        raise ValueError('This reconstruction requires exact zero initial modal velocity')
    if len(kernel.lag_s) < len(times) or not np.allclose(kernel.lag_s[:len(times)], times, rtol=0, atol=1e-12):
        raise ValueError('Surface kernel must cover the entire requested history at its saved cadence')
    result = v @ kernel.direct_velocity.T
    if kernel.direct_acceleration is not None and np.any(kernel.direct_acceleration != 0):
        if acceleration is None:
            raise ValueError('Actual modal accelerations are required for the identified instantaneous acceleration term')
        a = _finite(acceleration, v.shape, 'acceleration')
        result += a @ kernel.direct_acceleration.T
    h = kernel.values[:len(times)]
    n = next_fast_len(2*len(times)-1)
    product = np.sum(rfft(h, n=n, axis=0)*rfft(v, n=n, axis=0)[:, None, :], axis=2)
    full = irfft(product, n=n, axis=0)[:len(times)]
    # Endpoint half weights implement integral_0^t at each t. Both endpoints
    # coincide at t=0, where the integral is exactly zero.
    result += dt*(full-.5*(v @ h[0].T)-.5*np.sum(h*v[0][None, None, :], axis=2))
    return result


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


def reconstruct_total_surface(hydro, setup, motion, surface_time_s, *, radiation_settings=None):
    """Return packed original total-field grids; hull interiors are invalid water.

    Numerical policy: fixed mean-hull linear diffraction/radiation, positive
    sideband frequencies, intact finite-band offline complex response over the complete
    trajectory, exact incident at targets, explicit exterior-trace continuation
    of disturbances near/inside the mean footprint. This is an export of a
    selected linear reference approximation, not nonlinear contact dynamics.
    No body pose, velocity, load or hydrodynamic coefficient is modified.
    """
    require_surface_readiness()
    from scipy.spatial.transform import Rotation, Slerp
    from threadpoolctl import threadpool_limits
    from .visualization import RenderSettings

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
    radiation_settings = OfflineRadiationSettings() if radiation_settings is None else radiation_settings
    position = _finite(motion.position_m, (len(stamp), 3), 'world position')
    rotation = Rotation.from_quat(motion.orientation_wxyz, scalar_first=True)
    matrices = rotation.as_matrix()
    poses = Slerp(stamp, rotation)(surface_time).as_matrix()
    positions = np.column_stack([np.interp(surface_time, stamp, position[:, j]) for j in range(3)])
    body = setup.body
    local_hull = body.mesh.vertices-body.cg
    hull_world = np.einsum('tij,vj->tvi', matrices, local_hull)+position[:, None, :]
    embedding = _finite(setup.world_from_native, (3, 3), 'world_from_native')
    hull_native = hull_world @ embedding
    beam = float(np.ptp(body.mesh.vertices[:, 1]))
    if beam <= 0:
        raise ValueError('A resolved boat beam is required for surface export')
    lo, hi = hull_native[:, :, :2].min(axis=(0, 1))-beam, hull_native[:, :, :2].max(axis=(0, 1))+beam
    # Reproduce the renderer's fixed WORLD square at its actual sampled poses,
    # then transform its corners. A native-axis square alone is insufficient
    # at an oblique global heading or a Slerp orientation extremum.
    physical_ids = np.unique(setup.physical_faces)
    sampled = np.einsum('tij,vj->tvi', poses, local_hull[physical_ids])+positions[:, None, :]
    view_lo, view_hi = sampled[:, :, :2].min(axis=(0, 1)), sampled[:, :, :2].max(axis=(0, 1))
    view_center = (view_lo+view_hi)/2
    view_half = float(np.max(view_hi-view_lo)/2*RenderSettings().framing)
    corners = np.array([[view_center[0]+sx*view_half, view_center[1]+sy*view_half, 0.]
                        for sx in (-1, 1) for sy in (-1, 1)]) @ embedding
    lo = np.minimum(lo, corners[:, :2].min(axis=0))-.1875
    hi = np.maximum(hi, corners[:, :2].max(axis=0))+.1875
    center = (lo+hi)/2
    half = float(np.max(hi-lo)/2)
    xy, faces, axes = _regular_grid([center-half, center+half], .1875)
    mapped, continued = exterior_trace_points(hydro, xy, 0., allow_boundary=True)
    field_xy, inverse = np.unique(mapped, axis=0, return_inverse=True)
    print(f'Surface: evaluating total linear field at {len(field_xy)} exterior/boundary points; '
          f'{int(continued.sum())} target points use the exact mean-waterline disturbance trace.', file=sys.stderr)
    print('Surface numerical approximation: intact complex offline response; '
          f'temporal policy={"native unpadded" if radiation_settings.period_s is None else "explicit historical padded"}; '
          'finite frequency band. Raw startup is preserved. '
          'Whole-domain numerical accuracy and remaining-tail uncertainty require separate qualification.', file=sys.stderr)
    maximum = radiation_settings.maximum_omega_rad_s
    step = radiation_settings.omega_step_rad_s
    if not np.isfinite([maximum, step]).all() or min(maximum, step) <= 0:
        raise ValueError('Radiation frequency range and spacing must be positive and finite')
    if not isinstance(radiation_settings.chunk_points, int) or radiation_settings.chunk_points < 1:
        raise ValueError('Radiation spatial chunk size must be a positive integer')
    count = int(np.ceil(maximum/step))
    omega = np.linspace(0., maximum, count+1)
    radiation_evidence = []
    with threadpool_limits(1):
        progress = _StageProgress('offline radiation reconstruction', len(field_xy), 'field points')
        radiation = np.empty((len(surface_time), len(field_xy)))
        for start in range(0, len(field_xy), radiation_settings.chunk_points):
            end = min(start+radiation_settings.chunk_points, len(field_xy))
            # radiation_results retains the source bank across spatial chunks;
            # only native field evaluation repeats, never the body/source solve.
            transfer = radiation_velocity_transfer(hydro, field_xy[start:end], omega[1:], shared_green=True, allow_boundary=True)
            values = np.concatenate([np.zeros_like(transfer[:1]), transfer])
            response, evidence = offline_radiation_response(omega, values, stamp, velocity, acceleration,
                output_time_s=surface_time, settings=radiation_settings)
            radiation[:, start:end] = response
            evidence.update(field_point_start=start, field_point_end=end)
            radiation_evidence.append(evidence)
            print('SURFACE_NUMERICS '+json.dumps(evidence, allow_nan=False, separators=(',', ':')), file=sys.stderr)
            del transfer, values
            progress.update(end)
        progress.finish()
        _, scattered_spectrum = combined_diffraction(
            hydro, setup.train, field_xy, return_forces=False, return_incident=False, allow_boundary=True)
        progress = _StageProgress('field synthesis', 1, 'completed field')
        scattered = synthesize_spectrum(setup.train.omega, scattered_spectrum, surface_time)
        del scattered_spectrum
        incident = synthesize_spectrum(setup.train.omega, incident_spectrum(setup.train, xy), surface_time)
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
    directory = Path(hydro.directory).parent/'surface-recovery'
    cache = _save_completed_cache(directory, completed)
    return _export_completed_field(completed, directory, cache)


def _export_completed_field(completed, recovery_directory, cache_path, *, timings=None, analyze=False):
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


def _analyze_completed_field_impl(completed, recovery_directory, cache_path, *, timings=None):
    """Geometry and metrics only; all physical field samples are immutable."""
    from .surface_geometry import validate_opening_clearance, waterline_segments

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
    mean_rotation = _finite(completed['mean_rotation'], (3, 3), 'cached mean rotation')
    mean_cg = _finite(completed['mean_cg_m'], (3,), 'cached mean CG')
    initial_area = float(completed['initial_wetted_area_m2'])
    beam = float(completed['beam_m'])
    calm_limit = float(completed['calm_eta_limit_m'])
    world_vertical_offset = float(completed.get('world_vertical_offset_m', 0.))
    graph = WaterGraphPartition(xy, faces)
    packed_vertices, packed_faces = [], []
    vertex_offsets, face_offsets = [0], [0]
    moments = {name: [] for name in ('full', 'one_beam')}
    wet_first = wet_last = wet_worst_time = None
    wet_worst = 0.
    local_initial_peak = None
    progress = _StageProgress('surface geometry and diagnostics', len(surface_time), 'frames')
    for i, t in enumerate(surface_time):
        world_hull = local_hull @ poses[i].T+positions[i]
        native_hull = world_hull @ embedding
        native_points = np.c_[xy, total[i]]
        world_points = native_points @ embedding.T
        current_rotation = embedding.T @ poses[i]
        reference_normal = current_rotation @ mean_rotation.T @ np.array([0., 0., 1.])
        native_cg = positions[i] @ embedding
        reference_offset = mean_cg[2]-reference_normal @ native_cg
        evidence = dict(frame_index=np.array(i), time_s=np.array(t),
            water_vertices_world_m=world_points, water_vertices_native_m=native_points,
            water_faces=faces, hull_vertices_world_m=world_hull, hull_vertices_native_m=native_hull,
            hull_faces=full_faces, physical_faces=physical_faces, opening_rim=rim, artificial_faces=caps,
            reference_normal_native=reference_normal, reference_offset_native=np.array(reference_offset),
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
            stage = 'wet-set-area'
            wet_fraction, _ = measured(stage, graph.wet_set_area, native_hull, physical_faces,
                                      total[i], reference_normal, reference_offset, initial_area)
            if wet_fraction > .10:
                if wet_first is None:
                    wet_first = float(t)
                wet_last = float(t)
                if wet_fraction > wet_worst:
                    wet_worst, wet_worst_time = wet_fraction, float(t)
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
    if wet_first is not None:
        print(f'WEC-Sim applicability: represented wet-set symmetric-difference area exceeds 10% of '
              f'initial wetted area; first={wet_first:g}s, worst={100*wet_worst:.4g}% at '
              f'{wet_worst_time:g}s, last={wet_last:g}s. This is a fixed-mean linear-model limitation, '
              'not a nonlinear motion-error estimate.', file=sys.stderr)
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
