"""Native Capytaine coefficients and reviewed building blocks for surface export.

Capytaine is imported only when a producer is explicitly called.  All harmonic
quantities use Re(A exp(-i omega t)), fixed mean-world axes and the mean CG.
Coefficient production does not qualify transient field reconstruction. A
surface kernel requires explicit acceptance of its numerical approximation;
the selected causal projection retains a documented near-field phase error.
"""
from dataclasses import dataclass
from functools import wraps
from importlib.metadata import version
import json
import logging
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import trapezoid
from scipy.fft import rfft, irfft, next_fast_len
from scipy.spatial import ConvexHull, QhullError


MODES = ('Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw')


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


def exterior_field_points(hydro, xy):
    """Admit only points strictly outside the mean-waterline footprint.

    Capytaine's linear field is exterior to the fixed mean hull.  Moving-hull
    clipping cannot make an interior field evaluation physical.  The waterline
    convex hull is a conservative footprint: concave waterlines may reject
    extra points, but an interior point is never admitted by that approximation.
    Exact/on-boundary traces and near-hull continuation require a separately
    reviewed exterior-offset evaluation; this function never invents one.
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
    if np.any(np.all(distances <= tolerance, axis=1)):
        raise ValueError('Field points must lie strictly outside the mean hull footprint; interior/on-waterline Green-function evaluation is prohibited')
    return points


def exterior_trace_points(hydro, xy, offset_m):
    """Map a mean-footprint interior strip to an explicit exterior trace.

    Only disturbance fields may use these mapped points. The incident field
    must still be evaluated at the original XY. This is a bounded geometric
    continuation for export, not an exact nonlinear contact condition.
    """
    points = _finite(xy, (None, 2), 'xy')
    offset = float(offset_m)
    if not np.isfinite(offset) or offset <= 0:
        raise ValueError('Exterior trace offset must be positive and finite')
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
    exterior_field_points(hydro, mapped)
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
    the separate body-force table has a periodic full-heading seam.

    This function performs numerical work and must only be called by the sole
    coordinator-allocated numerical owner.  Its files belong in scratch space.
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
    # A coordinator-owned job-temp cache can be shared by producer instances.
    # None preserves the original Green-function settings and no global cache.
    cache = None if green_cache_dir is None else str(Path(green_cache_dir))
    if cache is not None:
        Path(cache).mkdir(parents=True, exist_ok=True)
    solver = cpt.BEMSolver(green_function=Delhommeau(tabulation_cache_dir=cache))
    radiation, diffraction, imported = {}, {}, []
    A = np.empty((len(rad_w), 6, 6)); B = np.empty_like(A)
    Ainf = np.empty((6, 6))
    H = np.empty((len(exc_w), len(angles), 6), dtype=complex)
    common = dict(body=body, rho=float(density), g=float(gravity), water_depth=depth)
    import_angles = (0., np.pi/4, np.pi/2)
    with threadpool_limits(1):
        for j, mode in enumerate(MODES):
            checkpoint()
            result = solver.solve(cpt.RadiationProblem(**common, omega=np.inf, radiating_dof=mode), keep_details=True, n_threads=1)
            Ainf[:, j] = [result.added_mass[m] for m in MODES]
            radiation[(np.inf, mode)] = result
        all_w = np.unique(np.r_[rad_w, exc_w])
        for w in all_w:
            checkpoint()
            rad_index = np.flatnonzero(rad_w == w)
            exc_index = np.flatnonzero(exc_w == w)
            if len(rad_index):
                for j, mode in enumerate(MODES):
                    checkpoint()
                    result = solver.solve(cpt.RadiationProblem(**common, omega=float(w), radiating_dof=mode), keep_details=True, n_threads=1)
                    radiation[(float(w), mode)] = result
                    imported.append(result)
                    A[rad_index[0], :, j] = [result.added_mass[m] for m in MODES]
                    B[rad_index[0], :, j] = [result.radiation_damping[m] for m in MODES]
            query_angles = angles[:-1] if len(exc_index) else np.asarray(import_angles)
            for a in query_angles:
                checkpoint()
                result = solver.solve(cpt.DiffractionProblem(**common, omega=float(w), wave_direction=float(a)), keep_details=True, n_threads=1)
                diffraction[(float(w), float(a))] = result
                if len(rad_index) and any(abs(a-b) < 1e-12 for b in import_angles):
                    imported.append(result)
                if len(exc_index):
                    j = int(np.flatnonzero(angles[:-1] == a)[0])
                    fk = froude_krylov_force(result.problem)
                    H[exc_index[0], j] = [result.forces[m] + fk[m] for m in MODES]
            if len(exc_index):
                H[exc_index[0], -1] = H[exc_index[0], 0]
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
    savemat(directory/'native.mat', dict(Ainf=Ainf, M=M, K=K, cg=center, cb=cb,
             volume=float(volume), rho=float(density), g=float(gravity)))
    return Hydrodynamics(directory, directory/'canonical.nc', directory/'native.mat', rad_w, exc_w,
                         angles, A, B, Ainf, H, M, K, center, cb, float(volume), body,
                         solver, radiation, diffraction, float(density), float(gravity), depth)


@_summarize_capytaine_warnings
def combined_diffraction(hydro, train, xy, *, return_forces=False, checkpoint=None, return_incident=True):
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
    xy = exterior_field_points(hydro, xy)
    omega = _grid(train.omega, 'train.omega', zero=True)
    C = _finite(train.C, (None, 2, len(omega)), 'train.C', complex_value=True)
    theta = _finite(train.theta, (len(C),), 'train.theta')
    incident = np.zeros((len(omega), len(xy)), dtype=complex)
    scattered = np.zeros_like(incident)
    force = np.zeros((len(omega), 6), dtype=complex)
    if omega[0] != 0 or np.any(C[:, :, 0] != 0):
        raise ValueError('A zero-DC signed-kappa train is required')
    common = dict(body=hydro.body, rho=hydro.density, g=hydro.gravity, water_depth=hydro.depth)
    for iw in range(1, len(omega)):
        checkpoint()
        if len(omega) > 100 and (iw % 100 == 0 or iw == len(omega)-1):
            print(f'Surface diffraction: {iw}/{len(omega)-1} positive sidebands', file=sys.stderr)
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
        result = hydro.solver.solve(problem, keep_details=True, n_threads=1)
        scattered[iw] = hydro.solver.compute_free_surface_elevation(xy, result)
        if return_forces:
            force[iw] = fk + np.array([result.forces[m] for m in MODES])
    return (incident, scattered, force) if return_forces else (incident, scattered)


def synthesize_spectrum(omega, coefficients, time_s):
    """Synthesize already quadrature-weighted positive-frequency coefficients."""
    omega = _grid(omega, 'omega', zero=True)
    values = _finite(coefficients, (len(omega), None), 'coefficients', complex_value=True)
    times = _finite(time_s, (None,), 'time_s')
    return (np.exp(-1j*times[:, None]*omega) @ values).real


@_summarize_capytaine_warnings
def radiation_velocity_transfer(hydro, xy, omega=None, *, checkpoint=None, shared_green=False):
    """Return H_eta,v=i H_eta,q/omega, shape (omega, point, six modes)."""
    import capytaine as cpt
    checkpoint = (lambda: None) if checkpoint is None else checkpoint
    xy = exterior_field_points(hydro, xy)
    frequencies = hydro.radiation_omega if omega is None else _grid(omega, 'omega')
    values = np.empty((len(frequencies), len(xy), 6), dtype=complex)
    for iw, w in enumerate(frequencies):
        if len(frequencies) > 100 and (iw % 100 == 0 or iw == len(frequencies)-1):
            print(f'Surface radiation: {iw+1}/{len(frequencies)} frequencies', file=sys.stderr)
        results = []
        for j, mode in enumerate(MODES):
            checkpoint()
            key = (float(w), mode)
            if key not in hydro.radiation_results:
                hydro.radiation_results[key] = hydro.solver.solve(cpt.RadiationProblem(
                    body=hydro.body, rho=hydro.density, g=hydro.gravity,
                    water_depth=hydro.depth, omega=float(w), radiating_dof=mode), keep_details=True, n_threads=1)
            results.append(hydro.radiation_results[key])
            if not shared_green:
                eta_q = hydro.solver.compute_free_surface_elevation(xy, hydro.radiation_results[key])
                values[iw, :, j] = 1j*eta_q/float(w)
        if shared_green:
            # This is precisely the public field evaluator's S @ sources,
            # batching the six source columns. No new BEM equation or solver.
            representative = results[0]
            S = hydro.solver.engine.build_S_matrix(np.c_[xy, np.zeros(len(xy))],
                hydro.body.mesh_including_lid, free_surface=representative.free_surface,
                water_depth=representative.water_depth, wavenumber=representative.encounter_wavenumber)
            sources = np.column_stack([r.sources for r in results])
            # i/omega * (i*omega/g) * potential = -potential/g.
            values[iw] = -(S @ sources)/hydro.gravity
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


@dataclass
class SurfaceKernel:
    lag_s: np.ndarray
    values: np.ndarray
    direct_velocity: np.ndarray
    accepted: bool = False
    direct_acceleration: np.ndarray | None = None


_SURFACE_POLICY_ACCEPTED = True


def require_surface_readiness():
    """Check availability of the coordinator-accepted provisional field policy.

    Acceptance is selected-case engineering qualification, not full-complex
    recovery or nonlinear wave/body accuracy. Runtime geometry/startup gates
    still apply independently to every produced result.
    """
    if not _SURFACE_POLICY_ACCEPTED:
        raise ValueError('WEC-Sim total-surface reconstruction awaits selected-case causal-projection pilot acceptance')


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
    result = []
    previous = points[-1]
    dp = float(previous @ normal+offset)
    for point in points:
        dc = float(point @ normal+offset)
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
    return float(np.linalg.norm(np.cross(relative, np.roll(relative, -1, axis=0)).sum(axis=0))/2)


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
        owners = {}
        for i, face in enumerate(self.faces):
            for a, b in zip(face, np.roll(face, -1)):
                key = tuple(sorted((int(a), int(b))))
                owners[key] = min(i, owners.get(key, i))
        self.owner = np.array([[owners[tuple(sorted((int(a), int(b))))]
                               for a, b in zip(face, np.roll(face, -1))] for face in self.faces])
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
                if len(polygon) < 3 or _polygon_area(polygon) <= eps**2:
                    continue
                wet = _clip_area_polygon(polygon, water_normal[i], water_offset[i], eps)
                mean = _clip_area_polygon(polygon, normal, reference_offset, eps)
                common = _clip_area_polygon(wet, normal, reference_offset, eps)
                current_area += _polygon_area(wet)
                mean_area += _polygon_area(mean)
                common_area += _polygon_area(common)
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


def _total_waterline(vertices, faces, hull, domain_bounds):
    """Actual clipped total-surface contact segments, excluding outer patch.

    Incidence-one edges on the known rectangular patch boundary are external.
    Every remaining boundary edge must lie on a support plane of the posed
    hull, with both endpoints inside all its halfspaces. No horizontal hull
    outline is substituted for a missing or invalid physical contact line.
    """
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    unique, count = np.unique(edges, axis=0, return_counts=True)
    boundary = vertices[unique[count == 1]]
    if not len(boundary):
        raise ValueError('Total surface has no resolved waterline boundary')
    scale = max(1., float(np.ptp(hull, axis=0).max()))
    tolerance = 1e-8*scale
    bounds = np.asarray(domain_bounds, dtype=float)
    if bounds.shape != (2, 2):
        raise ValueError('Waterline classification requires the known rectangular patch bounds')
    outer = np.zeros(len(boundary), dtype=bool)
    for axis in range(2):
        for value in bounds[:, axis]:
            outer |= np.all(abs(boundary[:, :, axis]-value) <= tolerance, axis=1)
    segments = boundary[~outer]
    if not len(segments):
        raise ValueError('Total surface has no resolved hull-contact waterline; near-band metric is undefined')
    planes = ConvexHull(hull).equations
    distances = np.einsum('evi,fi->evf', segments, planes[:, :3])+planes[:, 3]
    on_support = np.any(np.all(abs(distances) <= tolerance, axis=1), axis=1)
    inside = np.all(distances <= tolerance, axis=(1, 2))
    if not np.all(on_support & inside):
        raise ValueError('Surface inner boundary does not lie on the actual posed hull waterline')
    return segments


def _waterline_distance_xy(points, segments):
    a, b = segments[:, 0, :2], segments[:, 1, :2]
    edge = b-a
    length_squared = np.sum(edge*edge, axis=1)
    # Vertical waterline segments have a point projection; their horizontal
    # distance is still defined, without division by zero or discarding them.
    relative = points[:, None]-a[None]
    fraction = np.clip(np.sum(relative*edge[None], axis=2)/np.where(length_squared > 0, length_squared, 1.), 0, 1)
    distance = np.linalg.norm(points[:, None]-(a[None]+fraction[:, :, None]*edge[None]), axis=2).min(axis=1)
    return distance


def _delta_moments(vertices, faces, incident_at_vertices, hull, beam, domain_bounds):
    """Exact P1 delta moments; centroid classification for the one-beam band."""
    triangle = vertices[faces]
    edge1, edge2 = triangle[:, 1]-triangle[:, 0], triangle[:, 2]-triangle[:, 0]
    area = abs(edge1[:, 0]*edge2[:, 1]-edge1[:, 1]*edge2[:, 0])/2
    delta = vertices[:, 2]-incident_at_vertices
    q = delta[faces]
    first = area*q.mean(axis=1)
    square = area/6*(np.sum(q*q, axis=1)+q[:, 0]*q[:, 1]+q[:, 1]*q[:, 2]+q[:, 2]*q[:, 0])
    segments = _total_waterline(vertices, faces, hull, domain_bounds)
    near = _waterline_distance_xy(triangle[:, :, :2].mean(axis=1), segments) <= beam
    output = {}
    for name, mask in [('full', np.ones(len(faces), dtype=bool)), ('one_beam', near)]:
        if np.any(mask):
            subset = q[mask]
            it, iv = np.unravel_index(np.argmax(abs(subset)), subset.shape)
            peak = vertices[faces[mask][it, iv]]
            maximum = float(abs(subset[it, iv]))
        else:
            peak, maximum = np.zeros(3), 0.
        inc = incident_at_vertices[faces]
        incident_square = area/6*(np.sum(inc*inc, axis=1)+inc[:, 0]*inc[:, 1]
                                   +inc[:, 1]*inc[:, 2]+inc[:, 2]*inc[:, 0])
        output[name] = np.r_[area[mask].sum(), first[mask].sum(), square[mask].sum(),
                             maximum, incident_square[mask].sum(), peak]
    return output


def reconstruct_total_surface(hydro, setup, motion, surface_time_s):
    """Return packed total water-air meshes, clipped to each actual posed hull.

    Numerical policy: fixed mean-hull linear diffraction/radiation, positive
    sideband frequencies, finite-band causal field operator over the complete
    trajectory, exact incident at targets, explicit exterior-trace continuation
    of disturbances near/inside the mean footprint. This is an export of a
    selected linear reference approximation, not nonlinear contact dynamics.
    No body pose, velocity, load or hydrodynamic coefficient is modified.
    """
    require_surface_readiness()
    from scipy.spatial.transform import Rotation, Slerp
    from threadpoolctl import threadpool_limits
    from ._surface_export import export_surface
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
    trace_offset = min(.025, .01*beam)
    mapped, continued = exterior_trace_points(hydro, xy, trace_offset)
    field_xy, inverse = np.unique(mapped, axis=0, return_inverse=True)
    print(f'Surface: evaluating total linear field at {len(field_xy)} exterior points; '
          f'{int(continued.sum())} target points use an explicit {trace_offset:g} m disturbance trace.', file=sys.stderr)
    print('Surface numerical approximation: the selected causal finite-band projection showed '
          '42.7 mm maximum / 12.5 mm RMS difference from a signed-lag diagnostic at tested near/far points. '
          'That noncausal diagnostic is not a reference solution or an accuracy bound; '
          'full-complex near-field transfer recovery has not passed.', file=sys.stderr)
    # Delta omega < pi/T keeps the transform's periodic images outside the
    # requested +/-T window. It does not by itself certify causal projection.
    maximum = 24.
    step = min(.025, .98*np.pi/stamp[-1])
    count = int(np.ceil(maximum/step))
    omega = np.linspace(0., maximum, count+1)
    with threadpool_limits(1):
        transfer = radiation_velocity_transfer(hydro, field_xy, omega[1:], shared_green=True)
        radiation = np.empty((len(surface_time), len(field_xy)))
        for start in range(0, len(field_xy), 128):
            end = min(start+128, len(field_xy))
            values = np.concatenate([np.zeros_like(transfer[:1, start:end]), transfer[:, start:end]])
            candidate = cosine_kernel(omega, values, stamp,
                                      direct_velocity=np.zeros((end-start, 6)))
            # This is the fixed selected finite-band numerical policy. Its
            # actual-response sensitivity and causal-projection bias are stated
            # in the producer evidence; no full-complex-recovery pass is implied.
            candidate.accepted = True
            response = convolve_radiation(candidate, stamp, velocity)
            for j in range(end-start):
                radiation[:, start+j] = np.interp(surface_time, stamp, response[:, j])
        del transfer
        _, scattered_spectrum, force_spectrum = combined_diffraction(
            hydro, setup.train, field_xy, return_forces=True, return_incident=False)
        scattered = synthesize_spectrum(setup.train.omega, scattered_spectrum, surface_time)
        del scattered_spectrum
        incident = synthesize_spectrum(setup.train.omega, incident_spectrum(setup.train, xy), surface_time)
        total = incident+scattered[:, inverse]+radiation[:, inverse]
        direct_force = synthesize_spectrum(setup.train.omega, force_spectrum, stamp)
        prepared_force, _ = setup.train.excitation(hydro.excitation_omega, hydro.excitation_headings,
                                                 hydro.excitation, hydro.K, stamp)
    difference = np.max(abs(direct_force-prepared_force), axis=0)
    print('Surface: exact-sideband field/prepared-body load maximum differences '
          f'[N,N,N,Nm,Nm,Nm]={difference.tolist()}; the accepted body loads are unchanged.', file=sys.stderr)

    graph = WaterGraphPartition(xy, faces)
    mean_rotation = Rotation.from_quat(setup.equilibrium.orientation, scalar_first=True).as_matrix()
    packed_vertices, packed_faces = [], []
    vertex_offsets, face_offsets = [0], [0]
    moments = {name: [] for name in ('full', 'one_beam')}
    wet_first = wet_last = wet_worst_time = None
    wet_worst = 0.
    local_initial_peak = None
    for i, t in enumerate(surface_time):
        if len(surface_time) > 100 and (i % 100 == 0 or i == len(surface_time)-1):
            print(f'Surface export: {i+1}/{len(surface_time)} frames', file=sys.stderr)
        world_hull = local_hull @ poses[i].T+positions[i]
        native_hull = world_hull @ embedding
        native_points = np.c_[xy, total[i]]
        world_points = native_points @ embedding.T
        vertices, connectivity = export_surface(world_points, faces, world_hull, body.mesh.faces,
            opening_rim=body.mesh.opening_rim, artificial_faces=body.mesh.artificial_faces)
        # Metrics use the same represented incident P1 graph, interpolated at
        # any new cut vertices. They measure total-minus-incident, not a second
        # direct incident evaluation inconsistent with the represented graph.
        clipped_native = vertices @ embedding
        incident_cut = _grid_p1_values(axes, incident[i], clipped_native[:, :2])
        frame_moments = _delta_moments(clipped_native, connectivity, incident_cut, native_hull, beam,
                                      [[axes[0][0], axes[1][0]], [axes[0][-1], axes[1][-1]]])
        for name in moments:
            moments[name].append(frame_moments[name])
        current_rotation = embedding.T @ poses[i]
        reference_normal = current_rotation @ mean_rotation.T @ np.array([0., 0., 1.])
        native_cg = positions[i] @ embedding
        reference_offset = setup.mean_cg_m[2]-reference_normal @ native_cg
        wet_fraction, _ = graph.wet_set_area(native_hull, setup.physical_faces, total[i],
            reference_normal, reference_offset, setup.hydrostatic.wetted_area)
        if wet_fraction > .10:
            if wet_first is None:
                wet_first = float(t)
            wet_last = float(t)
            if wet_fraction > wet_worst:
                wet_worst, wet_worst_time = wet_fraction, float(t)
        if i == 0:
            # Exact extrema of the represented local P1 triangles, conservatively
            # selecting every triangle touching the initial hull XY bounding box.
            triangle = clipped_native[connectivity]
            hmin, hmax = native_hull[:, :2].min(axis=0), native_hull[:, :2].max(axis=0)
            region = np.all(triangle[:, :, :2].max(axis=1) >= hmin, axis=1) & np.all(triangle[:, :, :2].min(axis=1) <= hmax, axis=1)
            if not np.any(region):
                raise ValueError('Initial total surface lacks boat-region water coverage')
            local_initial_peak = float(np.max(abs(triangle[region, :, 2])))
            limit = .001*float(np.sum(setup.train.carriers.amplitudes))
            if local_initial_peak > limit+1e-12:
                raise ValueError(f'Initial represented total surface exceeds local calm eta limit: {local_initial_peak:g} > {limit:g} m')
        packed_vertices.append(vertices); packed_faces.append(connectivity)
        vertex_offsets.append(vertex_offsets[-1]+len(vertices))
        face_offsets.append(face_offsets[-1]+len(connectivity))
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
          'triangle-centroid horizontal distance to actual total-surface/hull waterline): '+json.dumps(summary), file=sys.stderr)
    if wet_first is not None:
        print(f'WEC-Sim applicability: represented wet-set symmetric-difference area exceeds 10% of '
              f'initial wetted area; first={wet_first:g}s, worst={100*wet_worst:.4g}% at '
              f'{wet_worst_time:g}s, last={wet_last:g}s. This is a fixed-mean linear-model limitation, '
              'not a nonlinear motion-error estimate.', file=sys.stderr)
    print(f'Surface startup: represented total local |eta|<={local_initial_peak:g}m. '
          'The separate full-hull incident eta/velocity bounds remain in force; '
          'no total disturbed-fluid velocity bound has been established.', file=sys.stderr)
    return (np.concatenate(packed_vertices), np.concatenate(packed_faces),
            np.asarray(vertex_offsets, dtype=np.int64), np.asarray(face_offsets, dtype=np.int64))
