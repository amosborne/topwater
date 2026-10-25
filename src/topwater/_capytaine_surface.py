"""Native Capytaine coefficients and reviewed building blocks for surface export.

Capytaine is imported only when a producer is explicitly called.  All harmonic
quantities use Re(A exp(-i omega t)), fixed mean-world axes and the mean CG.
Coefficient production does not qualify transient field reconstruction.  A
surface kernel is deliberately unusable until its endpoint treatment and full
complex transfer recovery have been accepted by the caller's validation gate.
"""
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.signal import fftconvolve
from scipy.spatial import ConvexHull, QhullError


MODES = ('Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw')


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


def combined_diffraction(hydro, train, xy, *, return_forces=False, checkpoint=None):
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


def radiation_velocity_transfer(hydro, xy, omega=None, *, checkpoint=None, shared_green=False):
    """Return H_eta,v=i H_eta,q/omega, shape (omega, point, six modes)."""
    import capytaine as cpt
    checkpoint = (lambda: None) if checkpoint is None else checkpoint
    xy = exterior_field_points(hydro, xy)
    frequencies = hydro.radiation_omega if omega is None else _grid(omega, 'omega')
    values = np.empty((len(frequencies), len(xy), 6), dtype=complex)
    for iw, w in enumerate(frequencies):
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


def require_surface_readiness():
    """Fail before body execution until the field pilot has qualified tails.

    Having Capytaine field-evaluation functions is not evidence that the finite
    band, causal surface reconstruction is physically ready.  Replace this
    gate only after the independent endpoint/transfer pilot has been accepted.
    """
    raise ValueError('WEC-Sim total-surface reconstruction awaits the retained-field endpoint/causal-transfer pilot')


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
    for j in range(6):
        for p in range(kernel.values.shape[1]):
            h = kernel.values[:len(times), p, j]
            # Endpoint half weights implement integral_0^t at each t.  Both
            # endpoints coincide at t=0, where the integral is exactly zero.
            full = fftconvolve(v[:, j], h)[:len(times)]
            result[:, p] += dt*(full-.5*h[0]*v[:, j]-.5*h*v[0, j])
    return result
