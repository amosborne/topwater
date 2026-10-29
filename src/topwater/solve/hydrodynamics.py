"""In-memory native coefficients; no reference-backend exchange artifacts."""
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from topwater.constants import GRAVITY, WATER_DENSITY


@dataclass
class Hydrodynamics:
    directory: Path
    radiation_omega: np.ndarray
    body: object
    solver: object
    M: np.ndarray
    K: np.ndarray
    cg: np.ndarray
    cb: np.ndarray
    volume: float
    density: float
    gravity: float
    depth: float
    radiation_results: dict = field(default_factory=dict)
    packet_response: object = None
    radiation_damping: np.ndarray = field(init=False)
    Ainf: np.ndarray = field(init=False)


def prepare_body(wet_triangles, cg, mass_matrix, workspace, *, hydrostatic_stiffness, cb, volume,
                 radiation_omega, excitation_omega, field_omega, density=WATER_DENSITY,
                 gravity=GRAVITY, depth=np.inf, green_cache_dir=None):
    import capytaine as cpt
    import xarray as xr
    from importlib.metadata import version
    from capytaine.green_functions.delhommeau import Delhommeau
    from topwater.solve.types import MODES
    from topwater.solve.validation import _finite, _grid
    from topwater.solve.progress import _StageProgress
    from topwater.solve.bem_checks import check_frequency_bank, PrecheckedBEMSolver
    if version('capytaine') != '3.0.0':
        raise ValueError('The reviewed coefficient producer requires Capytaine 3.0.0')
    triangles = _finite(wet_triangles, (None, 3, 3), 'wet_triangles')
    center = _finite(cg, (3,), 'cg')
    M = _finite(mass_matrix, (6, 6), 'mass_matrix')
    K = _finite(hydrostatic_stiffness, (6, 6), 'hydrostatic_stiffness')
    cb = _finite(cb, (3,), 'cb')
    if not len(triangles) or not np.isfinite([density, gravity, volume]).all() or min(density, gravity, volume) <= 0:
        raise ValueError('Positive density, gravity, displaced volume and wet triangles are required')
    if not (depth == np.inf or np.isfinite(depth) and depth > 0):
        raise ValueError('depth must be positive or infinite')
    rad_w = _grid(radiation_omega, 'radiation_omega')
    exc_w = _grid(excitation_omega, 'excitation_omega')
    field_w = _grid(field_omega, 'field_omega') if len(field_omega) else np.empty(0)
    directory = Path(workspace)
    directory.mkdir(parents=True, exist_ok=False)
    # Preserve the established triangle/collocation convention on both paths.
    starts = np.argmax(np.linalg.norm(triangles-triangles.mean(axis=1)[:, None], axis=2), axis=1)
    triangles = np.array([np.roll(t, -int(i), axis=0) for t, i in zip(triangles, starts)])
    vertices, ids = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
    mesh = cpt.Mesh(vertices=vertices, faces=ids.reshape(-1, 3), auto_clean=False, auto_check=False)
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=center),
        center_of_mass=center, mass=float(M[0, 0]), name='boat')
    labels = dict(influenced_dof=list(MODES), radiating_dof=list(MODES))
    body.inertia_matrix = xr.DataArray(M, dims=list(labels), coords=labels)
    body.hydrostatic_stiffness = xr.DataArray(K, dims=list(labels), coords=labels)
    checked = check_frequency_bank(body, np.r_[rad_w, exc_w, field_w, np.inf],
                                   density=density, gravity=gravity, depth=depth)
    cache = None if green_cache_dir is None else str(Path(green_cache_dir))
    if cache is not None:
        Path(cache).mkdir(parents=True, exist_ok=True)
    progress = _StageProgress('Green-function tabulation', 1, 'table')
    solver = PrecheckedBEMSolver(checked_frequencies=checked,
        green_function=Delhommeau(tabulation_cache_dir=cache))
    progress.finish()
    return Hydrodynamics(directory=directory, radiation_omega=rad_w, body=body, solver=solver,
        M=M, K=K, cg=center, cb=cb, volume=float(volume), density=float(density),
        gravity=float(gravity), depth=depth)


def prepare_radiation(hydro):
    import capytaine as cpt
    from topwater.solve.types import MODES
    from topwater.solve.progress import _StageProgress
    hydro.Ainf = np.empty((6, 6))
    hydro.radiation_damping = np.empty((len(hydro.radiation_omega), 6, 6))
    common = dict(body=hydro.body, rho=hydro.density, g=hydro.gravity, water_depth=hydro.depth)
    progress = _StageProgress('native radiation coefficients', 6*(len(hydro.radiation_omega)+1), 'BEM solves')
    count = 0
    with threadpool_limits(1):
        for omega in np.r_[np.inf, hydro.radiation_omega]:
            for j, mode in enumerate(MODES):
                result = hydro.solver.solve(cpt.RadiationProblem(**common,
                    omega=float(omega), radiating_dof=mode), keep_details=True)
                hydro.radiation_results[float(omega), mode] = result
                if np.isinf(omega):
                    hydro.Ainf[:, j] = [result.added_mass[m] for m in MODES]
                else:
                    i = np.searchsorted(hydro.radiation_omega, omega)
                    hydro.radiation_damping[i, :, j] = [result.radiation_damping[m] for m in MODES]
                count += 1
                progress.update(count)
    progress.finish()
    return hydro


def prepare_hydrodynamics(setup, workspace, options):
    from topwater.solve.fourier import surface_radiation_frequencies
    radiation = options.radiation_frequencies(setup.train.omega[-1])
    h = setup.hydrostatic
    hydro = prepare_body(setup.wet_triangles_m, setup.mean_cg_m, setup.mass_matrix, workspace,
        hydrostatic_stiffness=h.stiffness, cb=h.center_of_buoyancy, volume=h.displaced_volume,
        excitation_omega=setup.train.omega[1:], radiation_omega=radiation,
        field_omega=np.unique(np.r_[setup.train.omega[1:],
            surface_radiation_frequencies(setup.time_s, options.radiation_settings(1))[1:]]),
        green_cache_dir=Path.home()/'.cache/topwater/capytaine-3.0.0')
    return prepare_radiation(hydro)
