"""Validation-only Capytaine reference on exactly supplied wetted triangles."""
from pathlib import Path
from contextlib import contextmanager
from copy import deepcopy
from hashlib import sha256
from numbers import Integral
import platform
import numpy as np
from numpy.testing import assert_allclose
from threadpoolctl import threadpool_limits


MODES = ('Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw')


@contextmanager
def _prony_realization(green, seed):
    """Isolate a sequential reference call, without altering Capytaine's fit.

    Capytaine 3.0 uses a module-owned default_rng, not np.random's legacy RNG.
    This temporarily replaces that RNG; concurrent Capytaine calls in the same
    process are unsupported. No seed leaves the legacy path entirely unchanged.
    """
    if seed is None:
        yield None
        return
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError('Prony seed must be a nonnegative integer or None')
    import capytaine as cpt
    import scipy
    from capytaine.tools import prony_decomposition as prony
    from capytaine.green_functions import delhommeau

    sources = dict(prony=Path(prony.__file__), delhommeau=Path(delhommeau.__file__),
                   reference_helper=Path(__file__))
    core_path = getattr(getattr(green, 'fortran_core', None), '__file__', None)
    if core_path is not None:
        sources['compiled_green_function'] = Path(core_path)
    generator = np.random.default_rng(int(seed))
    metadata = dict(seed=int(seed), initial_rng_state=deepcopy(generator.bit_generator.state),
        fitted_decompositions=[], runtime=dict(python=platform.python_version(), capytaine=cpt.__version__, numpy=np.__version__,
            scipy=scipy.__version__),
        source_sha256={name: sha256(path.read_bytes()).hexdigest() for name, path in sources.items()},
        fit_policy='Unmodified Capytaine backend and tolerance; actual returned rows are lambda,a. '
            'Seed fixes this sequential call, not historical reference provenance.')
    original_rng = prony.RNG
    original_method = green.find_best_exponential_decomposition
    had_override = 'find_best_exponential_decomposition' in green.__dict__
    original_override = green.__dict__.get('find_best_exponential_decomposition')
    seen = set()

    def capture(kh, *, method=None):
        values = original_method(kh, method=method)
        key = (float(kh), method)
        if key not in seen:
            coefficients = np.asarray(values)
            metadata['fitted_decompositions'].append(dict(
                dimensionless_wavenumber='infinity' if np.isinf(kh) else float(kh),
                method=green.finite_depth_prony_decomposition_method if method is None else method,
                lambda_and_amplitude=coefficients.tolist(),
                coefficients_sha256=sha256(np.asarray(coefficients, dtype='<f8').tobytes()).hexdigest()))
            seen.add(key)
        return values

    try:
        prony.RNG = generator
        green.find_best_exponential_decomposition = capture
        yield metadata
    finally:
        metadata['final_rng_state'] = deepcopy(generator.bit_generator.state)
        prony.RNG = original_rng
        if had_override:
            green.find_best_exponential_decomposition = original_override
        else:
            green.__dict__.pop('find_best_exponential_decomposition', None)


def body_from_triangles(triangles, center, mass):
    """Exact supplied panels and world CG, without cleaning or mesh healing."""
    import capytaine as cpt

    assert cpt.__version__ == '3.0.0'
    vertices, indices = np.unique(triangles.reshape(-1, 3), axis=0,
                                  return_inverse=True)
    faces = indices.reshape(-1, 3)
    assert_allclose(vertices[faces], triangles, atol=0, rtol=0)
    mesh = cpt.Mesh(vertices=vertices, faces=faces, auto_clean=False, auto_check=False)
    normal = np.cross(triangles[:, 1]-triangles[:, 0],
                      triangles[:, 2]-triangles[:, 0])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    assert_allclose(mesh.faces_normals, normal, atol=1e-13, rtol=1e-12)
    return cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=center),
                            mass=mass, center_of_mass=center)


def instantaneous_inertia(triangles, center, mass):
    import capytaine as cpt

    body = body_from_triangles(triangles, center, mass)
    # Discover native BLAS/OpenMP only after loading the Green-function library.
    green = cpt.Delhommeau(tabulation_nr=0, tabulation_nz=0, tabulation_cache_dir=None)
    with threadpool_limits(1):
        solver = cpt.BEMSolver(green_function=green)
        reference, damping = np.empty((6, 6)), np.empty((6, 6))
        for column, name in enumerate(MODES):
            result = solver.solve(cpt.RadiationProblem(body=body, radiating_dof=name,
                omega=np.inf, rho=1000., g=9.81, water_depth=np.inf))
            reference[:, column] = [result.added_mass[row] for row in MODES]
            damping[:, column] = [result.radiation_damping[row] for row in MODES]
    return reference, damping


def frequency_response(triangles, center, mass, frequencies, headings, *, water_depth=np.inf,
                       prony_seed=None):
    """Finite-frequency reference only; excitation includes incident pressure.

    Unit wave elevation at the world origin, Re(amplitude*exp(-i*omega*t)),
    specified water depth (infinite by default), zero forward speed. Forces/moments use world axes about CG.
    Optional prony_seed isolates and records the actual finite-depth approximation;
    None retains the existing unseeded behavior. Seeded calls must be sequential.
    """
    import capytaine as cpt
    from capytaine.bem.airy_waves import froude_krylov_force

    body = body_from_triangles(triangles, center, mass)
    cache=Path('artifacts/capytaine-green-cache')
    cache.mkdir(parents=True,exist_ok=True)
    green = cpt.Delhommeau(tabulation_cache_dir=str(cache))
    added_mass = np.empty((len(frequencies), 6, 6))
    damping = np.empty_like(added_mass)
    diffraction = np.empty((len(frequencies), len(headings), 6), dtype=complex)
    incident = np.empty_like(diffraction)
    with _prony_realization(green, prony_seed) as realization, threadpool_limits(1):
        solver = cpt.BEMSolver(green_function=green)
        for f, omega in enumerate(frequencies):
            common = dict(body=body, omega=omega, rho=1000., g=9.81, water_depth=water_depth)
            for column, name in enumerate(MODES):
                result = solver.solve(cpt.RadiationProblem(**common, radiating_dof=name))
                added_mass[f, :, column] = [result.added_mass[row] for row in MODES]
                damping[f, :, column] = [result.radiation_damping[row] for row in MODES]
            for h, heading in enumerate(headings):
                problem = cpt.DiffractionProblem(**common, wave_direction=heading)
                result = solver.solve(problem)
                fk = froude_krylov_force(problem)
                diffraction[f, h] = [result.forces[row] for row in MODES]
                incident[f, h] = [fk[row] for row in MODES]
    output = dict(added_mass=added_mass, radiation_damping=damping,
                green_function_settings=green.exportable_settings,
                diffraction=diffraction, froude_krylov=incident,
                excitation=diffraction+incident)
    if realization is not None:
        output['prony_realization'] = realization
    return output
