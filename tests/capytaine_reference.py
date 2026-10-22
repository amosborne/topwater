"""Validation-only Capytaine reference on exactly supplied wetted triangles."""
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
from threadpoolctl import threadpool_limits


MODES = ('Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw')


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


def instantaneous_inertia(triangles, center, mass, *, density=1000., gravity=9.81):
    """Infinite-depth A-infinity on supplied physical panels, in SI units."""
    import capytaine as cpt

    if not np.isfinite([density, gravity]).all() or min(density, gravity) <= 0:
        raise ValueError('density and gravity must be positive finite scalars')
    body = body_from_triangles(triangles, center, mass)
    # Discover native BLAS/OpenMP only after loading the Green-function library.
    green = cpt.Delhommeau(tabulation_nr=0, tabulation_nz=0, tabulation_cache_dir=None)
    with threadpool_limits(1):
        solver = cpt.BEMSolver(green_function=green)
        reference, damping = np.empty((6, 6)), np.empty((6, 6))
        for column, name in enumerate(MODES):
            result = solver.solve(cpt.RadiationProblem(body=body, radiating_dof=name,
                omega=np.inf, rho=density, g=gravity, water_depth=np.inf))
            reference[:, column] = [result.added_mass[row] for row in MODES]
            damping[:, column] = [result.radiation_damping[row] for row in MODES]
    return reference, damping


def frequency_response(triangles, center, mass, frequencies, headings, *, water_depth=np.inf):
    """Finite-frequency reference only; excitation includes incident pressure.

    Unit wave elevation at the world origin, Re(amplitude*exp(-i*omega*t)),
    specified water depth (infinite by default), zero forward speed. Forces/moments use world axes about CG.
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
    with threadpool_limits(1):
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
    return dict(added_mass=added_mass, radiation_damping=damping,
                green_function_settings=green.exportable_settings,
                diffraction=diffraction, froude_krylov=incident,
                excitation=diffraction+incident)
