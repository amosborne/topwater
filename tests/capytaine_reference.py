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


def radiation_columns(triangles, center, mass, angular_frequency, modes, *, density, gravity):
    """Deep-water radiation on EXACT supplied physical panels.

    ``center`` is the WORLD CG/rotation center, not the material-frame body CG.
    Triangles have BODY-outward winding. All six influenced force/moment rows
    are returned for only the requested radiating columns. No wave excitation,
    hydrostatic restoring or rigid inertia is included.

    Capytaine uses Re(amplitude*exp(-i*omega*t)) and unit displacement. The
    returned force_per_velocity instead corresponds to unit velocity amplitude:
    F=(-B+i*omega*A)U, hence mean body power is -U*.B.U/2.
    At angular_frequency=inf only A-infinity and zero damping are returned;
    finite harmonic force amplitudes are undefined and their keys are omitted.
    """
    import capytaine as cpt

    modes = tuple(modes)
    if (not modes or len(set(modes)) != len(modes) or any(mode not in MODES for mode in modes)):
        raise ValueError('distinct known rigid radiation modes required')
    density, gravity, angular_frequency = map(float, (density, gravity, angular_frequency))
    if not np.isfinite([density, gravity, mass]).all() or min(density, gravity, mass) <= 0:
        raise ValueError('density, gravity and mass must be positive finite scalars')
    if np.isnan(angular_frequency) or angular_frequency <= 0:
        raise ValueError('angular_frequency must be positive finite or infinity')
    triangles, center = np.asarray(triangles, dtype=float), np.asarray(center, dtype=float)
    if (triangles.ndim != 3 or triangles.shape[1:] != (3, 3) or not len(triangles)
            or center.shape != (3,) or not np.isfinite(triangles).all() or not np.isfinite(center).all()
            or np.any(triangles[:, :, 2] > 0) or np.any(triangles[:, :, 2].mean(axis=1) >= 0)):
        raise ValueError('finite wetted triangles below z=0 and a finite world CG required')
    body = body_from_triangles(triangles, center, mass)
    infinite = np.isinf(angular_frequency)
    green = (cpt.Delhommeau(tabulation_nr=0, tabulation_nz=0, tabulation_cache_dir=None)
             if infinite else cpt.Delhommeau(tabulation_cache_dir=None))
    added, damping = np.empty((6, len(modes))), np.empty((6, len(modes)))
    forces = None if infinite else np.empty((6, len(modes)), dtype=complex)
    with threadpool_limits(1):
        solver = cpt.BEMSolver(green_function=green)
        for column, mode in enumerate(modes):
            result = solver.solve(cpt.RadiationProblem(body=body, radiating_dof=mode,
                omega=angular_frequency, rho=density, g=gravity, water_depth=np.inf))
            added[:, column] = [result.added_mass[row] for row in MODES]
            damping[:, column] = [result.radiation_damping[row] for row in MODES]
            if not infinite:
                forces[:, column] = [result.forces[row] for row in MODES]
    if (not np.isfinite(added).all() or not np.isfinite(damping).all()
            or (forces is not None and not np.isfinite(forces).all())):
        raise FloatingPointError('nonfinite Capytaine radiation reference')
    reference = dict(added_mass=added, radiation_damping=damping,
        angular_frequency='infinite' if infinite else angular_frequency,
        influenced_modes=MODES, radiating_modes=modes,
        density=float(density), gravity=float(gravity), water_depth='infinite',
        green_function_settings=green.exportable_settings)
    if not infinite:
        reference.update(force_per_displacement=forces, force_per_velocity=1j*forces/angular_frequency)
    return reference


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
