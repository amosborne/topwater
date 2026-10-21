"""Validation-only Capytaine reference on exactly supplied wetted triangles."""
import numpy as np
from numpy.testing import assert_allclose
from threadpoolctl import threadpool_limits


def instantaneous_inertia(triangles, center, mass):
    import capytaine as cpt

    assert cpt.__version__ == '3.0.0'
    names = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']
    # Discover native BLAS/OpenMP only after loading the Green-function library.
    green = cpt.Delhommeau(tabulation_nr=0, tabulation_nz=0, tabulation_cache_dir=None)
    with threadpool_limits(1):
        vertices, indices = np.unique(triangles.reshape(-1, 3), axis=0,
                                      return_inverse=True)
        faces = indices.reshape(-1, 3)
        # Exact shared coordinates, no geometry healing or disconnected faces.
        assert_allclose(vertices[faces], triangles, atol=0, rtol=0)
        mesh = cpt.Mesh(vertices=vertices, faces=faces, auto_clean=False, auto_check=False)
        normal = np.cross(triangles[:, 1]-triangles[:, 0],
                          triangles[:, 2]-triangles[:, 0])
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        assert_allclose(mesh.faces_normals, normal, atol=1e-13, rtol=1e-12)
        body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=center),
                                mass=mass, center_of_mass=center)
        solver = cpt.BEMSolver(green_function=green)
        reference, damping = np.empty((6, 6)), np.empty((6, 6))
        for column, name in enumerate(names):
            result = solver.solve(cpt.RadiationProblem(body=body, radiating_dof=name,
                omega=np.inf, rho=1000., g=9.81, water_depth=np.inf))
            reference[:, column] = [result.added_mass[row] for row in names]
            damping[:, column] = [result.radiation_damping[row] for row in names]
    return reference, damping
