"""Physical instantaneous added-mass correlation; finite-frequency work is future."""

import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
import pytest
from threadpoolctl import threadpool_limits

cpt = pytest.importorskip('capytaine')

from topwater import SimulationInputs
from topwater.hydrostatics import equilibrate_heave
from topwater.inertia import instantaneous_added_mass
from topwater.rankine import TrianglePanels
from .reference_hull import input_data


def compare_inertia(refinement):
    assert cpt.__version__ == '3.0.0'
    data = SimulationInputs.from_json(json.dumps(input_data(refinement)))
    rest, wet = equilibrate_heave(data.body, data.initial_state)
    panels = TrianglePanels(wet.wetted_triangles, -1)
    names = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']
    # Load the native kernel before discovering process BLAS/OpenMP controls.
    green = cpt.Delhommeau(tabulation_nr=0, tabulation_nz=0, tabulation_cache_dir=None)
    with threadpool_limits(1):
        actual = instantaneous_added_mass(panels, rest.position)
        # Clipping was independently correlated separately. Here both methods
        # use identical physical triangles, no artificial cap, CG and fluid.
        vertices, indices = np.unique(panels.triangles.reshape(-1, 3), axis=0, return_inverse=True)
        faces = indices.reshape(-1, 3)
        # Share exactly equal coordinates without moving/healing any point.
        # Disjoint vertex indices would make Capytaine's topology checks see
        # every triangle as a separate body component.
        assert_allclose(vertices[faces], panels.triangles, atol=0, rtol=0)
        mesh = cpt.Mesh(vertices=vertices, faces=faces, auto_clean=False, auto_check=False)
        assert_allclose(mesh.faces_normals, -panels.normals, atol=1e-13, rtol=1e-12)
        body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=rest.position),
                                mass=data.body.mass, center_of_mass=rest.position)
        solver = cpt.BEMSolver(green_function=green)
        reference = np.empty((6, 6))
        damping = np.empty((6, 6))
        for column, name in enumerate(names):
            result = solver.solve(cpt.RadiationProblem(body=body, radiating_dof=name,
                omega=np.inf, rho=1000., g=9.81, water_depth=np.inf))
            reference[:, column] = [result.added_mass[row] for row in names]
            damping[:, column] = [result.radiation_damping[row] for row in names]
    scale = np.array([1., 1., 1., 3.5, 3.5, 3.5])
    normalized = actual.added_mass/scale[:, None]/scale[None, :]/(1000*3.5**3)
    expected = reference/scale[:, None]/scale[None, :]/(1000*3.5**3)
    floor = max(1e-10, .005*np.max(np.abs(expected)))
    tolerance = .05*np.abs(expected)+floor
    assert np.all(np.abs(normalized-expected) <= tolerance)
    assert np.all(np.isfinite(actual.acceleration_potential))
    assert np.min(np.linalg.eigvalsh((normalized+normalized.T)/2)) > 0
    assert np.max(np.abs(normalized-normalized.T)) < .005*np.max(np.abs(normalized))
    assert_allclose(damping, 0, atol=0, rtol=0)
    report = dict(capytaine=cpt.__version__, comparison='instantaneous_added_mass',
        refinement=refinement, wetted_panel_count=len(panels), center_of_mass=rest.position.tolist(),
        added_mass=actual.added_mass.tolist(), capytaine_added_mass=reference.tolist(),
        normalized_max_absolute_error=float(np.max(np.abs(normalized-expected))),
        normalized_reciprocity_defect=float(np.max(np.abs(normalized-normalized.T))),
        normalized_energy_eigenvalues=np.linalg.eigvalsh((normalized+normalized.T)/2).tolist(),
        family_absolute_floor=floor, diagnostics=actual.diagnostics,
        scope='same-mesh infinite-frequency correlation; not geometry convergence or finite-frequency validation')
    destination = Path('artifacts')
    destination.mkdir(exist_ok=True)
    (destination/f'capytaine-inertia-{refinement}.json').write_text(
        json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return report


@pytest.mark.parametrize('refinement', [3, 4, 5])
def test_instantaneous_six_motion_inertia(refinement):
    compare_inertia(refinement)
