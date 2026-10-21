"""Independent physical calm-water correlation; no wave/BEM accuracy claim."""

import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
import pytest

cpt = pytest.importorskip('capytaine')

from topwater import InitialState, SimulationInputs
from topwater.frames import rotation_matrix
from topwater.hydrostatics import equilibrate_heave
from .reference_hull import input_data


@pytest.mark.parametrize('refinement', [2, 3, 4])
@pytest.mark.parametrize('yaw', [0., .4])
def test_canoe_hydrostatic_correlation(refinement, yaw):
    assert cpt.__version__ == '3.0.0'
    data = SimulationInputs.from_json(json.dumps(input_data(refinement)))
    state = InitialState([0, 0, 0], [np.cos(yaw/2), 0, 0, np.sin(yaw/2)])
    rest, result = equilibrate_heave(data.body, state)
    vertices = ((data.body.mesh.vertices-data.body.cg)
                @ rotation_matrix(rest.orientation).T+rest.position)
    # Independent clipping of the original physical faces. No Topwater clipped
    # triangles are passed to the reference, and no artificial closure is wet.
    mesh = cpt.Mesh(vertices=vertices,
        faces=data.body.mesh.faces[data.body.mesh.physical_faces],
        auto_clean=False, auto_check=False, quadrature_method='Gauss-Legendre 2')
    body = cpt.FloatingBody(mesh=mesh, mass=data.body.mass,
        center_of_mass=rest.position,
        dofs=cpt.rigid_body_dofs(rotation_center=rest.position)).immersed_part()
    # Capytaine's default centroid rule is not exact for quadratic moments.
    # Its 2x2 Gauss rule integrates these planar-face polynomials exactly.
    assert body.mesh.quadrature_method == 'Gauss-Legendre 2'
    order = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']
    stiffness = body.compute_hydrostatic_stiffness(rho=1000., g=9.81).sel(
        influenced_dof=order, radiating_dof=order).values
    assert_allclose(result.displaced_volume, body.disp_volume, rtol=1e-10, atol=1e-12)
    assert_allclose(result.wetted_area, body.wet_surface_area, rtol=1e-10, atol=1e-12)
    assert_allclose(result.waterplane_area, body.waterplane_area, rtol=1e-10, atol=1e-12)
    assert_allclose(result.center_of_buoyancy, body.center_of_buoyancy, rtol=1e-10, atol=1e-12)
    assert_allclose(result.stiffness, stiffness, rtol=1e-10, atol=1e-8)
    report = dict(capytaine=cpt.__version__, refinement=refinement, yaw=yaw,
        displaced_volume=result.displaced_volume, waterplane_area=result.waterplane_area,
        center_of_buoyancy=result.center_of_buoyancy.tolist(),
        stiffness_max_absolute_error=float(np.max(np.abs(result.stiffness-stiffness))))
    destination = Path('artifacts')
    destination.mkdir(exist_ok=True)
    (destination/f'capytaine-hydrostatics-{refinement}-{yaw}.json').write_text(
        json.dumps(report, indent=2)+'\n', encoding='utf-8')
