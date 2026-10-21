"""Physical geometry/weak-load checks for the contained coupled-domain candidate."""
import json

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from topwater import InitialState, SimulationInputs
from topwater.hydrostatics import equilibrate_heave
from topwater.rankine import TrianglePanels
from .reference_hull import input_data, surface_integrals
from .experiments.free_surface.canoe import body_load, fluid_domain, solve_inertia


@pytest.mark.parametrize('refinement', [3, 4])
@pytest.mark.parametrize('rotvec,offset', [([0, 0, 0], [0, 0, 0]),
    ([0, 0, .4], [.173, -.219, 0]), ([.08, .025, .4], [.173, -.219, 0])])
def test_canoe_hole_closes_and_weak_load_recovers_pressure(refinement, rotvec, offset):
    data = SimulationInputs.from_json(json.dumps(input_data(refinement)))
    pose = InitialState(offset, np.roll(Rotation.from_rotvec(rotvec).as_quat(), 1))
    rest, wet = equilibrate_heave(data.body, pose)
    mesh, body, free = fluid_domain(wet)
    panels = TrianglePanels.from_mesh(mesh)
    # HullMesh already verifies exact shared-edge incidence and orientation.
    # These independent integrals reject overlapping/missing surface regions.
    volume, _ = surface_integrals(mesh.vertices, mesh.faces)
    assert_allclose(volume, 32*np.pi**2-wet.displaced_volume, rtol=1e-13)
    assert_allclose(panels.areas[free].sum(), 8*np.pi**2-wet.waterplane_area, rtol=1e-13)
    assert_allclose(panels.triangles[body, ::-1], wet.wetted_triangles, rtol=0, atol=0)
    load = body_load(mesh, panels, body, rest.position)
    # A linear hydrostatic pressure, integrated through the same nodal load
    # map as motion forcing, must recover independent physical force/moment.
    assert_allclose(load.T@mesh.vertices[:, 2],
                    -np.r_[wet.force, wet.moment]/9810, atol=1e-13, rtol=1e-12)
    # Uniform potential tests the complete load map, including waterline rows.
    expected = np.r_[0., 0., wet.waterplane_area, wet.stiffness[2, 3:5]/9810, 0.]
    assert_allclose(load.sum(axis=0), expected, atol=1e-13, rtol=1e-12)


def test_same_domain_affine_flow_and_loads_refine_with_canoe_present():
    data = SimulationInputs.from_json(json.dumps(input_data(4)))
    rest, wet = equilibrate_heave(data.body, data.initial_state)
    rows = []
    coupling = []
    with threadpool_limits(1):
        for resolution in ((8, 4, 4), (16, 8, 4)):
            mesh, body, free = fluid_domain(wet, resolution)
            _, affine, diagnostics = solve_inertia(mesh, body, free, rest.position)
            rows.append(affine)
            coupling.append(diagnostics['coupling'])
    # phi=z satisfies zero FS potential and q=n_z on every other face,
    # including the bottom. Exact Dirichlet energy (twice unit-density kinetic
    # energy) is the actual fluid volume.
    for metric in ('potential_max_depth_normalized', 'energy_relative_error',
                   'exact_trace_energy_relative_error'):
        assert abs(rows[1][metric]) < abs(rows[0][metric])/2
    assert rows[1]['potential_max_depth_normalized'] < .001
    assert abs(rows[1]['energy_relative_error']) < .0001
    expected = np.r_[wet.force, wet.moment]/9810
    scale = np.r_[np.ones(3), np.full(3, 3.5)]*wet.displaced_volume
    errors = [np.linalg.norm((np.array(row['pressure_load'])-expected)/scale) for row in rows]
    assert errors[1] < errors[0]/2
    assert errors[1] < .01
    # Nonzero FS potential and body/wall forcing share one solve. Independently
    # known fields x,y,z,1 have elevation rates 0,0,1,0 over the WHOLE surface.
    for row in coupling:
        assert_allclose(row['integrated_flux'], row['expected_integrated_flux'], atol=1e-11, rtol=0)
        assert row['whole_surface_rate_rms'][3] < 1e-12
        assert row['potential_max_length_normalized'][3] < 1e-12
        assert row['exact_trace_load_max_length_normalized'] < 1e-13
    coarse = np.array(coupling[0]['whole_surface_rate_rms'][:3])
    fine = np.array(coupling[1]['whole_surface_rate_rms'][:3])
    assert np.all(fine < coarse)
    assert np.max(fine) < .002
    assert coupling[1]['solved_load_max_length_normalized'] < coupling[0]['solved_load_max_length_normalized']/2
