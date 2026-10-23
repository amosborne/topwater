"""Canonical-boat ALE pressure under a nonzero exact co-moving flow."""
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from topwater import InitialState
from topwater.frames import rotation_matrix
from .reference_boat import GRAVITY, WATER_DENSITY, equilibrium, reference_body
from .experiments.free_surface.boat_geometry_template import evaluate_template
from .experiments.free_surface.boat_world_boundary import construct
from .experiments.free_surface.body_pressure import evaluate_body_pressure


def test_uniform_horizontal_coflow_has_only_hydrostatic_pressure(tmp_path):
    body = reference_body()
    pose, hydrostatic = equilibrium()
    # Avoid the symmetric hull centerline coinciding with parent grid edges;
    # the manufactured field is translation invariant.
    pose = InitialState(pose.position+np.array([.123, .087, 0.]), pose.orientation)
    eta = np.zeros(15)
    _, arrays, _ = construct(tmp_path/'coflow.json', pose=pose, parent_eta=eta,
        parent_psi=np.zeros(15), export_template=True)
    velocity = np.array([.3, -.2, 0.])
    h = 1e-5
    moved = [evaluate_template(arrays, pose=InitialState(pose.position+dt*velocity, pose.orientation),
        parent_eta=eta, parent_psi=np.zeros(15))[0] for dt in (h, -h)]
    for current in moved:
        assert_array_equal(current['faces'], arrays['faces'])
    mesh_velocity = (moved[0]['points']-moved[1]['points'])/(2*h)
    rotation = rotation_matrix(pose.orientation)
    # phi=a.x-|a|^2*t/2 has exact uniform velocity a. Its Eulerian
    # time derivative cancels kinetic pressure for a co-moving hull.
    # Rates follow actual sliding chart nodes, not rigid material labels.
    potential = arrays['points']@velocity
    rate = mesh_velocity@velocity-velocity@velocity/2
    result = evaluate_body_pressure(arrays['points'], arrays['faces'], arrays['body_faces'],
        pose.position, body.mass, rotation@body.inertia_cg@rotation.T, np.r_[velocity, np.zeros(3)],
        potential, rate, mesh_velocity, np.zeros((len(potential), 6)), density=WATER_DENSITY, gravity=GRAVITY)
    # Independent hydrostatic clipping/integration, not another call to the
    # same pressure routine or an assertion of rigid-body equation balance.
    assert_allclose(result.known_pressure_load, np.r_[hydrostatic.force, hydrostatic.moment], atol=2e-8, rtol=1e-12)
    assert_allclose(result.acceleration, 0., atol=1e-10)
