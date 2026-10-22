"""Canonical-canoe pressure, acceleration reaction and world mechanics."""
import json

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.linalg import solve
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from topwater import InitialState, SimulationInputs
from topwater.frames import rotation_matrix
from topwater.hydrostatics import equilibrate_heave
from topwater.rankine import TrianglePanels
from .reference_hull import input_data
from .experiments.free_surface.body_pressure import evaluate_body_pressure
from .experiments.free_surface.canoe import body_load, fluid_domain
from .experiments.free_surface.operator import assemble_energy


@pytest.fixture(scope='module')
def canoe():
    data = SimulationInputs.from_json(json.dumps(input_data(2)))
    rotation = Rotation.from_rotvec([.025, -.017, .31])
    pose = InitialState([.073, -.109, 0.], np.roll(rotation.as_quat(), 1))
    rest, wet = equilibrate_heave(data.body, pose)
    mesh, body, free = fluid_domain(wet, (4, 2, 2))
    panels = TrianglePanels.from_mesh(mesh)
    load = body_load(mesh, panels, body, rest.position)
    fixed = np.unique(mesh.faces[free])
    unknown = np.setdiff1d(np.arange(len(mesh.vertices)), fixed)
    with threadpool_limits(1):
        energy, _ = assemble_energy(mesh.faces, panels)
        auxiliary = np.zeros_like(load)
        auxiliary[unknown] = solve(energy[np.ix_(unknown, unknown)],
                                   load[unknown], assume_a='pos')
    matrix = rotation_matrix(rest.orientation)
    return dict(data=data, rest=rest, wet=wet, mesh=mesh, body=body,
                auxiliary=auxiliary, load=load, rotation=matrix,
                inertia=matrix@data.body.inertia_cg@matrix.T)


def evaluate(canoe, *, velocity=None, potential=None, rate=None, moving=None,
             auxiliary=None, **changes):
    mesh = canoe['mesh']; count = len(mesh.vertices)
    parameters = dict(points=mesh.vertices, faces=mesh.faces, body_faces=canoe['body'],
        center=canoe['rest'].position, mass=canoe['data'].body.mass,
        inertia_world=canoe['inertia'],
        velocity_world=np.zeros(6) if velocity is None else velocity,
        potential=np.zeros(count) if potential is None else potential,
        potential_rate_ale=np.zeros(count) if rate is None else rate,
        mesh_velocity=np.zeros((count, 3)) if moving is None else moving,
        acceleration_potentials=canoe['auxiliary'] if auxiliary is None else auxiliary)
    parameters.update(changes)
    return evaluate_body_pressure(**parameters)


def test_canonical_release_has_one_hydrostatic_load_and_correct_acceleration_reaction(canoe):
    result = evaluate(canoe)
    body = canoe['data'].body
    hydrostatic = np.r_[canoe['wet'].force, canoe['wet'].moment]
    # Independent clipping/pressure moments supply the static pressure term.
    assert_allclose(result.known_pressure_load, hydrostatic, rtol=1e-12, atol=2e-11)
    # The existing exact P1 weak rigid-motion load supplies the reaction map.
    expected_added = 1000*canoe['load'].T@canoe['auxiliary']
    assert_allclose(result.added_mass, expected_added, rtol=1e-11, atol=2e-12)
    rigid = np.zeros((6, 6)); rigid[:3, :3] = body.mass*np.eye(3)
    rigid[3:, 3:] = canoe['inertia']
    expected_acceleration = np.linalg.solve(rigid+expected_added,
        hydrostatic+np.array([0., 0., -body.mass*9.81, 0., 0., 0.]))
    assert_allclose(result.acceleration, expected_acceleration, rtol=1e-11, atol=2e-13)
    assert_allclose(result.hydrodynamic_load,
                    hydrostatic-expected_added@expected_acceleration, atol=2e-11)


@pytest.mark.parametrize('at_quadrature', [False, True])
def test_uniform_translation_ale_transport_preserves_hydrostatic_pressure(canoe, at_quadrature):
    points = canoe['mesh'].vertices
    translation = np.array([.37, -.23, 0.])
    # Exact uniform-flow potential: phi=U.x-|U|^2 t/2. Following moving
    # nodes gives phi_dot_ALE=|U|^2/2, whose transport restores Bernoulli.
    result = evaluate(canoe, velocity=np.r_[translation, np.zeros(3)],
        potential=points@translation+2.7,
        rate=np.full((len(canoe['body']), 4) if at_quadrature else len(points),
                     .5*translation@translation),
        moving=np.broadcast_to(translation, points.shape))
    static = evaluate(canoe)
    assert_allclose(result.known_pressure_load, static.known_pressure_load,
                    rtol=1e-12, atol=2e-11)
    assert_allclose(result.acceleration, static.acceleration, atol=2e-13)


def test_rotating_canoe_obeys_body_euler_equation_and_world_frame_covariance(canoe):
    points = canoe['mesh'].vertices
    center = canoe['rest'].position
    velocity = np.array([.12, -.07, .03, .13, -.09, .17])
    potential = canoe['auxiliary']@velocity
    moving = velocity[:3]+np.cross(velocity[3:], points-center)
    result = evaluate(canoe, velocity=velocity, potential=potential, moving=moving)
    rotation = canoe['rotation']; inertia = canoe['data'].body.inertia_cg
    omega_body = rotation.T@velocity[3:]
    expected_alpha_body = np.linalg.solve(inertia,
        rotation.T@result.hydrodynamic_load[3:]-np.cross(omega_body, inertia@omega_body))
    assert_allclose(rotation.T@result.acceleration[3:], expected_alpha_body, atol=2e-13)
    body = canoe['data'].body
    assert_allclose(body.mass*result.acceleration[:3],
        result.hydrodynamic_load[:3]+[0., 0., -body.mass*9.81], atol=2e-11)

    yaw = Rotation.from_rotvec([0., 0., .43]).as_matrix()
    transform = np.zeros((6, 6)); transform[:3, :3] = yaw; transform[3:, 3:] = yaw
    offset = np.array([-.219, .173, 0.])
    transformed = evaluate(canoe, points=points@yaw.T+offset,
        center=center@yaw.T+offset, inertia_world=yaw@canoe['inertia']@yaw.T,
        velocity=transform@velocity, potential=potential, moving=moving@yaw.T,
        auxiliary=canoe['auxiliary']@transform.T)
    assert_allclose(transformed.added_mass, transform@result.added_mass@transform.T,
                    rtol=1e-11, atol=2e-12)
    assert_allclose(transformed.hydrodynamic_load, transform@result.hydrodynamic_load,
                    rtol=1e-11, atol=2e-11)
    assert_allclose(transformed.acceleration, transform@result.acceleration, atol=2e-13)
