"""Different convex open-box boat inputs for generic solver smoke controls.

This deliberately simple physical input is not an alternative production hull
or a qualified seakeeping case. No solver imports this fixture.
"""
import numpy as np

from topwater import HullMesh, InitialState, RigidBody, SimulationInputs, sea_states
from topwater.constants import GRAVITY, WATER_DENSITY


def box_inputs():
    length, beam, depth, mass = 6.5, 2.4, 1.2, 2000.
    x, y = length/2, beam/2
    vertices = np.array([[-x,-y,0], [x,-y,0], [x,y,0], [-x,y,0],
                         [-x,-y,depth], [x,-y,depth], [x,y,depth], [-x,y,depth]])
    faces = np.array([[0,2,1], [0,3,2], [0,1,5], [0,5,4], [1,2,6], [1,6,5],
                      [2,3,7], [2,7,6], [3,0,4], [3,4,7], [4,5,6], [4,6,7]])
    mesh = HullMesh(vertices, faces, [4,5,6,7], [10,11])
    body = RigidBody(mesh, mass, [0.,0.,depth/2], np.diag([1200.,7300.,8000.]))
    # Independent rectangular displacement balance in mean water.
    draft = mass/(WATER_DENSITY*length*beam)
    state = InitialState(position=[.13,.21,depth/2-draft],
        linear_velocity=[.03,-.02,.01], angular_velocity=[.01,-.015,.02])
    waves = sea_states.realize(beam_m=beam, gravity=GRAVITY, seed=42)
    return SimulationInputs(body, waves, state)
