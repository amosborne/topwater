"""Physical checks on the shared canoe, independently of the BEM operators."""

from dataclasses import replace
import json

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.spatial.transform import Rotation

from topwater import InitialState, SimulationInputs
from topwater.hydrostatics import (
    FloodingError, LossOfContactError, NoEquilibriumError,
    equilibrate_heave, hydrostatics,
)
from .reference_hull import SEMI_AXES, WATERLINE_HEIGHT, input_data


def reference(refinement=3):
    return SimulationInputs.from_json(json.dumps(input_data(refinement)))


def orientation(rotvec):
    return np.roll(Rotation.from_rotvec(rotvec).as_quat(), 1)


def test_independent_ellipsoid_geometry_converges():
    a, b, c = SEMI_AXES
    z = WATERLINE_HEIGHT
    volume = np.pi*a*b*(z-z**3/(3*c*c)+2*c/3)
    first_z = np.pi*a*b*(z*z/2-z**4/(4*c*c)-c*c/4)
    area = np.pi*a*b*(1-z*z/(c*c))
    waterplane_second = np.pi*a*b*(1-z*z/(c*c))**2*np.array([b*b, a*a])/4
    errors = []
    for refinement in (3, 4, 5):
        data = reference(refinement)
        result = hydrostatics(data.body, data.initial_state)
        expected_cb_z = first_z/volume + data.initial_state.position[2]-data.body.cg[2]
        expected_restoring = 1000*9.81*(waterplane_second+first_z-volume*data.body.cg[2])
        errors.append(np.abs(np.array([
            result.displaced_volume/volume-1,
            result.waterplane_area/area-1,
            result.center_of_buoyancy[2]/expected_cb_z-1,
            result.stiffness[3, 3]/expected_restoring[0]-1,
            result.stiffness[4, 4]/expected_restoring[1]-1,
        ])))
    assert np.all(np.diff(errors, axis=0) < 0)
    assert np.max(errors[-1]) < .01


def test_heave_equilibrium_balances_loaded_canoe():
    data = reference(4)
    rest, result = equilibrate_heave(data.body, data.initial_state)
    assert_allclose(result.force, [0, 0, data.body.mass*9.81], rtol=1e-11, atol=1e-10)
    assert_allclose(result.moment, 0, atol=1e-10)
    assert_allclose(rest.position[2], data.initial_state.position[2], atol=.001)
    assert_allclose(result.stiffness, result.stiffness.T, atol=1e-10)
    assert not np.isin(result.source_faces, data.body.mesh.artificial_faces).any()


def test_pressure_and_world_load_derivative_at_heeled_pose():
    data = reference()
    state = InitialState(data.initial_state.position+[1.2, -.7, .01],
                         orientation([.08, .025, .4]))
    result = hydrostatics(data.body, state)
    expected_force = np.array([0., 0., 1000*9.81*result.displaced_volume])
    assert_allclose(result.force, expected_force, atol=1e-10)
    assert_allclose(result.moment, np.cross(result.center_of_buoyancy-state.position,
                                          expected_force), atol=1e-10)
    step = 1e-6
    derivative = np.empty((6, 6))
    rotation = Rotation.from_quat(np.roll(state.orientation, -1))
    for axis in range(6):
        loads = []
        for sign in (-1, 1):
            offset = np.eye(3)[axis % 3]*step*sign
            moved = (replace(state, position=state.position+offset) if axis < 3 else
                     replace(state, orientation=np.roll(
                         (Rotation.from_rotvec(offset)*rotation).as_quat(), 1)))
            value = hydrostatics(data.body, moved)
            loads.append(np.r_[value.force, value.moment])
        derivative[:, axis] = -(loads[1]-loads[0])/(2*step)
    assert_allclose(result.stiffness, derivative, atol=2e-5, rtol=1e-7)


def test_failure_boundaries_and_overloading():
    data = reference()
    with pytest.raises(FloodingError):
        hydrostatics(data.body, InitialState([0, 0, data.body.cg[2]]))
    with pytest.raises(LossOfContactError):
        hydrostatics(data.body, InitialState([0, 0, 2]))
    with pytest.raises(NoEquilibriumError):
        equilibrate_heave(replace(data.body, mass=1000))


def test_genuinely_closed_submerged_hull_has_no_heave_restoring():
    from .reference_hull import surface_integrals
    data = reference()
    mesh = replace(data.body.mesh, opening_rim=(), artificial_faces=())
    body = replace(data.body, mesh=mesh)
    state = InitialState([0, 0, -2])
    result = hydrostatics(body, state)
    volume, center = surface_integrals(mesh.vertices, mesh.faces)
    assert_allclose(result.displaced_volume, volume, rtol=1e-12)
    assert_allclose(result.center_of_buoyancy, center-body.cg+state.position, atol=1e-12)
    assert result.waterplane_area == 0
    assert result.stiffness[2, 2] == 0
    assert result.rim_clearance is None


def test_equilibrium_is_independent_of_world_water_altitude():
    data = reference()
    rest, result = equilibrate_heave(data.body)
    shifted, value = equilibrate_heave(data.body, water_level=1e9)
    # The remaining error is representation of the final world coordinate,
    # rather than an absolute-altitude-dependent root stopping criterion.
    assert_allclose(shifted.position[2]-1e9, rest.position[2], atol=1e-7, rtol=0)
    assert_allclose(value.force, result.force, atol=.001, rtol=0)


def test_waterline_crossing_mesh_vertices_preserves_volume_derivative():
    data = reference()
    # Cross an actual horizontal row of hull vertices, not a generic smooth cut.
    levels = np.unique(data.body.mesh.vertices[:, 2])
    level = levels[len(levels)//2]
    state = InitialState([0, 0, data.body.cg[2]-level])
    center = hydrostatics(data.body, state)
    step = 1e-7
    below = hydrostatics(data.body, replace(state, position=state.position-[0, 0, step]))
    above = hydrostatics(data.body, replace(state, position=state.position+[0, 0, step]))
    assert_allclose((below.displaced_volume-above.displaced_volume)/(2*step),
                    center.waterplane_area, rtol=1e-6)
