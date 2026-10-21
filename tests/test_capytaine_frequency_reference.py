"""An independent pressure-limit check for finite-frequency reference plumbing."""

import json
import numpy as np
import pytest

pytest.importorskip('capytaine')

from topwater import InitialState, SimulationInputs
from topwater.hydrostatics import equilibrate_heave
from .capytaine_reference import frequency_response
from .reference_hull import input_data


def test_excitation_approaches_independent_quasistatic_water_level_load():
    data = SimulationInputs.from_json(json.dumps(input_data(2)))
    # Heel the same canoe to exercise moments and the rotation reference.
    state = InitialState(data.initial_state.position, [np.cos(.04), np.sin(.04), 0., 0.])
    rest, wet = equilibrate_heave(data.body, state)
    response = frequency_response(wet.wetted_triangles, rest.position, data.body.mass,
                                  [.04, .02], [0., np.pi/2])
    # As wavelength -> infinity, unit incident elevation shifts the water level
    # uniformly. Its total pressure load tends to the hydrostatic heave column.
    # Diffraction alone cannot pass this check; signs and CG moments matter.
    scale = np.array([1., 1., 1., 3.5, 3.5, 3.5])
    expected = wet.stiffness[:, 2]/scale
    errors = np.linalg.norm(response['excitation']/scale-expected, axis=-1)
    assert np.all(errors[1] < .3*errors[0])
    assert np.max(errors[1])/np.linalg.norm(expected) < 1e-4
