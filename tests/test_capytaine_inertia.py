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
from .capytaine_reference import instantaneous_inertia as reference_inertia


def compare_inertia(refinement):
    assert cpt.__version__ == '3.0.0'
    data = SimulationInputs.from_json(json.dumps(input_data(refinement)))
    rest, wet = equilibrate_heave(data.body, data.initial_state)
    panels = TrianglePanels(wet.wetted_triangles, -1)
    with threadpool_limits(1):
        actual = instantaneous_added_mass(panels, rest.position)
    reference, damping = reference_inertia(panels.triangles, rest.position, data.body.mass)
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
