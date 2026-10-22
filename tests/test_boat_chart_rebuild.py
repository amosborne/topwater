"""A real flat-water chine crossing, without a volume/fluid operator."""
import numpy as np
import pytest

from topwater import InitialState
from .experiments.free_surface.boat_chart_rebuild import diagnostics
from .experiments.free_surface.boat_contact_preflight import make_chart, pose
from .experiments.free_surface.canoe_contact import ContactTopologyError


def test_rebuild_handles_recorded_chine_transition_and_preserves_normal_work():
    body, rest, frozen = make_chart()
    position, orientation = pose(rest, np.zeros(3), np.array([0., .01, 0.]))
    with pytest.raises(ContactTopologyError):
        frozen.evaluate(position, orientation, np.zeros(len(frozen.free_nodes)))
    state = InitialState(position, orientation)
    _, _, rebuilt = make_chart(state=state)
    assert set(rebuilt.body_source_faces) != set(frozen.body_source_faces)
    result = diagnostics(body, state, rebuilt)
    assert abs(result['fluid_volume_error']) < 1e-10
    assert result['wet_parent_area_max_error'] < 1e-12
    assert abs(result['normal_work_error']) < 1e-12
    assert abs(result['constant_gauge_work_error']) < 1e-12
    assert result['rim_clearance_flat_water'] > 0
