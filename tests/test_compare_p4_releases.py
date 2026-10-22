"""Offline endpoint conventions and genuine temporal-control ownership."""
from copy import deepcopy

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from tests.experiments.free_surface.compare_p4_releases import compare


def histories():
    report = dict(source_npz_sha256='same-initial-stage', reference_points_sha256='points',
        elements_sha256='elements', free_nodes_sha256='surface', eta_modes=2, duration=1.,
        derivative_step=1e-4, runner_sha256='runner', source_stage_metadata={'case': 'same'},
        rtol=1e-4, atol=1e-10, max_step=.01)
    state = np.zeros((3, 17))
    state[:, :3] = np.arange(3)[:, None]*np.array([.01, -.02, .03])
    state[:, 3:7] = Rotation.from_rotvec(np.arange(3)[:, None]*np.array([.02, .03, .04])).as_quat(scalar_first=True)
    state[:, 7:13] = np.arange(3)[:, None]*np.arange(1., 7.)/100
    data = dict(times=np.array([0., .5, 1.]), states=state, canonical_states=np.zeros((3, 16)),
        hamiltonian=np.array([.3, .301, .3005]), fluid_volume=np.array([10., 10.+1e-8, 10.+2e-8]),
        rim_clearance_lower_bound=np.array([.2, .19, .18]),
        weak_normal_relative_residual=np.full(3, 1e-12), fluid_work_identity_absolute_residual=np.full(3, 1e-14))
    tight_report = dict(report, rtol=1e-5, max_step=.005)
    return report, data, tight_report, deepcopy(data)


def test_quaternion_sign_is_not_a_physical_orientation_or_initial_state_error():
    loose_report, loose, tight_report, tight = histories()
    tight['states'][:, 3:7] *= -1
    output = compare(loose_report, loose, tight_report, tight)
    endpoint = output['exact_endpoint_comparison']
    assert endpoint['orientation_geodesic_difference_rad'] == pytest.approx(0., abs=1e-15)
    assert endpoint['initial_orientation_geodesic_difference_rad'] == 0
    assert endpoint['initial_nonquaternion_component_maximum_absolute_difference'] == 0
    np.testing.assert_allclose(output['loose']['maximum_world_relative_rotation_components_rad'], [.04, .06, .08])
    assert output['qualified'] is False


def test_raw_accepted_drift_and_exact_endpoint_differences_are_retained():
    loose_report, loose, tight_report, tight = histories()
    tight['states'][-1, 2] += 3e-5
    tight['hamiltonian'][1:] = [.3001, .30005]
    output = compare(loose_report, loose, tight_report, tight)
    assert output['loose']['maximum_accepted_hamiltonian_change_joules'] == pytest.approx(.001)
    assert output['tight']['maximum_accepted_hamiltonian_change_joules'] == pytest.approx(.0001)
    assert output['exact_endpoint_comparison']['position_difference_m'][2] == pytest.approx(-3e-5)
    np.testing.assert_allclose(output['loose']['maximum_world_velocity_change'], np.arange(1., 7.)/50)


@pytest.mark.parametrize('change', ['stage', 'shape-step', 'looser-control', 'no-tightening'])
def test_different_case_or_nonrefinement_is_rejected(change):
    loose_report, loose, tight_report, tight = histories()
    if change == 'stage':
        tight_report['source_npz_sha256'] = 'different'
    elif change == 'shape-step':
        tight_report['derivative_step'] *= .5
    elif change == 'looser-control':
        tight_report['atol'] *= 10
    else:
        tight_report = deepcopy(loose_report)
    with pytest.raises(ValueError):
        compare(loose_report, loose, tight_report, tight)
