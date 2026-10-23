"""Current canonical chart derivatives, checked against owned geometry replay."""
import json

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.spatial.transform import Rotation

from topwater import InitialState
from tests.reference_boat import equilibrium
from .experiments.free_surface.boat_chart_tangent import chart_tangent
from .experiments.free_surface.boat_geometry_template import InvalidTemplate, evaluate_template
from .reference_boat import boundary_snapshot


@pytest.fixture(scope='module')
def generic_chart(tmp_path_factory):
    base = equilibrium()[0]
    pose = InitialState(np.array([.123, .087, base.position[2]]), base.orientation)
    _, reference, _ = boundary_snapshot(tmp_path_factory.mktemp('chart-tangent')/'generic.json',
        pose=pose, export_template=True)
    return reference, pose


def _mapping(current):
    result = np.zeros((len(current['free_nodes']), len(current['parent_nodes'])))
    for slot in range(3):
        ids = current['free_parent_indices'][:, slot]
        selected = ids >= 0
        result[np.flatnonzero(selected), ids[selected]] += current['free_parent_weights'][selected, slot]
    return result


def _moved(reference, pose, direction, step):
    orientation = (Rotation.from_rotvec(step*direction[3:6])
        *Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
    return evaluate_template(reference,
        pose=InitialState(pose.position+step*direction[:3], orientation),
        parent_eta=reference['parent_nodes'][:, 2]+step*direction[6:],
        parent_psi=reference['parent_psi'])[0]


def test_all_21_columns_and_mixed_direction_match_owned_h_h2(generic_chart):
    reference, pose = generic_chart
    current, tangent = chart_tangent(reference, pose=pose,
        parent_eta=reference['parent_nodes'][:, 2], parent_psi=reference['parent_psi'])
    assert_allclose(current['points'], reference['points'], atol=1e-12, rtol=0)
    assert tangent['point_jacobian'].shape == (*current['points'].shape, 21)
    assert tangent['parent_map_jacobian'].shape == (len(current['free_nodes']), 15, 21)
    rng = np.random.default_rng(71)
    directions = [*np.eye(21), rng.normal(size=21)/np.sqrt(21)]
    errors = []
    for direction in directions:
        for halving in range(12):
            h = 1e-5/2**halving
            try:
                moved = [_moved(reference, pose, direction, offset)
                         for offset in (h, -h, h/2, -h/2)]
            except InvalidTemplate:
                continue
            break
        else:
            pytest.fail('no valid owned finite-difference stencil')
        for key in ('faces', 'node_chart_key', 'free_parent_indices'):
            assert all(np.array_equal(part[key], current[key]) for part in moved)
        expected_points = tangent['point_jacobian']@direction
        expected_map = tangent['parent_map_jacobian']@direction
        for level in range(2):
            plus, minus = moved[2*level:2*level+2]
            step = h/2**level
            point_fd = (plus['points']-minus['points'])/(2*step)
            map_fd = (_mapping(plus)-_mapping(minus))/(2*step)
            assert_allclose(point_fd, expected_points, atol=2e-7, rtol=2e-6)
            assert_allclose(map_fd, expected_map, atol=2e-7, rtol=2e-6)
        errors.append(float(np.max(abs(point_fd-expected_points))))
    assert_allclose(tangent['parent_map_jacobian'].sum(axis=1), 0., atol=1e-12)
    print('all-column/mixed maximum fine-stencil coordinate error:', max(errors))


def test_rejects_default_symmetry_derivative_separation(tmp_path):
    _, reference, _ = boundary_snapshot(tmp_path/'default.json', export_template=True)
    pose = InitialState(reference['pose_position'], reference['pose_orientation'])
    # The center itself replays; its separating all-direction tangent must reject.
    evaluate_template(reference, pose=pose, parent_eta=reference['parent_nodes'][:, 2],
        parent_psi=reference['parent_psi'])
    with pytest.raises(InvalidTemplate, match='derivatives separated'):
        chart_tangent(reference, pose=pose, parent_eta=reference['parent_nodes'][:, 2],
            parent_psi=reference['parent_psi'])


def test_rejects_relevant_incidence_derivative_violation(generic_chart, monkeypatch):
    from .experiments.free_surface import boat_chart_tangent as module
    reference, pose = generic_chart
    signs, relevant = reference['template_reference_signs'], reference['template_relevant_planes']
    triples = reference['template_owner_triples']
    candidates = [(row, int(plane)) for row in range(len(triples))
                  for plane in np.flatnonzero((signs[row] == 0)&relevant[row])
                  if plane not in triples[row]]
    assert candidates, 'actual chart must carry overdetermined relevant incidence'
    plane = candidates[0][1]
    aliases = json.loads(str(reference['template_plane_aliases']))
    selected = {json.dumps(label): orientation for label, orientation in aliases[plane]}
    original = module._plane_tangent
    def changed(label, *args):
        result = original(label, *args)
        key = json.dumps(label)
        if key in selected:
            result[3, 6] += selected[key]
        return result
    monkeypatch.setattr(module, '_plane_tangent', changed)
    with pytest.raises(InvalidTemplate, match='relevant incident-plane derivatives separated'):
        chart_tangent(reference, pose=pose, parent_eta=reference['parent_nodes'][:, 2],
            parent_psi=reference['parent_psi'])


def test_retains_base_chart_ownership_guard(generic_chart):
    reference, pose = generic_chart
    altered = dict(reference, faces=reference['faces'][:, ::-1])
    with pytest.raises(InvalidTemplate, match='reindexed'):
        chart_tangent(altered, pose=pose, parent_eta=reference['parent_nodes'][:, 2],
            parent_psi=reference['parent_psi'])
