"""Offline provenance/subdivision and full coupling checks; no Capytaine."""
from copy import deepcopy
from hashlib import sha256

import numpy as np
import pytest

from tests.experiments.free_surface.moving_canoe import make_domain, geometry
from tests.experiments.free_surface.capytaine_reference_stability import compare, FREQUENCIES
from tests.experiments.capytaine_mesh_cost import subdivide
from tests.experiments.free_surface.validate_radiation import MODES


@pytest.fixture(scope='module')
def inputs():
    domain = make_domain(12, False, 2)
    points, _, _ = geometry(domain, domain.rest.position, domain.rest.orientation,
                             np.zeros(len(domain.free_nodes)))
    triangles = points[domain.local_mapping][domain.contact.faces[domain.contact.body_faces]][:, ::-1]
    fine_triangles = subdivide(triangles)
    geometry_hash = lambda t: sha256(np.asarray(t, dtype='<f8').tobytes()).hexdigest()
    # Sentinel coefficients are built independently in physical units to test
    # parser/scaling/coupling behavior, not to claim a hydrodynamic solution.
    d = np.array([1., 1., 1., 3.5, 3.5, 3.5])
    units = 1000*3.5**3*np.outer(d, d)
    a = units*np.diag(np.arange(1., 7.))
    a[0, 4] = a[4, 0] = .3*units[0, 4]
    cases = [dict(omega=w, capytaine=dict(added_mass=a.tolist(),
             radiation_damping=(.2*w*a).tolist())) for w in FREQUENCIES]
    coarse = dict(mode_order=list(MODES), water_depth=4., cases=cases,
        rotation_center=domain.rest.position.tolist(), wetted_panels=704,
        wetted_triangles_sha256=geometry_hash(triangles), wet_subdivisions=2,
        capytaine_green_function_settings={'green_function': 'sentinel'})
    origin = deepcopy(coarse)
    origin.update(wetted_panels=2816, wetted_triangles_sha256=geometry_hash(fine_triangles), wet_subdivisions=3)
    fine = dict(frequencies=FREQUENCIES.copy(), water_depth=4., panels=2816,
        wetted_triangles_sha256=origin['wetted_triangles_sha256'],
        green_function_settings=coarse['capytaine_green_function_settings'],
        added_mass=[deepcopy(case['capytaine']['added_mass']) for case in cases],
        radiation_damping=[deepcopy(case['capytaine']['radiation_damping']) for case in cases])
    return coarse, fine, origin, triangles, fine_triangles


def test_actual_canonical_subdivision_and_complete_moment_scaling(inputs):
    result = compare(*inputs)
    assert result['qualified'] is False
    assert result['reference_stability_gate_passed']
    assert result['geometry']['ordered_subdivision_exact']
    for case in result['cases']:
        np.testing.assert_allclose(np.diag(case['reference_added_mass_nondimensional']), np.arange(1., 7.))
        assert case['reference_added_mass_nondimensional'][0][4] == pytest.approx(.3)


def test_one_off_diagonal_at_last_frequency_is_not_hidden_by_matching_diagonals(inputs):
    coarse, fine, origin, a, b = deepcopy(inputs)
    coarse['cases'][2]['capytaine']['added_mass'][0][4] *= -1
    result = compare(coarse, fine, origin, a, b)
    assert not result['reference_stability_gate_passed']
    gate = result['gates']['added_mass']
    assert (gate['worst_omega'], gate['worst_force_mode'], gate['worst_radiating_mode']) == (6., 'surge', 'pitch')
    assert gate['maximum_tolerance_ratio'] > 1


def test_hash_consistent_but_physically_changed_child_mesh_is_rejected(inputs):
    coarse, fine, origin, a, b = deepcopy(inputs)
    b[0, 0, 0] += 1e-4
    changed = sha256(np.asarray(b, dtype='<f8').tobytes()).hexdigest()
    fine['wetted_triangles_sha256'] = origin['wetted_triangles_sha256'] = changed
    with pytest.raises(ValueError, match='subdivision children within roundoff'):
        compare(coarse, fine, origin, a, b)


@pytest.mark.parametrize('key,value', [('rotation_center', [0., 0., .3]),
                                      ('water_depth', 5.), ('mode_order', list(reversed(MODES)))])
def test_inconsistent_origin_metadata_is_rejected(inputs, key, value):
    coarse, fine, origin, a, b = deepcopy(inputs)
    origin[key] = value
    with pytest.raises(ValueError, match='metadata differs'):
        compare(coarse, fine, origin, a, b)
