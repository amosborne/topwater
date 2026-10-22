"""Offline radiation family checks include couplings and dimensional moments."""
from copy import deepcopy
import numpy as np
import pytest
from tests.experiments.free_surface.validate_radiation import aggregate, MODES, SETTINGS


def reports(omega=4.):
    # Build physical matrices from known nondimensional entries independently.
    d = np.diag([1., 1., 1., 3.5, 3.5, 3.5])
    a = 1000*3.5**3*(d@np.diag(np.arange(1., 7.))@d)
    a[0, 4] = a[4, 0] = 1000*3.5**4*.3
    b = omega*a*.2
    result = []
    for column, mode in enumerate(MODES):
        row = {key: 1 for key in SETTINGS}
        row.update(mode=mode, mode_order=list(MODES), omega=omega,
            amplitude_unit='m' if column < 3 else 'rad', periods=6., ramp_periods=2., fit_periods=3,
            n=24, material_subdivisions=2, water_depth=4.,
            column_harmonic_fits=[dict(added_mass=a[i, column], radiation_damping=b[i, column]) for i in range(6)],
            capytaine=dict(added_mass=a.tolist(), radiation_damping=b.tolist()))
        result.append(row)
    return result


def test_complete_matrices_scale_moments_and_use_global_frequency_family_floor():
    inputs = reports()+reports(2.)
    # The larger second frequency controls the floor at BOTH frequencies.
    for report in inputs[6:]:
        for key in ('added_mass', 'radiation_damping'):
            report['capytaine'][key] = (10*np.array(report['capytaine'][key])).tolist()
        for row in report['column_harmonic_fits']:
            row['added_mass'] *= 10
            row['radiation_damping'] *= 10
    output = aggregate(inputs)
    assert output['coefficient_gates_passed']
    assert output['qualified'] is False
    first = output['cases'][1]
    np.testing.assert_allclose(np.diag(first['added_mass_nondimensional']), np.arange(1., 7.))
    assert first['added_mass_nondimensional'][0][4] == pytest.approx(.3)
    assert output['gates']['added_mass']['family_floor'] == pytest.approx(.3)
    assert output['unresolved']


def test_gross_off_diagonal_failure_is_not_hidden_by_matching_diagonals():
    inputs = reports()
    # Surge force from pitch, with all six diagonals left exactly correct.
    inputs[4]['column_harmonic_fits'][0]['added_mass'] *= -1
    result = aggregate(inputs)
    assert not result['coefficient_gates_passed']
    gate = result['gates']['added_mass']
    assert (gate['worst_force_mode'], gate['worst_radiating_mode']) == ('surge', 'pitch')
    assert gate['maximum_tolerance_ratio'] > 1


@pytest.mark.parametrize('failure', ['missing', 'duplicate', 'settings', 'reference', 'fit'])
def test_incomplete_or_inconsistent_runs_are_rejected(failure):
    inputs = reports()
    if failure == 'missing':
        inputs.pop()
    elif failure == 'duplicate':
        inputs.append(deepcopy(inputs[0]))
    elif failure == 'settings':
        inputs[3]['sponge_strength'] = 9.
    elif failure == 'reference':
        inputs[2]['capytaine']['added_mass'][1][0] = 1.
    else:
        inputs[5]['column_harmonic_fits'] = None
    with pytest.raises(ValueError):
        aggregate(inputs)


def test_reference_only_smoke_and_entirely_zero_family_floor():
    inputs = reports()
    reference = deepcopy(inputs[0]); reference['periods'] = .1
    for report in inputs:
        report.pop('capytaine')
        for row in report['column_harmonic_fits']:
            row['radiation_damping'] = 0.
    reference['capytaine']['radiation_damping'] = np.zeros((6, 6)).tolist()
    result = aggregate(inputs, references=[reference])
    assert result['gates']['damping']['passed']
    assert result['gates']['damping']['family_floor'] == 1e-10
    inputs[0]['column_harmonic_fits'][4]['radiation_damping'] = 1.
    assert not aggregate(inputs, references=[reference])['gates']['damping']['passed']
