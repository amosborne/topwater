"""Regression gates for the tolerance/physics/contention comparison mistake."""
from copy import deepcopy
import numpy as np
import pytest

from .experiments.free_surface.boat_benchmark_contract import performance_comparison, array_identity
from .experiments.free_surface.boat_integration_controls import BODY_ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE


def report():
    return dict(complete=True, achieved_duration_seconds=1., benchmark_seconds=40.,
        benchmark_contract=dict(physics='evolving sea and body', duration_seconds=1., state_count=43,
            input_sha256='owned-physical-input', integrator='RK23', max_step_seconds=.2, first_step_seconds=.2,
            error_norm='RK23 RMS over43states', state_layout=['body13','eta15','psi15'],
            solve_controls=dict(residual=1e-9),
            tolerances=dict(rtol=RELATIVE_TOLERANCE,
                absolute_by_component=np.r_[BODY_ABSOLUTE_TOLERANCE, np.full(15,.001), np.full(15,.01)].tolist()),
            spatial=dict(grid=129)),
        measurement_environment=dict(resource_status='uncontended', hardware=dict(cpu='targetCPU', gpu='targetGPU'),
            software=dict(python='same'), timing_scope='setup plus integration', startup_treatment='same',
            thread_pools=[dict(num_threads=1)], thread_environment=dict(OMP_NUM_THREADS='1'),
            output_policy='same', cache_policy='prewarmed current producers'))


@pytest.mark.parametrize('change', ['rtol', 'atol', 'physics', 'duration', 'state_count', 'spatial',
                                     'incomplete', 'short_run', 'contended', 'unknown_load', 'hardware', 'invalid_time'])
def test_changed_comparison_condition_never_produces_speedup(change):
    reference = report(); candidate = deepcopy(reference); candidate['benchmark_seconds'] = 1.
    contract = candidate['benchmark_contract']
    if change == 'rtol': contract['tolerances']['rtol'] = 1e-5
    elif change == 'atol': contract['tolerances']['absolute_by_component'][0] = 1e-7
    elif change == 'physics': contract['physics'] = 'fixed forced boundary and body'
    elif change == 'duration': contract['duration_seconds'] = .1; candidate['achieved_duration_seconds'] = .1
    elif change == 'state_count': contract['state_count'] = 13
    elif change == 'spatial': contract['spatial']['grid'] = 257
    elif change == 'incomplete': candidate['complete'] = False
    elif change == 'short_run': candidate['achieved_duration_seconds'] = .1
    elif change == 'contended': candidate['measurement_environment']['resource_status'] = 'contended'
    elif change == 'unknown_load': candidate['measurement_environment']['resource_status'] = 'uncontrolled'
    elif change == 'hardware': candidate['measurement_environment']['hardware']['gpu'] = 'differentGPU'
    elif change == 'invalid_time': candidate['benchmark_seconds'] = float('nan')
    result = performance_comparison(reference, candidate)
    assert not result['comparable'] and result['reasons'] and 'speedup' not in result


def test_only_matched_reports_allow_ratio_and_implementation_hashes_may_change():
    a = report(); b = deepcopy(a)
    a['source_sha256'] = 'reference'; b['source_sha256'] = 'optimized'
    b['benchmark_seconds'] = 20.
    assert performance_comparison(a, b)['speedup'] == 2.


def test_old_reports_are_not_reconstructed_into_eligible_benchmarks():
    a = report(); b = deepcopy(a); del b['benchmark_contract']
    assert not performance_comparison(a, b)['comparable']


def test_array_identity_owns_values_and_shape_without_archive_metadata():
    a = dict(nodes=np.arange(6.).reshape(2, 3), length=64.)
    assert array_identity(a) == array_identity(dict(reversed(list(a.items()))))
    assert array_identity(a) != array_identity(dict(nodes=np.arange(6.).reshape(3, 2), length=64.))
    assert array_identity(a) != array_identity(dict(nodes=np.arange(6.).reshape(2, 3), length=32.))
    assert array_identity(dict(x=np.array([0.]))) != array_identity(dict(x=np.empty((1,0))))
    assert array_identity(dict(x=64.)) != array_identity(dict(x=np.array([64.])))


@pytest.mark.parametrize('field,value', [('rtol', None), ('rtol', -1.), ('rtol', float('nan')),
                                      ('atol', 0.), ('atol', float('inf')), ('count', True), ('count', 0)])
def test_identically_invalid_contracts_do_not_allow_ratio(field, value):
    a = report()
    if field == 'count': a['benchmark_contract']['state_count'] = value
    elif field == 'rtol': a['benchmark_contract']['tolerances']['rtol'] = value
    else: a['benchmark_contract']['tolerances']['absolute_by_component'][0] = value
    assert not performance_comparison(a, deepcopy(a))['comparable']


def test_undeclared_cache_state_cannot_masquerade_as_controlled_startup():
    a = report(); a['measurement_environment']['cache_policy'] = 'uncontrolled persistent caches'
    assert not performance_comparison(a, deepcopy(a))['comparable']
