"""Actual Capytaine approximation controls; no panels, BEM or fluid solve."""
from copy import deepcopy
import json

import numpy as np
import pytest

from .capytaine_reference import _prony_realization

cpt = pytest.importorskip('capytaine')


def green_without_tabulation():
    # Exercise the installed real cached decomposition, without constructing its
    # large Green table or evaluating any Green function/PDE.
    green = cpt.Delhommeau.__new__(cpt.Delhommeau)
    green._hash = object.__hash__(green)
    green.finite_depth_prony_decomposition_method = 'python'
    return green


def test_actual_seeded_prony_repeats_coefficients_and_preserves_caller_rng():
    from capytaine.tools import prony_decomposition as module
    original = module.RNG
    state = deepcopy(original.bit_generator.state)
    results = []
    for _ in range(2):
        green = green_without_tabulation()
        with _prony_realization(green, 0) as metadata:
            finite = green.find_best_exponential_decomposition(1.2)
            infinite = green.find_best_exponential_decomposition(np.inf)
            # A second request is genuinely served by Capytaine's existing cache.
            np.testing.assert_array_equal(green.find_best_exponential_decomposition(1.2), finite)
        assert module.RNG is original and original.bit_generator.state == state
        assert 'find_best_exponential_decomposition' not in green.__dict__
        assert len(metadata['fitted_decompositions']) == 2
        assert metadata['initial_rng_state'] != metadata['final_rng_state']
        assert metadata['fitted_decompositions'][1]['dimensionless_wavenumber'] == 'infinity'
        np.testing.assert_array_equal(metadata['fitted_decompositions'][0]['lambda_and_amplitude'], finite)
        np.testing.assert_array_equal(metadata['fitted_decompositions'][1]['lambda_and_amplitude'], infinite)
        json.dumps(metadata, allow_nan=False)
        results.append((finite, infinite, metadata))
    np.testing.assert_array_equal(results[0][0], results[1][0])
    np.testing.assert_array_equal(results[0][1], results[1][1])
    assert results[0][2] == results[1][2]


def test_seeded_prony_restores_rng_and_method_when_reference_fails():
    from capytaine.tools import prony_decomposition as module
    green = green_without_tabulation()
    original, state = module.RNG, deepcopy(module.RNG.bit_generator.state)
    failure = RuntimeError('later reference calculation failed')
    with pytest.raises(RuntimeError) as caught:
        with _prony_realization(green, 17) as metadata:
            green.find_best_exponential_decomposition(1.2)
            raise failure
    assert caught.value is failure
    assert module.RNG is original and module.RNG.bit_generator.state == state
    assert 'find_best_exponential_decomposition' not in green.__dict__
    assert 'final_rng_state' in metadata and len(metadata['fitted_decompositions']) == 1


def test_absent_seed_keeps_legacy_rng_and_method_untouched():
    from capytaine.tools import prony_decomposition as module
    green = green_without_tabulation()
    original, method = module.RNG, green.find_best_exponential_decomposition
    with _prony_realization(green, None) as metadata:
        assert module.RNG is original
        assert green.find_best_exponential_decomposition == method
        assert metadata is None
    assert 'find_best_exponential_decomposition' not in green.__dict__


@pytest.mark.parametrize('seed', [-1, True, 1.5])
def test_invalid_seed_rejected_without_rng_change(seed):
    from capytaine.tools import prony_decomposition as module
    green = green_without_tabulation()
    original = module.RNG
    with pytest.raises(ValueError):
        with _prony_realization(green, seed):
            pytest.fail('invalid seed entered reference calculation')
    assert module.RNG is original
