from types import SimpleNamespace

import capytaine as cpt
import numpy as np
import pytest

from topwater.solve import bem_checks


def test_preflight_checks_unique_complete_bank_once_without_solving(monkeypatch):
    events = []
    monkeypatch.setattr(bem_checks, 'LinearPotentialFlowProblem',
                        lambda **kwargs: SimpleNamespace(**kwargs))
    def observe(name):
        def check(problems):
            events.append((name, [problem.omega for problem in problems]))
        return check
    monkeypatch.setattr(bem_checks, '_check_wavelength_and_mesh_resolution', observe('mesh'))
    monkeypatch.setattr(bem_checks, '_check_wavelength_and_irregular_frequencies', observe('irregular'))
    monkeypatch.setattr(bem_checks, '_check_wavelength_and_water_depth',
                        lambda _: pytest.fail('infinite depth needs no finite-depth advice'))
    monkeypatch.setattr(cpt.BEMSolver, 'solve', lambda *a, **kw: pytest.fail('preflight ran BEM'))
    bank = bem_checks.check_frequency_bank(object(), [48., .1, .1, 12., np.inf],
                                           density=1025., gravity=9.81, depth=np.inf)
    expected = [.1, 12., 48., np.inf]
    assert events == [('mesh', expected), ('irregular', expected)]
    assert bank == frozenset(expected)


def test_prechecked_solver_skips_only_completed_checks_and_rejects_new_frequency(monkeypatch):
    calls = []
    def solve(self, problem, **kwargs):
        calls.append(kwargs)
        return 'result'
    monkeypatch.setattr(cpt.BEMSolver, 'solve', solve)
    solver = bem_checks.PrecheckedBEMSolver.__new__(bem_checks.PrecheckedBEMSolver)
    solver.checked_frequencies = frozenset([1., np.inf])
    for _ in range(6):
        assert solver.solve(SimpleNamespace(omega=1.)) == 'result'
    assert len(calls) == 6
    assert all(call['_check_wavelength'] is False for call in calls)
    assert all(call['keep_details'] is True for call in calls)
    with pytest.raises(ValueError, match='initial preflight bank'):
        solver.solve(SimpleNamespace(omega=2.))
    assert len(calls) == 6


@pytest.mark.parametrize('frequencies', [[], [0.], [-1.], [np.nan], [[1., 2.]]])
def test_invalid_preflight_bank_is_rejected(frequencies):
    with pytest.raises(ValueError, match='positive frequencies'):
        bem_checks.check_frequency_bank(object(), frequencies, density=1025., gravity=9.81, depth=np.inf)
