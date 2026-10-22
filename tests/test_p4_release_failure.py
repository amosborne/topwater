"""Accepted-history preservation with cheap analytic ODEs, no fluid solve."""
import json

import numpy as np
import pytest
from scipy.integrate import RK23

from tests.experiments.free_surface.p4_canoe_release import _advance, _check_failure_output


def test_later_failure_preserves_only_accepted_endpoints_and_reraises_original(tmp_path):
    error = RuntimeError('synthetic later physical failure')

    def rate(time, state):
        if time > .25:
            raise error
        return np.ones_like(state)

    solver = RK23(rate, 0., [0.], 1., first_step=.1, max_step=.1)
    times, states = [], []
    path = tmp_path/'release-failed.json'

    def record():
        times.append(solver.t); states.append(solver.y.copy())

    def save_failure(caught):
        np.savez(path.with_suffix('.npz'), times=times, states=states)
        path.write_text(json.dumps(dict(incomplete=True, last_accepted_time=times[-1],
                                       error_type=type(caught).__name__, error_message=str(caught))))

    with pytest.raises(RuntimeError) as caught:
        _advance(solver, record, save_failure)
    assert caught.value is error
    metadata = json.loads(path.read_text())
    with np.load(path.with_suffix('.npz')) as saved:
        np.testing.assert_allclose(saved['times'], [0., .1, .2])
        np.testing.assert_allclose(saved['states'][:, 0], saved['times'])
        assert metadata['last_accepted_time'] == saved['times'][-1]
    assert metadata['incomplete'] is True
    assert metadata['error_message'] == str(error)


def test_output_failure_adds_context_without_masking_original():
    original = ValueError('original failure')

    def record():
        raise original

    def save_failure(error):
        raise OSError('disk full')

    with pytest.raises(ValueError) as caught:
        _advance(None, record, save_failure)
    assert caught.value is original
    assert 'disk full' in original.__notes__[0]


def test_no_history_initial_failure_has_nothing_to_save(tmp_path):
    path = tmp_path/'initial-failed.json'
    history = []

    def record():
        raise RuntimeError('initial endpoint unavailable')

    def save_failure(error):
        if history:
            path.write_text('must not be written')

    with pytest.raises(RuntimeError, match='initial endpoint unavailable'):
        _advance(None, record, save_failure)
    assert not path.exists()


@pytest.mark.parametrize('suffix', ['.json', '.npz'])
def test_existing_failure_artifacts_reject_before_any_run_and_are_preserved(tmp_path, suffix):
    output = tmp_path/'release-failed.json'
    existing = output.with_suffix(suffix)
    existing.write_bytes(b'previous accepted history')
    with pytest.raises(FileExistsError, match='failure artifact already exists'):
        _check_failure_output(output)
    assert existing.read_bytes() == b'previous accepted history'
