"""Small real sparse-system controls; no canoe assembly or physics solve."""
import json

import numpy as np
import pytest
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import splu

from .experiments.free_surface import p4_fluid_ordering as module


def saved_system(path):
    line = diags([-np.ones(6), 2*np.ones(7), -np.ones(6)], [-1, 0, 1], format='csc')
    matrix = (kron(eye(7), line)+kron(line, eye(7))+.2*eye(49)).tocsc()
    rows = np.arange(1., 50.)[:, None]
    known = np.sin(rows*np.arange(1., 7.)[None, :]/11)
    load = matrix@known
    baseline = splu(matrix, permc_spec='COLAMD').solve(load)
    velocity = np.arange(1., 7.)/100
    rigid = np.diag(np.arange(2., 8.))
    potential = baseline@velocity
    momentum = rigid@velocity+1000*load.T@potential
    coupled = dict(rigid_mass=rigid, body_map=np.eye(6), body_momentum=momentum,
        psi=np.zeros(3), fluid_surface_rhs=np.zeros(49), p_effective=momentum,
        reference_velocity=velocity, reference_potential=potential)
    module.save_input(path, matrix, load, baseline, 1000*load.T@baseline,
        np.zeros(9), np.arange(49), np.arange(49, 52),
        dict(case='49-row real 2D elliptic control; not a canoe qualification'), coupled=coupled)
    return known, load


def test_saved_exact_matrix_mmd_recovers_known_six_column_solution(tmp_path):
    path, output = tmp_path/'input.json', tmp_path/'comparison.json'
    known, load = saved_system(path)
    report = module.compare(path, output)
    assert report['complete'] and not report['qualified']
    assert report['relative_residual'] < 3e-14
    with np.load(output.with_suffix('.npz')) as saved:
        np.testing.assert_allclose(saved['response'], known, atol=3e-13, rtol=3e-13)
        np.testing.assert_allclose(saved['added_mass'], 1000*load.T@known, atol=1e-9, rtol=3e-13)
    assert report['response_relative_difference'] < 3e-14
    assert report['factor_wall_seconds'] > 0 and report['peak_rss_bytes_after'] > 0
    with pytest.raises(FileExistsError):
        module.compare(path, output)


def test_mutated_saved_matrix_rejected_before_factor(tmp_path, monkeypatch):
    path = tmp_path/'input.json'
    saved_system(path)
    with path.with_suffix('.npz').open('ab') as stream:
        stream.write(b'changed archive')
    monkeypatch.setattr(module, 'splu', lambda *a, **kw: pytest.fail('changed input reached factorization'))
    with pytest.raises(ValueError, match='differs'):
        module.compare(path, tmp_path/'changed.json')


def test_factor_failure_preserves_original_error_and_report(tmp_path, monkeypatch):
    path, output = tmp_path/'input.json', tmp_path/'failed.json'
    saved_system(path)
    error = MemoryError('factor memory exhausted')
    def failed(*args, **kwargs):
        raise error
    monkeypatch.setattr(module, 'splu', failed)
    with pytest.raises(MemoryError) as caught:
        module.compare(path, output)
    assert caught.value is error
    report = json.loads(output.read_text())
    assert not report['complete'] and report['failure']['type'] == 'MemoryError'
    assert report['factor_wall_seconds'] >= 0 and report['peak_rss_bytes_after'] > 0
