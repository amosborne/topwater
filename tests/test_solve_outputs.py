"""The solve pair stays coherent through preparation and publication failures."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from topwater import solve_outputs
from topwater.workflow_data import Result, Sea
from tests.test_workflow_data import result as tiny_result


def previous_pair(tmp_path):
    directory = tmp_path/'case'
    directory.mkdir()
    archive, plot = directory/'result.npz', directory/'result.png'
    archive.write_bytes(b'previous archive')
    plot.write_bytes(b'previous plot')
    return directory, archive, plot


def fake_plot(result, path):
    Path(path).write_bytes(b'new plot')


def assert_previous_pair(directory, archive, plot):
    assert archive.read_bytes() == b'previous archive'
    assert plot.read_bytes() == b'previous plot'
    assert {p.name for p in directory.iterdir()} == {'result.npz', 'result.png'}
    assert list(directory.parent.iterdir()) == [directory]


def test_success_publishes_real_png_and_matching_archive(tmp_path):
    import matplotlib.image as mpimg

    directory, archive, plot = previous_pair(tmp_path)
    data = tiny_result()
    position = data.position_m.copy()
    position[:, 0] += [0., .1, .3]
    position[:, 2] += [0., -.05, .12]
    data = replace(data, position_m=position)
    paths = solve_outputs.save_solve_outputs(data, directory)
    assert paths == (archive, plot)
    loaded = Result.load(archive)
    np.testing.assert_array_equal(loaded.time_s, data.time_s)
    np.testing.assert_array_equal(loaded.position_m, data.position_m)
    assert plot.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    assert mpimg.imread(plot).shape[0] > 100
    assert {p.name for p in directory.iterdir()} == {'result.npz', 'result.png'}
    assert list(tmp_path.iterdir()) == [directory]


def test_original_sea_is_explicit_plot_input_and_never_added_to_result_archive(tmp_path, monkeypatch):
    data = tiny_result()
    sea = Sea(np.array([.3]), np.array([1.]), np.array([.4]), np.array([.7]))
    recorded = []
    def capture(result, path, *, sea):
        recorded.append((result, sea))
        Path(path).write_bytes(b'complete combined plot')
    monkeypatch.setattr(solve_outputs, 'save_motion_plot', capture)
    archive, _ = solve_outputs.save_solve_outputs(data, tmp_path/'case', sea=sea)
    assert recorded == [(data, sea)]
    with np.load(archive, allow_pickle=False) as saved:
        assert set(saved.files) == set(Result.members)
        assert not any(name.startswith('sea_') for name in saved.files)


def test_plot_failure_preserves_both_previous_outputs(tmp_path, monkeypatch):
    directory, archive, plot = previous_pair(tmp_path)
    def fail(*args):
        raise ValueError('injected plot failure')
    monkeypatch.setattr(solve_outputs, 'save_motion_plot', fail)
    with pytest.raises(ValueError, match='plot failure'):
        solve_outputs.save_solve_outputs(tiny_result(), directory)
    assert_previous_pair(directory, archive, plot)


def test_archive_save_failure_preserves_both_previous_outputs(tmp_path, monkeypatch):
    directory, archive, plot = previous_pair(tmp_path)
    monkeypatch.setattr(solve_outputs, 'save_motion_plot', fake_plot)
    def fail(*args, **kwargs):
        raise OSError('injected archive failure')
    monkeypatch.setattr(Result, 'save', fail)
    with pytest.raises(OSError, match='archive failure'):
        solve_outputs.save_solve_outputs(tiny_result(), directory)
    assert_previous_pair(directory, archive, plot)


def test_second_replace_failure_rolls_back_previous_pair(tmp_path, monkeypatch):
    directory, archive, plot = previous_pair(tmp_path)
    monkeypatch.setattr(solve_outputs, 'save_motion_plot', fake_plot)
    original_replace = solve_outputs.os.replace
    def fail_second(source, destination):
        if Path(destination) == plot:
            raise PermissionError('injected second publication failure')
        return original_replace(source, destination)
    monkeypatch.setattr(solve_outputs.os, 'replace', fail_second)
    with pytest.raises(PermissionError, match='second publication'):
        solve_outputs.save_solve_outputs(tiny_result(), directory)
    assert_previous_pair(directory, archive, plot)


def test_second_replace_failure_removes_new_archive_if_no_previous_pair(tmp_path, monkeypatch):
    directory = tmp_path/'new-case'
    monkeypatch.setattr(solve_outputs, 'save_motion_plot', fake_plot)
    original_replace = solve_outputs.os.replace
    def fail_second(source, destination):
        if Path(destination) == directory/'result.png':
            raise PermissionError('injected publication failure')
        return original_replace(source, destination)
    monkeypatch.setattr(solve_outputs.os, 'replace', fail_second)
    with pytest.raises(PermissionError, match='publication failure'):
        solve_outputs.save_solve_outputs(tiny_result(), directory)
    assert list(directory.iterdir()) == []
    assert list(tmp_path.iterdir()) == [directory]


def test_backup_failure_is_before_publication(tmp_path, monkeypatch):
    directory, archive, plot = previous_pair(tmp_path)
    monkeypatch.setattr(solve_outputs, 'save_motion_plot', fake_plot)
    def fail(*args):
        raise OSError('injected backup failure')
    monkeypatch.setattr(solve_outputs.shutil, 'copyfile', fail)
    with pytest.raises(OSError, match='backup failure'):
        solve_outputs.save_solve_outputs(tiny_result(), directory)
    assert_previous_pair(directory, archive, plot)


def test_failed_os_rollback_retains_previous_bytes_outside_canonical_folder(tmp_path, monkeypatch):
    directory, archive, plot = previous_pair(tmp_path)
    monkeypatch.setattr(solve_outputs, 'save_motion_plot', fake_plot)
    original_replace = solve_outputs.os.replace
    def fail_publication_and_rollback(source, destination):
        if Path(destination) == plot or Path(source).name.startswith('previous-'):
            raise PermissionError('injected publication/rollback lock')
        return original_replace(source, destination)
    monkeypatch.setattr(solve_outputs.os, 'replace', fail_publication_and_rollback)
    with pytest.raises(RuntimeError, match='backups retained at'):
        solve_outputs.save_solve_outputs(tiny_result(), directory)
    recovery = list(tmp_path.glob('.topwater-solve-*'))
    assert len(recovery) == 1
    assert (recovery[0]/'previous-result.npz').read_bytes() == b'previous archive'
    assert (recovery[0]/'previous-result.png').read_bytes() == b'previous plot'
    assert plot.read_bytes() == b'previous plot'
    assert {p.name for p in directory.iterdir()} == {'result.npz', 'result.png'}
