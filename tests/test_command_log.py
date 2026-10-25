"""Real-file command logging controls; no solver or encoder execution."""

from io import StringIO
import re
import sys

import pytest

from topwater.command_log import CommandLog


@pytest.fixture
def consoles(monkeypatch):
    # Install during the test call: pytest replaces capture streams between
    # fixture setup and the call, so setup-time replacements are not observed
    # by CommandLog under ordinary pytest capture.
    def install():
        stdout, stderr = StringIO(), StringIO()
        monkeypatch.setattr(sys, 'stdout', stdout)
        monkeypatch.setattr(sys, 'stderr', stderr)
        return stdout, stderr
    return install


def records(path):
    lines = path.read_text(encoding='utf-8').splitlines()
    for line in lines:
        assert re.fullmatch(r'\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3}Z '
                            r'(?:DEBUG|INFO|WARNING|ERROR|CRITICAL) solve .+', line), line
    return lines


def test_tee_preserves_console_once_and_formats_every_file_line(tmp_path, consoles):
    consoles = consoles()
    stdout, stderr = consoles
    path = tmp_path/'example'/'result.log'
    with CommandLog(path, 'solve', details='duration=60 solver=wecsim') as log:
        print('ordinary progress')
        print('stderr progress', file=sys.stderr)
        log.event('INFO', 'preparing hydrodynamics')
        log.published(tmp_path/'example'/'result.npz')
    assert sys.stdout is stdout and sys.stderr is stderr
    assert stdout.getvalue() == 'ordinary progress\n'
    assert stderr.getvalue() == 'stderr progress\n'
    lines = records(path)
    assert any(' INFO solve ordinary progress' in line for line in lines)
    assert any(' INFO solve stderr progress' in line for line in lines)
    assert any(' INFO solve preparing hydrodynamics' in line for line in lines)
    assert 'status=success' in lines[-1] and 'output_published=true' in lines[-1]
    assert 'warnings=0' in lines[-1]


@pytest.mark.parametrize('message', [
    'warning [calm-restoring-force] at 1s: departure',
    'warning summary [calm-restoring-force]: 3 samples',
    'WARNING Capytaine mesh resolution',
    'Capytaine: Mesh resolution warning applies to 8 solves',
    'WEC-Sim applicability: wet-set exceeds 10%',
    'Surface numerical approximation: provisional field',
    '\x1b[33mWarning: native model approximation\x1b[0m',
])
def test_warning_markers_have_consistent_severity(tmp_path, consoles, message):
    consoles = consoles()
    path = tmp_path/'result.log'
    with CommandLog(path, 'solve') as log:
        print(message, file=sys.stderr)
    lines = records(path)
    warnings = [line for line in lines if ' WARNING solve ' in line]
    assert len(warnings) == 1
    assert '\x1b' not in path.read_text(encoding='utf-8')
    assert log.warning_count == 1
    assert 'warnings=1' in lines[-1] and 'status=success' in lines[-1]
    assert consoles[1].getvalue() == message+'\n'


def test_partial_lines_flush_and_streams_restore_on_exception(tmp_path, consoles):
    consoles = consoles()
    path = tmp_path/'result.log'
    previous_output = tmp_path/'result.npz'
    previous_output.write_bytes(b'previous result')
    primary = ValueError('invalid contact topology')
    with pytest.raises(ValueError) as caught:
        with CommandLog(path, 'solve'):
            sys.stdout.write('partial progress')
            sys.stderr.write('warning [contact]: partial warning')
            raise primary
    assert caught.value is primary
    assert sys.stdout is consoles[0] and sys.stderr is consoles[1]
    assert previous_output.read_bytes() == b'previous result'
    lines = records(path)
    assert any(' INFO solve partial progress' in line for line in lines)
    assert any(' WARNING solve warning [contact]: partial warning' in line for line in lines)
    assert any(' ERROR solve ValueError: invalid contact topology' in line for line in lines)
    assert 'status=failed' in lines[-1] and 'output_published=false' in lines[-1]


def test_handled_cli_fatal_is_not_mislabelled_success(tmp_path, consoles):
    consoles = consoles()
    path = tmp_path/'result.log'
    with CommandLog(path, 'solve'):
        print('topwater solve: native coupled dynamics are unavailable', file=sys.stderr)
    lines = records(path)
    assert any(' ERROR solve topwater solve:' in line for line in lines)
    assert 'status=failed' in lines[-1]


def test_interrupt_is_logged_and_not_suppressed(tmp_path, consoles):
    consoles = consoles()
    path = tmp_path/'result.log'
    with pytest.raises(KeyboardInterrupt):
        with CommandLog(path, 'solve'):
            raise KeyboardInterrupt
    assert sys.stdout is consoles[0] and sys.stderr is consoles[1]
    assert 'status=interrupted' in records(path)[-1]


def test_latest_attempt_truncates_previous_log(tmp_path, consoles):
    consoles = consoles()
    path = tmp_path/'result.log'
    with CommandLog(path, 'solve'):
        print('old attempt')
    with CommandLog(path, 'solve'):
        print('new attempt')
    text = path.read_text(encoding='utf-8')
    assert 'old attempt' not in text and 'new attempt' in text
    assert text.count('BEGIN') == 1 and text.count('END status=') == 1


def test_explicit_flush_emits_partial_line_and_ansi_is_removed(tmp_path, consoles):
    consoles = consoles()
    path = tmp_path/'result.log'
    with CommandLog(path, 'solve'):
        sys.stdout.write('\x1b[32mprogress\x1b[0m')
        sys.stdout.flush()
        assert ' INFO solve progress' in path.read_text(encoding='utf-8')
    assert consoles[0].getvalue() == '\x1b[32mprogress\x1b[0m'
    assert '\x1b' not in path.read_text(encoding='utf-8')


def test_open_failure_leaves_streams_and_output_untouched(tmp_path, consoles):
    consoles = consoles()
    path = tmp_path/'result.log'
    path.mkdir()
    previous_output = tmp_path/'result.npz'
    previous_output.write_bytes(b'previous result')
    with pytest.raises(OSError):
        with CommandLog(path, 'solve'):
            pytest.fail('physics must not begin when the log cannot open')
    assert sys.stdout is consoles[0] and sys.stderr is consoles[1]
    assert previous_output.read_bytes() == b'previous result'
