"""Real-file command logging controls; no solver or encoder execution."""

from io import StringIO
import logging
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
    'warning [solver] at 1s: numerical issue',
    'warning summary [solver]: 3 samples',
    'WARNING Capytaine mesh resolution',
    'Capytaine: Mesh resolution warning applies to 8 solves',
    'Surface numerical approximation: provisional field',
    '\x1b[33mWarning: native model approximation\x1b[0m',
    'C:\\case\\result_plot.py:46: UserWarning: Six-motion Euler labels are ambiguous',
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


def test_explicit_flush_preserves_partial_line_until_exit_and_strips_ansi(tmp_path, consoles):
    consoles = consoles()
    path = tmp_path/'result.log'
    with CommandLog(path, 'solve'):
        sys.stdout.write('\x1b[32mprogress\x1b[0m')
        sys.stdout.flush()
        assert ' INFO solve progress' not in path.read_text(encoding='utf-8')
        assert consoles[0].getvalue() == '\x1b[32mprogress\x1b[0m'
    assert ' INFO solve progress' in path.read_text(encoding='utf-8')
    assert consoles[0].getvalue() == '\x1b[32mprogress\x1b[0m'
    assert '\x1b' not in path.read_text(encoding='utf-8')


def test_write_flush_segments_form_one_real_warning_line(tmp_path, consoles):
    consoles = consoles()
    path = tmp_path/'result.log'
    fragments = ['[1:54:49 PM]', ' ', '\x1b[33mWARNING', '\x1b[0m ',
                 'Mesh resolution for ', 'RadiationProblem', '(', 'body', '=', 'boat', ')']
    with CommandLog(path, 'solve') as log:
        for fragment in fragments:
            sys.stdout.write(fragment)
            sys.stdout.flush()
        assert len(records(path)) == 1  # BEGIN only; flush is not a newline.
        sys.stdout.write('\r')
        sys.stdout.flush()
        sys.stdout.write('\n')  # Split Windows newline must not create a blank row.
        sys.stdout.flush()
    assert consoles[0].getvalue() == ''.join(fragments)+'\r\n'
    lines = records(path)
    assert len(lines) == 3
    assert lines[1].endswith('WARNING solve [1:54:49 PM] WARNING Mesh resolution for RadiationProblem(body=boat)')
    assert log.warning_count == 1


def test_rich_continuation_levels_reset_and_do_not_cross_streams(tmp_path, consoles):
    consoles = consoles()
    path = tmp_path/'result.log'
    prefix = '[13:54:49] WARNING  '
    padding = ' '*len(prefix)
    with CommandLog(path, 'solve'):
        print(prefix+'first warning row')
        print(padding+'wrapped warning row')
        print(padding+'other stream progress', file=sys.stderr)
        print('progress resumed')
        print(padding+'indented progress after reset')
        print(prefix+'another warning')
        print()
        print(padding+'indented progress after blank')
        print(prefix+'third warning')
        print('[13:54:50] INFO     explicit progress')
        print(padding+'wrapped ordinary progress')
    lines = records(path)
    assert any(' WARNING solve '+padding+'wrapped warning row' in line for line in lines)
    for message in ('other stream progress', 'indented progress after reset',
                    'indented progress after blank', 'wrapped ordinary progress'):
        assert any(' INFO solve '+padding+message in line for line in lines)


def test_actual_rich_handler_remains_unchanged_with_segment_flushes(tmp_path, consoles):
    rich_console = pytest.importorskip('rich.console')
    rich_logging = pytest.importorskip('rich.logging')
    consoles = consoles()
    path = tmp_path/'result.log'

    class SegmentFlushingStream:
        # Exercise the per-styled-segment flush pattern of Rich's Windows
        # renderer while retaining the actual RichHandler layout on every OS.
        def __init__(self, target):
            self.target = target

        def write(self, text):
            for part in re.findall(r'\s+|\S+', text):
                self.target.write(part)
                self.target.flush()
            return len(text)

        def flush(self):
            self.target.flush()

    def make_handler(stream):
        console = rich_console.Console(file=stream, width=58, force_terminal=False,
                                       legacy_windows=False, color_system=None)
        handler = rich_logging.RichHandler(console=console, show_path=False,
                                           log_time_format='[%X]')
        handler.setFormatter(logging.Formatter('%(message)s'))
        handler.addFilter(logging.Filter())
        return handler

    warning = logging.LogRecord('capytaine.test', logging.WARNING, __file__, 1,
        'Mesh resolution may be insufficient near the hull; check the largest '
        'panel radius before accepting this field.\nSecond warning paragraph '
        'continues the same diagnostic without repeating its level.', (), None)
    progress = logging.LogRecord('capytaine.test', logging.INFO, __file__, 2,
                                 'Coefficient generation continues', (), None)
    reference = StringIO()
    reference_handler = make_handler(reference)
    reference_handler.handle(warning)
    warning_rows = [line for line in reference.getvalue().splitlines() if line.strip()]
    reference_handler.handle(progress)

    with CommandLog(path, 'solve') as log:
        handler = make_handler(SegmentFlushingStream(sys.stdout))
        logger = logging.Logger('capytaine.test', logging.DEBUG)
        logger.addHandler(handler)
        original = (tuple(logger.handlers), tuple(handler.filters), handler.formatter,
                    handler.level, logger.level, logger.propagate)
        logger.handle(warning)
        logger.handle(progress)
        assert original == (tuple(logger.handlers), tuple(handler.filters), handler.formatter,
                            handler.level, logger.level, logger.propagate)
    assert consoles[0].getvalue() == reference.getvalue()
    assert consoles[1].getvalue() == ''
    lines = records(path)
    assert [line.split(' solve ', 1)[1] for line in lines[1:-1]] == [
        line for line in reference.getvalue().splitlines() if line.strip()]
    assert len(warning_rows) > 2  # Actual wrapped/multiline Rich output.
    assert all(' WARNING solve ' in line for line in lines[1:1+len(warning_rows)])
    assert all(' INFO solve ' in line for line in lines[1+len(warning_rows):-1])
    assert log.warning_count == len(warning_rows)


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
