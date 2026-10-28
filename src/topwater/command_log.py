"""Latest-attempt command logs, with Python console text retained unchanged.

Use one context per serial CLI invocation. The context captures writes through
the current Python stdout/stderr, including backend text printed by Python.
Native file-descriptor writes, subprocess output not relayed through Python,
binary stream writes, and logging handlers bound to other streams are outside
this capture. Rich-style continuation levels are inferred from rendered text;
the original LogRecord is not intercepted. This is an execution log, not a
physical-accuracy assessment.
"""

from datetime import datetime, timezone
import logging
from pathlib import Path
import re
import sys
import threading
import time


_ANSI = re.compile(r'\x1b(?:\[[0-?]*[ -/]*[@-~]|\][^\x07\x1b]*(?:\x07|\x1b\\)|[@-_])')
_WARNING = re.compile(
    r'^\s*(?:warning(?:\s+summary)?\s*(?:\[|:)|'
    r'WEC-Sim applicability:|Surface numerical approximation:)', re.IGNORECASE)
_ERROR = re.compile(r'^\s*(?:ERROR|CRITICAL|FATAL)\b', re.IGNORECASE)
# Capytaine's Rich handler renders a time column (or spaces when repeated),
# a level column, and then the message. Match that header before inferring the
# severity of immediately following rows at the identical message column.
_RICH_HEADER = re.compile(
    r'^(?:\[[^\]\r\n]*\d:\d[^\]\r\n]*\]| {2,})[ \t]+'
    r'(?P<level>DEBUG|INFO|WARNING|ERROR|CRITICAL)[ \t]+(?=\S)')


class _UTCFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        return datetime.fromtimestamp(record.created, timezone.utc).isoformat(
            timespec='milliseconds').replace('+00:00', 'Z')


class _FileHandler(logging.FileHandler):
    def handleError(self, record):
        # Default logging error reporting writes to stderr, which would recurse
        # through the tee. A log-write failure must be visible to the caller.
        failure = sys.exc_info()[1]
        if failure is not None:
            raise failure
        raise OSError('command log could not write a record')


class _ConsoleTee:
    def __init__(self, owner, console):
        self.owner = owner
        self.console = console
        self.pending = ''
        self.previous_cr = False
        self.continuation_level = None
        self.continuation_column = None

    def write(self, text):
        if not isinstance(text, str):
            raise TypeError('console text must be a string')
        with self.owner._lock:
            self.console.write(text)
            self.pending += text
            while True:
                match = re.search(r'[\r\n]', self.pending)
                if match is None:
                    break
                line = self.pending[:match.start()]
                self.pending = self.pending[match.end():]
                if match.group() == '\n' and self.previous_cr and not line:
                    self.previous_cr = False
                    continue
                self.previous_cr = match.group() == '\r'
                self.owner._console_line(line, self)
        return len(text)

    def flush(self):
        with self.owner._lock:
            self.console.flush()
            self.owner._handler.flush()

    def drain(self):
        """Publish an unterminated tail only when the command ends."""
        with self.owner._lock:
            if self.pending:
                line, self.pending = self.pending, ''
                self.owner._console_line(line, self)

    def __getattr__(self, name):
        return getattr(self.console, name)


class CommandLog:
    """Tee console text into a timestamped UTF-8 log for the latest attempt.

    Opening the context truncates its log and creates its parent directory;
    failures opening it occur before either Python stream is replaced.
    ``event`` writes only to the file. Call ``published`` only after the caller
    has validated and published an output. Warning counts count WARNING records,
    including repeated-warning summaries, rather than distinct physical events.
    Exceptions are logged and propagate. An ERROR record also marks a handled
    failure when a command catches its exception inside this context.
    Console flushes do not delimit log lines; unterminated tails are recorded
    on exit. Rich warning/error continuation levels use a conservative rendered
    column heuristic, not interception or changes to the caller's log handlers.
    """

    def __init__(self, path, command, details=''):
        if not isinstance(command, str) or not command.strip() or '\n' in command or '\r' in command:
            raise ValueError('command must be a nonempty single-line label')
        self.path = Path(path)
        self.command = command.strip()
        self.details = str(details)
        self.warning_count = 0
        self._errors = 0
        self._outputs = []
        self._active = False
        self._used = False
        self._lock = threading.RLock()

    def __enter__(self):
        if self._used:
            raise RuntimeError('a CommandLog context cannot be reused')
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handler = _FileHandler(self.path, mode='w', encoding='utf-8')
        handler.setFormatter(_UTCFormatter('%(asctime)s %(levelname)s %(command)s %(message)s'))
        self._logger = logging.Logger('topwater.command', level=logging.DEBUG)
        self._logger.propagate = False
        self._logger.addHandler(handler)
        self._handler = handler
        self._stdout, self._stderr = sys.stdout, sys.stderr
        self._tees = (_ConsoleTee(self, self._stdout), _ConsoleTee(self, self._stderr))
        self._started = time.monotonic()
        self._active = self._used = True
        try:
            self.event('INFO', 'BEGIN' + (' ' + self.details if self.details else ''))
        except BaseException:
            self._active = False
            handler.close()
            raise
        sys.stdout, sys.stderr = self._tees
        return self

    def event(self, level, message):
        """Write file-only stage information using a standard logging level."""
        if not self._active:
            raise RuntimeError('command log is not active')
        if isinstance(level, str):
            levels = {'DEBUG': logging.DEBUG, 'INFO': logging.INFO,
                      'WARNING': logging.WARNING, 'WARN': logging.WARNING,
                      'ERROR': logging.ERROR, 'CRITICAL': logging.CRITICAL}
            try:
                level = levels[level.upper()]
            except KeyError as error:
                raise ValueError('unsupported command-log level') from error
        if level not in (logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL):
            raise ValueError('unsupported command-log level')
        with self._lock:
            for line in _ANSI.sub('', str(message)).splitlines():
                if not line.strip():
                    continue
                self._logger.log(level, line, extra={'command': self.command})
                if level == logging.WARNING:
                    self.warning_count += 1
                elif level >= logging.ERROR:
                    self._errors += 1

    def published(self, path):
        """Record an output whose validation/publication the caller completed."""
        self._outputs.append(str(Path(path)))
        self.event('INFO', f'PUBLISHED {self._outputs[-1]}')

    def _console_line(self, line, tee):
        line = _ANSI.sub('', line)
        if not line.strip():
            tee.continuation_level = tee.continuation_column = None
            return
        header = _RICH_HEADER.match(line)
        column = tee.continuation_column
        if header is not None:
            level = header.group('level')
            if level in ('WARNING', 'ERROR', 'CRITICAL'):
                tee.continuation_level = level
                tee.continuation_column = header.end()
            else:
                tee.continuation_level = tee.continuation_column = None
        elif _ERROR.search(line) or line.lstrip().startswith(f'topwater {self.command}:'):
            level = 'ERROR'
            tee.continuation_level = tee.continuation_column = None
        elif (column is not None and line.startswith(' '*column)
              and len(line) > column and not line[column].isspace()):
            level = tee.continuation_level
        elif (_WARNING.search(line) or re.search(r'\b(?:WARNING|WARN)\b', line)
              or re.search(r':\d+:\s+(?:[A-Za-z]+Warning|Warning):', line)
              or (line.lstrip().startswith('Capytaine:') and re.search(r'\bwarning\b', line, re.IGNORECASE))):
            level = 'WARNING'
            tee.continuation_level = tee.continuation_column = None
        else:
            level = 'INFO'
            tee.continuation_level = tee.continuation_column = None
        self.event(level, line)

    def __exit__(self, kind, error, traceback):
        cleanup_failure = None
        try:
            for tee in self._tees:
                tee.flush()
                tee.drain()
            if error is not None:
                self.event('ERROR', f'{type(error).__name__}: {error}'.rstrip(': '))
            status = ('interrupted' if isinstance(error, KeyboardInterrupt)
                      else 'failed' if error is not None or self._errors else 'success')
            self.event('ERROR' if status != 'success' else 'INFO',
                f'END status={status} elapsed_s={time.monotonic()-self._started:.6f} '
                f'output_published={str(bool(self._outputs)).lower()} warnings={self.warning_count}')
        except BaseException as failure:
            cleanup_failure = failure
        finally:
            sys.stdout, sys.stderr = self._stdout, self._stderr
            self._active = False
            try:
                self._handler.close()
            except BaseException as failure:
                if cleanup_failure is None:
                    cleanup_failure = failure
        if cleanup_failure is not None:
            if error is None:
                raise cleanup_failure
            message = f'Command log finalization failed: {cleanup_failure}'
            error.add_note(message)
            try:
                self._stderr.write(message+'\n')
                self._stderr.flush()
            except Exception:
                pass
        return False
