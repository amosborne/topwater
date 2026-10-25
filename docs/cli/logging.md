# CLI execution logs

Each command writes a log beside its data output:

| Command | Data output | Log |
|---|---|---|
| `generate-sea` | `sea.npz` | `sea.log` |
| `solve` | `result.npz`, `result.png` | `result.log` |
| `render` | `video.mp4` | `video.log` |

A future `generate-boat` command will use `boat.log`. The current example has
three NPZ files, one PNG, one MP4 and these three text logs. No log is embedded in an NPZ.

Logs describe the latest invocation and are replaced when that command starts.
A failed attempt therefore replaces its log while preserving the previous data
output until a complete replacement is validated. Other commands' logs remain
untouched. Argument parsing and help occur before a log directory is resolved;
those messages appear only in the terminal.

Every nonempty line has the same UTF-8 format:

```text
2026-10-25T19:00:00.123Z INFO solve BEGIN directory=examples/simple-boat duration=60.0 solver=wecsim
2026-10-25T19:00:00.124Z INFO solve STAGE validate solve inputs and backend
2026-10-25T19:01:00.000Z WARNING solve warning [category]: diagnostic message
2026-10-25T19:02:00.000Z INFO solve PUBLISHED examples/simple-boat/result.npz
2026-10-25T19:02:00.001Z INFO solve END status=success elapsed_s=120.001000 output_published=true warnings=1
```

Timestamps use UTC. Stages, elapsed time, publication status and warning records
help distinguish a completed output from a failed or interrupted attempt. The
warning count counts log records, including warning summaries; affected time
samples and BEM solves are reported separately by the producing diagnostic.
A missing `END` record indicates an incomplete log, not successful completion.

Python stdout/stderr are copied to the terminal unchanged and recorded with
levels; ANSI escapes are removed from the file. Ordinary stderr progress remains
INFO. Existing warning prefixes are classified as WARNING, and exceptions as
ERROR. Console flushes do not end log records: Rich may flush after each styled
word or punctuation mark. Only newline/carriage-return boundaries end records;
an unterminated final line is recorded when the command exits, including on an
exception. A partially written progress line can therefore be visible in the
terminal before it appears in the file.

Rich time/level columns are recognized. Immediately following rows starting at
the same message column inherit a warning/error level on that stream. A blank
row, different indentation or a new level header resets that continuation.
This is a rendered-text heuristic, not interception of the original LogRecord:
other layouts can lose continuation severity, and unrelated text at exactly
the same indentation can inherit it. Existing logging handlers and their
formatting are unchanged. The warning count counts rendered WARNING rows,
including wrapped continuation rows, rather than unique warning events.

MATLAB output relayed by the backend is included. Native descriptor
writes, binary streams and unrelated logging handlers bound to other streams
are outside this capture. Failure to open the log stops the command before its
physical work starts.

Warnings need interpretation. [Warning triage](warning-triage.md) records the
current case's expected messages, evaluated approximations, outstanding accuracy
questions and defects requiring correction. An expected warning is not an
accuracy certificate or permission to suppress it.
