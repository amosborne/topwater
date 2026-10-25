# CLI execution logs

Each command writes a log beside its data output:

| Command | Data output | Log |
|---|---|---|
| `generate-sea` | `sea.npz` | `sea.log` |
| `solve` | `result.npz` | `result.log` |
| `render` | `video.mp4` | `video.log` |

A future `generate-boat` command will use `boat.log`. The current example has
three NPZ files, one MP4 and these three text logs. No log is embedded in an NPZ.

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
ERROR. MATLAB output relayed by the backend is included. Native descriptor
writes, binary streams and unrelated logging handlers bound to other streams
are outside this capture. Failure to open the log stops the command before its
physical work starts.

Warnings need interpretation. [Warning triage](warning-triage.md) records the
current case's expected messages, evaluated approximations, outstanding accuracy
questions and defects requiring correction. An expected warning is not an
accuracy certificate or permission to suppress it.
