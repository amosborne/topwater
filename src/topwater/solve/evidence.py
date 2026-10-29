"""Foundational completed-field error evidence, independent of replay and analysis."""
from pathlib import Path
import json
import sys
import traceback
import numpy as np


class SurfaceExportError(RuntimeError):
    """Primary geometry/metric failure with a private replayable field cache."""
    def __init__(self, message, recovery_directory, *, stage, frame_index, time_s):
        super().__init__(message)
        self.recovery_directory = Path(recovery_directory)
        self.stage = stage
        self.frame_index = int(frame_index)
        self.time_s = None if time_s is None else float(time_s)


def _save_completed_cache(directory, data):
    """Private disposable intermediate, never placed in the example folder."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    path = directory/'completed-field.npz'
    with path.open('xb') as stream:
        np.savez(stream, **data)
    return path


def _save_frame_failure(directory, evidence, error, stage, index, stamp, cache_path):
    """Best-effort secondary diagnostics must never replace the primary error."""
    directory = Path(directory)
    report = dict(stage=stage, frame_index=int(index), time_s=(None if stamp is None else float(stamp)),
                  error_type=type(error).__name__, error_message=str(error),
                  traceback=traceback.format_exc(),
                  completed_cache=str(cache_path))
    try:
        directory.mkdir(parents=True, exist_ok=True)
        with (directory/'failure.json').open('x', encoding='utf-8') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
        with (directory/'failure-frame.npz').open('xb') as stream:
            np.savez(stream, **evidence)
    except Exception as diagnostic_error:
        print(f'Surface failure evidence could not be fully saved: {diagnostic_error}; '
              'the original failure remains authoritative.', file=sys.stderr)
