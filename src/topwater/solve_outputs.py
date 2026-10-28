"""Publish the final solve archive and six-DOF plot as a recoverable pair."""

import os
from pathlib import Path
import shutil
import tempfile

from .result_plot import save_motion_plot


def save_solve_outputs(result, directory, *, sea=None):
    """Prepare both files before publishing; roll back handled replace failures.

    Return (result.npz path, result.png path). Same-volume scratch and backups
    stay outside the canonical directory. Backups use filesystem copies, never
    copies of large result arrays. The two replaces are individually atomic;
    abrupt process termination between them is not globally atomic. If the OS
    also prevents rollback, retain the backups and report their recovery path.
    The optional original Sea adds its static spectrum to result.png only;
    it is never stored in result.npz. Production CLI solve supplies this input.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.topwater-solve-', dir=directory.parent))
    targets = (directory/'result.npz', directory/'result.png')
    ready = (staging/'result.npz', staging/'result.png')
    retain_backups = False
    try:
        if os.stat(staging).st_dev != os.stat(directory).st_dev:
            raise OSError('solve outputs require staging on the destination filesystem')
        if sea is None:
            save_motion_plot(result, ready[1])
        else:
            save_motion_plot(result, ready[1], sea=sea)
        result.save(ready[0])
        if not all(path.is_file() and path.stat().st_size > 0 for path in ready):
            raise OSError('both solve outputs must be complete files before publication')
        backups = {}
        for target in targets:
            if target.exists():
                backup = staging/('previous-'+target.name)
                shutil.copyfile(target, backup)
                backups[target] = backup
        published = []
        try:
            for source, target in zip(ready, targets):
                os.replace(source, target)
                published.append(target)
        except BaseException as publication_error:
            rollback_errors = []
            for target in reversed(published):
                try:
                    if target in backups:
                        os.replace(backups[target], target)
                    else:
                        target.unlink()
                except OSError as rollback_error:
                    rollback_errors.append(rollback_error)
            if rollback_errors:
                retain_backups = True
                raise RuntimeError(f'could not restore previous solve outputs; backups retained at {staging}') from publication_error
            raise
        return targets
    finally:
        if not retain_backups:
            shutil.rmtree(staging, ignore_errors=True)
