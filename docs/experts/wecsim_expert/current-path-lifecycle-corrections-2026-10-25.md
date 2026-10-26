# Current-path lifecycle corrections

Verified checkpoint, 2026-10-25. The read-only WEC audit identified two
executed-path defects. These changes do not modify physics, arrays, the fork
pin, body model, hydrodynamic producer or public result schema. No numerical or
static test was executed by the author. Independent source review passed;
root's serial Windows/R2025b checks passed 21 Python backend tests and 13 MATLAB
tests (two settings tests plus eleven finite-coefficient tests). Old-source
failures were reproduced. Evidence: `artifacts/workflow/current-path-checks-attempt2/`.
Exact processes/frozen identities were verified terminal and the lease released.

## Primary error plus engine-restoration failure

Previously `_run_matlab` raised a native-work error with buffered MATLAB output,
then a failure in its `finally` restoration replaced the displayed error.
Exception context existed internally, but CLI and `CommandLog` record the final
exception's string, hiding the original failure and transcript.

The combined error remains `EngineCleanupError`, so scratch retention and the
engine-contamination signal remain intact. Its message includes the primary
error, MATLAB stdout/stderr and cleanup failure. `primary_error`/`cleanup_error`
retain both exception objects, and explicit chaining preserves the cause.
Cleanup-only failure retains its original exception identity. Normal restoration
does not change the primary exception. A local captured failure is used rather
than ambient `sys.exc_info`, avoiding accidental capture of a caller's unrelated
handled exception. Interruptions still attempt restoration; if that restoration
also fails, the combined cleanup error includes the interruption and transcript.

Minimal executable pre/post witness, no MATLAB needed:

```powershell
uv run pytest -q tests/test_wecsim_backend.py -k matlab_cleanup_failure_keeps_primary_transcript_and_recovery
```

The fake engine emits `MATLAB_PRIMARY_MARKER`, `MATLAB_STDOUT_MARKER` and
`MATLAB_STDERR_MARKER`; restoration raises `ENGINE_CLEANUP_MARKER`. Before the
fix, the primary markers are missing from the displayed/logged error. After
the fix, all markers must remain in `str(error)` and the actual failed command
log, with `EngineCleanupError`, retained scratch and both exception objects.
Controls cover cleanup after successful native work and invocation inside a
caller's unrelated exception handler. The old implementation failed both combined-
failure cases; the corrected implementation passes all four parameterized cases.

## Newer-model load policy

Pinned `initializeWecSim.m:541` sets root `ErrorIfLoadNewModel` to `off`.
`WorkflowEngineGuard` previously omitted it from its snapshot. The guard now
captures its exact value, independently restores it on every cleanup attempt,
and reads it back before reporting success. Other global/model settings are
unchanged. Its existing destructor fallback still restores the setting after a
forced error, without closing any unowned model or globally clearing classes.

Minimal executable pre/post witness in an isolated MATLAB environment:

```matlab
results = runtests('tests/matlab/workflowEngineSettingTest.m');
assertSuccess(results);
```

For each initial `on`/`off` value, the test changes the setting to its opposite,
then checks restoration after normal cleanup and after an intentional exception.
Before the fix the opposite value remains. After the fix the original value
must return. Test teardown restores the incoming setting/path; no model,
simulation or named engine is required. This tests the guard boundary rather
than rerunning WEC preprocessing. Both tests fail against the old guard and pass
against the corrected guard on Windows/R2025b.

## Integration and limits

The backend resource-copy tuple also includes the separately owned
`validate_workflow_hydrodynamics.m` helper. That is a packaging integration for
the BEMIO owner's finite-boundary check; no edits were made here to the owner's
helper or `run_workflow_case.m`.

Root owns independent review, serial execution, evidence and checkpointing.
Broader session workspace ownership, physics accuracy, inertia export/native
probes, public resimulation and upstream work are outside this correction.
