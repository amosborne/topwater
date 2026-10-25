# Shared WEC-Sim engine lifecycle verification

Status: **PASS for shared-engine lifecycle only**, 2026-10-25. This is not a
body-dynamics, total-surface, accuracy, or performance qualification.

## Frozen execution and evidence

Coordinator `/root` ran the serial checkpoint at
`2026-10-25T17:43:55.068695+00:00`. Exact producer PID 34608, owner launcher PID
24052, terminal exit 0, complete and frozen-input identities verified. Root
completed terminal reconciliation and lease release and reported independent
reviewer approval before this documentation update.

- [Frozen manifest](../../../artifacts/workflow/engine-guard-attempt4/manifest.json),
  SHA256 `cf7f247fd9e28ea1d283d0cdeb8a049d38a0ae99df9a8011b15db2ab8d3f1569`.
- [Combined output](../../../artifacts/workflow/engine-guard-attempt4/output.log),
  SHA256 `0d9289cd40e2ba57e394df9a5ef23ee6e300a7cffc1f7c4abaeefbf0f49e0b11`.
- [Terminal receipt](../../../artifacts/workflow/engine-guard-attempt4/terminal.json),
  SHA256 `e0d65f2833b085ac4baa9cb94f0a5ba13d232a8109fd0a95099a27e057d37842`.

The manifest owns exact source/dependency hashes, including the packaged
`WorkflowEngineGuard.m`, Python wrapper, tests, private probe, MATLAB engine
package, and MATLAB executable. Earlier failed attempts remain preserved under
their original outcomes.

## Checks that passed

- Thirteen backend tests passed in 0.86 s, including cleanup failures in settings
  restoration, JSON decoding, guard destruction, and final search-path restoration.
- The MATLAB skill's `m_check.py` reported all four packaged `.m` files clean.
- MATLAB `checkcode` found no parse/syntax rejection. It returned seven other
  diagnostics for `run_workflow_case.m`, two for `WorkflowEngineGuard.m`, and zero
  for the other two files; these counts are not represented as zero warnings.
- The private probe performed two cycles using distinct temporary class
  directories. One was normal; the other raised an intentional MATLAB error.
  Each created a uniquely owned empty model and changed diary, file-generation,
  warning, thread, figure-default, and directory state without simulating dynamics.
- After each cycle, exact Python-side `pwd` and search-path strings matched the
  originals. MATLAB assertions checked diary filename/on-off state, warnings,
  thread count, figure default, and cache/codegen folder paths. The owned model
  was closed and the temporary class directory could be deleted.
- The probe attached to the existing hidden/minimized `CAPYTAINE_REFERENCE`
  engine. It started no engine, opened no model editor or Mechanics Explorer,
  and did not restore/activate the desktop.

## Ownership and failure behavior

The guard is registered before preprocessing or model construction. It restores
non-path settings first while retaining its class folder. Python then explicitly
deletes/releases the guard, clears only that class, and restores/verifies the
original search path last. File-generation restoration can remove the cache
folder from the path, so the guard explicitly keeps its class reachable through
destruction.

Cleanup operations are independent: a settings/JSON/destructor error cannot skip
the later release and final path attempts. Unverified reasons are aggregated in
`EngineCleanupError`. Scratch is retained in system temporary storage for recovery
if restoration cannot be verified; an existing public result stays intact.

## Limits and next checkpoint

The probe exercises lifecycle and parser behavior, not BEMIO import, compilation
of the full body model, time integration, or total free-surface reconstruction.
The file-generation folder-structure enum is copied/restored by source but was
not separately asserted by this probe. No convergence, energy, passivity, or
nonlinear moving-contact claim follows from this pass.

The Capytaine owner continues the independently reviewed total-field work.
Production body/surface execution still requires separate coordinator allocation
and its physical acceptance checks. This follow-up changed documentation only.
