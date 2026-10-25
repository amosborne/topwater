# BEMIO sidecar handles and scratch-error preservation

Status: **PASS for importer handle lifecycle and numeric-array identity**.
Maintained-fork commit `f07a8dfe1d81d9430a22c299a83d91252afb3a33` is published on
the user-owned `origin/topwater/reference-v7.1.0`; the coordinator updated the
backend and fork pin. No coefficients, body loads, or physics policies change.

## Executed checkpoint

[Regression report](../../../artifacts/workflow/bemio-handles-attempt3/report.json),
SHA256 `437dacee9c4a77d73df88884653850ef82e9b9148397199ffae9c44380207564`,
records 1.672 s reported wall time, with exact source/dependencies in the adjacent
[manifest](../../../artifacts/workflow/bemio-handles-attempt3/manifest.json) and
[terminal receipt](../../../artifacts/workflow/bemio-handles-attempt3/terminal.json).

- The original importer reproduced one KH probe leak; only that fixture-owned
  handle was closed.
- Two corrected valid imports preserved all 32 numeric/logical arrays exactly
  on each call.
- Malformed KH and Hydrostatics inputs failed as intended, with unchanged
  handle inventories. The initial inventory was empty; a nonempty unrelated-user
  handle fixture was not exercised by this run.
- All 15 copied fixture files were renamed and deleted on Windows.
- The private MATLAB helper passed the house checker. The maintained upstream
  importer retains 17 pre-existing style findings; these were unchanged from
  baseline, whose extracted copy had one additional BOM-only finding.
- No body simulation, BEM solve, or new MATLAB engine was involved. The existing
  engine remained hidden/minimized.

The related recovery/interface suite passed 209 tests with one optional skip in
13.05 s. These passes qualify resource/error handling and interface behavior,
not the pending second full canonical pipeline attempt.

## Observed failure and ownership

The first canonical CLI run terminated with exit 2 after 1989.285 s and published
no result/video. It reached surface export after completing native body and
frequency-field work. A Windows sharing error deleting `hydrodynamics/KH.dat`
then masked the original export exception.
[Run output](../../../artifacts/workflow/canonical-pipeline-attempt1/output.log).

Read-only source inspection found that maintained-fork `readCAPYTAINE.m` opened
`KH.dat` (or `KH_0.dat`) as a presence probe, reopened the same file inside the
per-body loop, and closed only the second handle. The coordinator independently
confirmed MATLAB file ID 3 referred to the exact owned temporary KH file, closed
only that handle, and verified the opened-file inventory was empty before source
editing resumed. No `fclose('all')` or unrelated-file closure is proposed.

## Source correction

- The maintained importer closes the successful KH presence-probe handle before
  entering the per-body loop.
- A small local `read_lines(fileID)` helper owns each actual sidecar read through
  `onCleanup`. Both normal return and a read exception close that handle before
  parsing continues. The original textscan format and numeric parsing remain
  unchanged. This also protects the Hydrostatics sidecar read.
- The package scratch context re-raises the same primary exception. A subsequent
  removal failure is printed and attached as an exception note rather than
  replacing it. Removal failure without a primary error is itself fatal.
- `EngineCleanupError`, or a surface-export exception carrying
  `recovery_directory`, retains the context's own system-temporary workspace.
  This supports private export replay without changing the public result format.

The coordinator completed review, tests, maintained-fork commit/pin integration,
and publication to the owned fork only. A dirty or mismatched maintained fork
still fails runtime readiness; this correction does not bypass that gate.

## Regression method

Use a private copy of an already produced Capytaine netCDF and its Hydrostatics/
KH sidecars; no Capytaine solve or body rerun is needed.

1. Record `openedFiles`, preserving any unrelated pre-existing handles. Import
   the valid fixture twice and require the same handle inventory after each call.
2. Compare all common imported numeric arrays with the retained accepted import:
   frequencies/headings, mass-related metadata, stiffness, added mass, damping,
   and complex FK/scattered/total excitation. Require unchanged values for the
   same fixture; field names/path labels are not numerical results.
3. In a fresh fixture copy, truncate KH rows to force a parse error after opening.
   Require failure plus unchanged handle inventory. Repeat with malformed
   Hydrostatics rows. This checks read-success/parse-failure cleanup; source review
   covers the guard around a textscan exception itself.
4. Confirm each fixture can be renamed/deleted after the call on Windows. Leave
   unrelated original fixtures and handles untouched.
5. Run package tests for an original stage error plus locked-file cleanup,
   cleanup-only failure, and intentional surface-replay retention. The original
   exception object/message must survive.

These checks establish resource lifecycle and array identity only. The original
surface export error remains unknown until its retained/replayed stage reports it.

The executed private harness is
`artifacts/workflow/check_bemio_handles.py` with its adjacent `.m` helper. It
copies the accepted fixture into normal system TEMP before importer calls,
reproduces/closes only the original importer’s owned KH leak, compares every
numeric/logical array recursively, and tests malformed KH/Hydrostatics inputs.
It uses `openedFiles`, since MATLAB R2025b no longer accepts `fopen('all')`.
The coordinator allocated and executed the checkpoint; this evidence update is
documentation-only. Earlier failed attempts retain their original outcomes.
