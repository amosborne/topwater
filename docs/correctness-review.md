# Correctness review of the current WEC-Sim workflow

This review follows the code used by `topwater solve --solver=wecsim`, including
supported boat inputs whose inertia differs from the example. Capytaine, BEMIO
and WEC-Sim experts traced their respective stages. It does not cover unused
Drag/Flex branches, and it is not a claim of physical
accuracy for the completed baseline.

The authorized transition gives WEC-Sim ownership of excitation as well as motion.
Multiple directional imported-wave branches therefore belong to the intended
codepath and require history, summation and actual applied-signal checks. The
saved baseline used Topwater-precomputed forces; it does not verify that transition.

The baseline ran successfully. The review found integration checks and cleanup
that need strengthening; it did not establish a new native Capytaine sign or
scaling bug.

## Findings and required examples

| Finding | Why it matters | Small before/after check | Status |
|---|---|---|---|
| Cleanup can hide the original MATLAB error | A failed run can report only the later cleanup failure, losing the useful diagnosis and transcript | Force a primary error and a cleanup error together; the corrected message must retain both | Fixed; old-code failure reproduced and corrected backend suite passes |
| One MATLAB setting is not restored | WEC-Sim changes `ErrorIfLoadNewModel`, which can affect later work in the shared session | Start with a different setting, exercise success/failure, and verify it is restored | Fixed; both old-code witnesses fail and corrected native setting tests pass |
| Imported-coefficient checks can overlook NaNs | MATLAB's `max` omits NaNs, so a comparison can pass despite an invalid coefficient | Inject one NaN into otherwise valid imported data; the corrected boundary must reject it clearly while a finite control passes | Fixed; eleven standalone guard tests pass |
| The active File Solid needs its inertia products reordered | The release supplies the values in a different order from the native block; other boats or adjusted added mass can expose it | Apply independent forces/torques to the active rigid Solid and infer its matrix from acceleration, before and after the reorder | Current pin contains the fix; measured dry and adjusted fixtures pass to roundoff. Minimal fork rebuild remains outstanding |
| Capytaine importer leaves a file open | Windows cannot remove/move an input file after import | One sphere import leaves one file open before the fix and zero afterward | Fixed in current pin; measured example and tests in fork PR #6 |
| Imported paths are interpreted as formatting | Logs can warn or print a truncated Windows filename | Print the same path before and after; corrected output is literal and warning-free | Fixed in current pin; example and tests in fork PR #7 |

The three new fixes passed independent source review and focused Windows/R2025b
checks: 21 Python backend tests and 13 MATLAB tests. Old-source error-masking and
setting failures were reproduced; finite-data witnesses evaluate the exact old
predicate before checking the production guard. Evidence:
`artifacts/workflow/current-path-checks-attempt2/`. These are integration checks,
not a new BEM/body run or evidence that the saved canonical result was invalid.

## What the source review supports

- Capytaine radiation, free-surface and complex-exponential conventions match
  the wrapper at zero forward speed. Diffraction and Froude-Krylov are added once.
- Exact facet hydrostatics and direct native infinite-frequency added mass are
  the consumed values. Rounded sidecars and BEMIO's finite-band added-mass estimate
  are deliberately overwritten; they are not the physical authority here.
- The inspected BEMIO units, phase conversion and H5 permutations are consistent
  in source for this single-body six-mode path. The richer coefficient handoff
  now passes 19 MATLAB controls, including official H5 readback and preservation
  of the original full-band radiation kernel and infinite-frequency added mass.
  These synthetic controls do not establish actual boat force equivalence.
- The saved baseline compares its applied excitation bus with the Topwater force
  history. The replacement must instead verify WEC-Sim's computed total against
  the actual applied bus, plus an independent offline force comparison. Neither
  force nor motion equivalence of the replacement is established yet.
- MATLAB's HDF5 identifiers have automatic close callbacks. Missing explicit
  `H5T.close` is not itself a demonstrated leak; no writer fix is proposed.

The native coefficient exchange and existing surface fixtures pass 45 Python
tests on the reviewed producer. Complex phase, axis ordering, MATLAB exchange,
periodic closure, FK/diffraction summation and invalid inputs are checked without
a new BEM solve. Evidence: `artifacts/workflow/native-producer-tests-attempt1/`.
The [measured inertia proof](experts/wecsim_expert/inertia-adjustment-proof-2026-10-25.md)
covers dry and adjusted native matrices plus an algebraic residual; actual
H5/CIC loading, dynamic residual forces and full boat motion remain outside it.

## Physical questions still open

[Radiation phase #31](https://github.com/amosborne/topwater/issues/31),
[mesh and irregular-frequency uncertainty #32](https://github.com/amosborne/topwater/issues/32),
[field/load differences #33](https://github.com/amosborne/topwater/issues/33), and
[linear applicability #34](https://github.com/amosborne/topwater/issues/34)
remain unresolved. Saved-sample force balance and kernel reconstruction establish
implementation consistency, not time-step, memory-length or mesh convergence.
The rotating adjusted Solid and fixed-world residual also need their stated
finite-attitude limitations preserved.

The first diagnostics should reuse stored data or small synthetic controls:
the existing phase discriminator, distinct matrix/kernel sentinels, a restoring
perturbation, and an independent orientation check. No new full boat/BEM run or
performance campaign is allocated by this document.

## Scope and evidence

Reviewed identities: Topwater implementation `6dbc7f6`, WEC-Sim fork
`9952634f6f2e2746a1626484e27c3003d1e1d5a9` over release v7.1.0 `a23d299`,
Capytaine 3.0.0, MATLAB R2025b Update 1 and matlabengine 25.2.2. New fixes are
being prepared against those sources. The execution ledger owns subsequent runs.

See the [handoff](checkpoint-handoff.md),
[Capytaine current facts](experts/capytaine_expert/current-state.md),
[warning triage](cli/warning-triage.md), and
[fork regression evidence](experts/wecsim_expert/fork-rebuild-regressions-2026-10-25.md).
Every fix needs the same small input before/after, actual outcomes, independent
review and a regression test. Fork proposals also state whether they are an
upstream candidate and why. Upstream publication and merges require separate
authorization; the Topwater pin remains on hold.
