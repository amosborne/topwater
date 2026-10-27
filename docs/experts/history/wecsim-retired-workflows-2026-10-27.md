# Retired WEC-Sim integration studies

Captured 2026-10-27T05:00:02.198731+00:00. This record preserves the abandoned imported-elevation/convolution route and the directional-weight duplicate. It does not claim those unused workflows were fixed by the spectral replacement.

Native archived PRs preserve their own review history. Failed numerical studies retain their original policies and results under artifacts/workflow; this record changes no acceptance.

Original issue links inside the preserved quotations below may be retired. Use the local records: [Issue1](#issue-1), [Issue2](#issue-2), [Issue11](#issue-11), [Issue13](#issue-13).

## Associated pull requests

- PR8: Preserve imported wave history and force timing; head `805bb29b1aa5479988b61d5187f52809656cc82c`, base `25b58bda96ab0272cbfb4b0755539d5bb8ca227e`; source branch `fix/imported-excitation-history`.
- PR10: Apply all imported wave directions exactly once; head `eb722540da07c3cba737193f43712bf395ff15fa`, base `25b58bda96ab0272cbfb4b0755539d5bb8ca227e`; source branch `fix/native-imported-excitation`.

## Issue 1

Imported wave forces can lose input history or use the wrong time offsets

Forces calculated from an imported wave record can be wrong near the start and end of a simulation. WEC-Sim first cuts the record to the requested output times, discarding nearby wave samples that the force calculation still needs. The force response can also be shifted when its time offsets are not centered on zero.

## How this relates to Topwater

The working Topwater workflow uses Capytaine and WEC-Sim:

1. Capytaine calculates how waves load the boat's hull.
2. Before the motion simulation starts, Topwater uses those results to calculate the wave-driven pushes and twists at each time in the requested simulation.
3. WEC-Sim then advances the boat motion using those prepared forces and its own rigid-body and hydrodynamic calculations. Topwater does not update forces inside WEC-Sim's time-stepping loop.
4. After the motion simulation finishes, Topwater uses the saved motion to reconstruct the water surface. That does not feed back into the motion.

WEC-Sim also has its own calculation that turns a file of wave measurements into forces. During initialization, our workflow calls it with a file containing zeros, then replaces its result with the forces Topwater prepared earlier.

**We are now changing the workflow so WEC-Sim calculates the wave forces as well as the boat motion.** Topwater will prepare the wave inputs and hydrodynamic data, but will stop calculating the applied forces itself.

This makes WEC-Sim's built-in calculation part of the intended path for the real directional wave measurements, rather than just zeros during initialization. PR #8 remains open while we check whether its history/timing correction is needed to reproduce the original forces and startup behavior. Separate integration changes will remove Topwater's force preparation; this PR alone does not perform that switch.

This is unrelated to the Topwater-native solver, which remains a nonworking prototype.

## Small example

Simulate times 0 to 10 with a wave pulse at time -1. That pulse is present in the input file and should still contribute to the force at time 0. The synthetic response below makes the expected contribution easy to check.

Run from the repository root in separate fresh MATLAB processes, first on v7.1.0 and then with PR #8. For the old-version run, copy the test helper `tests/helpers/ImportedExcitationHistoryProbe.m` from that PR. It exercises the real wave/force preparation while skipping unrelated mass and radiation setup.

```matlab
addpath(genpath('source'),fullfile('tests','helpers'));
t = (-3:13)';
file = [tempname '.txt'];
writematrix([t double(t == -1)],file,'Delimiter','tab');
cleanup = onCleanup(@() delete(file));
waves = waveClass('elevationImport');
waves.elevationFile = file;
waves.direction = 0;
waves.setup([1 2],inf,0,1,10,0:10,1,1);
lags = -2:2;
c = (1:6)';
kernel = c*(2+lags)+c.^2*(1+.1*lags.^2);
body = ImportedExcitationHistoryProbe(lags,kernel);
force = body.prepare(waves,1,0);
fprintf('Channel 1 at time 0: %.1f',force(1,1));
```

| Version | Measured force at time 0 |
|---|---|
| Before the fix | 0.0: the earlier pulse was discarded |
| After the fix | 4.1: matches the independent reference calculation |

This is the first pulse case in the executed `testFullFileHistorySixChannels` regression. It is a synthetic force calculation, not a boat-motion simulation.

## Cause and fix

Use the complete supplied wave record and the response's actual time offsets. Beyond the ends of the file, continue assuming zero wave elevation, but warn that the available history is incomplete. Keep the existing gradual turn-on of waves, when enabled, and do not apply it twice.

## Acceptance

- Supplied samples before and after the output interval contribute correctly.
- Uneven and one-sided response windows use their actual time offsets.
- Short files warn and use the documented zero extension; invalid data fail clearly.
- Existing calling code, gradual wave turn-on and ordinary regular-wave behavior stay compatible.


The fix is in #8. This is a candidate for an upstream pull request because it corrects general imported-wave force preparation without Topwater-specific behavior. Upstream maintainers should review the assumptions beyond the ends of the input file and the preserved gradual-turn-on behavior.

Comment by amosborne, 2026-10-26T03:04:07Z:

The replacement Topwater handoff is now committed at [23b091c](https://github.com/amosborne/topwater/commit/23b091c): WEC-Sim receives physical directional elevation recordings and computes the applied forces. The former Python force calculation is retained only as an offline comparison reference.

Current component checks pass (36 Python backend/comparison checks and two clean MATLAB startup controls); fresh genuine coefficients were generated on the unchanged canonical inputs. The first actual force comparison stopped before MATLAB because of an LF/CRLF integrity-guard mismatch. A narrow correction has passed independent source review; actual force and motion equivalence are still unverified.

PR #8 remains the needed full-history/timing dependency. Its code/head and eleven previously measured regression passes are unchanged. Coupled directional summation, once-only application, compatible library packaging and the full comparison are tracked in the updated #2. The held dependency pin and public example outputs remain unchanged.


Comment by amosborne, 2026-10-26T12:28:53Z:

Superseded for the Topwater integration by the user's decision to pursue WEC-Sim's native spectral excitation instead of imported-elevation convolution.

The replacement work is tracked in [#12](https://github.com/amosborne/WEC-Sim/issues/12) (preserve supplied directional phases) and [#13](https://github.com/amosborne/WEC-Sim/issues/13) (consistent directional force weights), followed by exact component/coefficient handoff and complete force/startup/motion verification.

The original defect reports, commits, tests and numerical evidence remain preserved. Closing this item does not disprove the generic history/aggregation/performance findings. Its branch is retained; it is no longer a proposed change for the active Topwater route. No production pin, reference-branch merge or upstream submission is included in this transition.

## Issue 2

Sum directional imported-wave forces and apply the total once

Scope: directional imported-wave excitation for `topwater solve --solver=wecsim`, in our `amosborne/WEC-Sim` fork.

## Problem

When two wave recordings arrive from different directions, WEC-Sim's release preprocessing uses only the first recording to calculate the wave-driven force. A second problem appears when the recordings are summed: the rigid-body model can create several branches that each apply that same total, multiplying the force.

The intended behavior is simple: calculate each recording's contribution, sum those contributions, and apply the total to the body exactly once. Summation and application therefore belong in one fix. The separate history/timing correction is tracked in #1 and PR #8.

## Why Topwater needs this

Previously, Topwater calculated wave-driven force histories from Capytaine coefficients and handed those forces to WEC-Sim for the motion calculation. We are replacing that arrangement: Topwater supplies coefficients and directional water-elevation recordings; WEC-Sim calculates both the forces and the motion. Topwater reconstructs the water surface afterward.

The replacement Topwater source is committed at [23b091c](https://github.com/amosborne/topwater/commit/23b091c). It supplies no precomputed applied wave-force history. Phases, directions, time zero and the initially calm/resting condition must be preserved, with no artificial wave ramp.

## Evidence and current status

The earlier minimal two-recording reproducer omitted the second contribution, with a maximum six-channel discrepancy of 32 in that synthetic fixture. Existing selected-path corrected evidence is preserved; it does not qualify the new clean-reference candidate.

The new source candidate includes directional summation, once-only rigid dispatch and safe repeated initialization. Before deleting a generated clone, it checks internal wiring as well as block parameters; a rewired clone must be rejected and retained. The four focused group/dispatch tests and eleven history controls now pass at eb722540 on Windows/R2025b. The same-input example improves from maximum discrepancy 18 to 0; the native applied-total test also passes. Draft PR #10 contains this candidate and explicitly depends on PR #8.

Topwater component checks now pass: 36 Python backend/comparison checks, two clean MATLAB startup checks, 19 BEMIO/H5 exchange checks, and 45 producer/surface checks. A fresh run completed all 1,800 genuine Capytaine coefficient solves on the unchanged canonical inputs. These component passes do not prove force or motion equivalence.

The first actual force-comparison attempt stopped before MATLAB launched because the integrity guard compared LF source bytes with a CRLF resource copy. The narrow guard correction passed independent review and three byte-only controls. A subsequent actual comparison and independent stored-kernel replay are reported in the latest measured-progress comment below; the current discrepancy points to Topwater kernel-generation settings, not a newly demonstrated WEC-Sim convolution bug.

Private supported-API library preparation also exposed a packaging problem: a normal save migrates the release's R2020b library to R2025b and changes unrelated archive content. That whole-library save failed the strict audit and remains unaccepted. A separately reviewed narrow extraction now passes the unchanged zero/one-header audits and fresh API reopening, with no unrelated model changes.

## Acceptance

- [x] Calculate the correct six-channel total from two distinct recordings/headings, including one-recording and ordering controls.
- [x] Verify the actual applied-force bus receives that total exactly once.
- [x] Check matching/mismatching time grids, repeated initialization, regular/import transitions and rejection of modified clones without deletion.
- [x] Preserve the release library format with a narrowly audited change.
- [ ] Compare actual WEC-Sim forces against the unchanged offline reference on identical inputs, then compare all six motion histories and total-surface outputs.

The proposal covers the selected rigid imported-wave group. Mixed-group behavior must be explicit; Drag/Flex and stationary directional-spectrum weighting are outside this fix. No Topwater-prepared force override may stand in for the native applied-bus test.

Target: `amosborne/WEC-Sim:topwater/reference`; runtime pin `9952634` and canonical outputs remain unchanged. Draft PR #10 is tested at its current head; dependency integration and actual-case qualification remain open. The resulting combined library must be reconciled and retested after PR #8/#9 integration.

**Upstream candidate: conditional.** Correct summation/application is generic, but the supported group semantics, compatible library patch and independent native tests must be established first.


Comment by amosborne, 2026-10-25T23:57:21Z:

Read-only design review for this fork work item is complete; implementation remains queued after the inertia checkpoint.

Recommended representation: one aggregate `hydroForce.hf#.userDefinedFe` history plus one active imported excitation branch, preserving the original wave array. Evaluate each component's own signal/direction with PR #8 history/ramp semantics and sum fresh arrays on the common simulation output grid. Different raw-file grids are valid; do not require identical raw timestamps.

Compatibility correction: official v7.1.0/a23d299 `initializeWecSim.m` lines322-327 ALREADY rejects mixed wave types. Moving that check earlier preserves the release contract. The same release excludes imported elevation with nonlinear hydro, Morison and variable hydrodynamics; this proposal adds no support for those combinations.

Acceptance gap: historical dispatch receipt has `body_simulated:false`; it proves graph structure, not the applied force. New tests must SIMULATE the actual copied native excitation parent and capture its downstream sum/applied-force signal. Two distinct analytic six-channel components must yield F1+F2 once, excluding F1,2F1 and2(F1+F2), without Topwater force-table rebinding. Exercise ordering, repeat initialization, one/two/one components, regular/import transitions, incompatible output grids, protected custom wiring and an unchanged regular-wave control. Changed Flex dispatch needs its own >6DOF excitation-parent witness; no full Flex dynamics claim follows.

Only recognized native clones and their exact positive sum ports may be reconciled; preserve unrelated/custom blocks and reject ambiguous modified clones. No code/model/numerical job, pin change, merge or upstream write occurred in this review.


Comment by amosborne, 2026-10-26T03:46:49Z:

The actual force comparison now has a measured result. It used held fork `9952634`, including the full imported-wave history correction, and unchanged canonical inputs/current genuine coefficients.

With the exported excitation kernel limited to +/-20 s, dt=0.02 s and frequency step=0.005 rad/s, force RMS differences relative to the spectral reference were 1.052%, 2.703%, 0.252%, 2.446%, 1.301% and 4.353% for surge through yaw. Initial heave was -37.257 N versus -5.962 N; the unchanged 24.525 N force-norm calm limit rejected startup. No motion was run.

A separate stored-kernel replay at ten selected times, including startup and each discrepancy peak, matches WEC-Sim to at most 1.10e-10 N. Endpoint weighting changes initial heave by only 1.079 N, versus the 31.295 N gap, and does not explain it.

This localizes the current discrepancy to the exported finite-kernel representation/numerical settings versus the direct spectral operator. It does not establish a new WEC-Sim convolution defect. The next diagnostic varies Topwater's kernel window and frequency quadrature separately, without changing inputs or relaxing startup limits. The existing history correction and separate aggregation/once-only dispatch proposal remain distinct.

Evidence: `artifacts/workflow/native-force-comparison-attempt2/` and `native-force-discriminator-attempt1/`; exact process trees are terminal, cleanup/frozen identities verified, and leases released. Numerical equivalence, full motion and the replacement fork remain unqualified.


Comment by amosborne, 2026-10-26T06:37:56Z:

The clean-reference candidate at `eb722540` now passes all four directional-group/native-bus tests and all eleven history controls on Windows/R2025b.

The minimal example uses the same two recordings and independent six-channel expected table before and after. At t=1, the baseline prepared force is `[2, 4, 6, 8, 10, 12]`; the correct total is `[20, 19, 18, 17, 16, 15]`. The maximum discrepancy across all samples/channels falls from **18 to 0** after the fix.

The native imported source and Add-block fixture also passes its applied-total check (`1e-10` tolerance), with one imported branch, two regular-wave graph branches and one branch after returning to imported waves. The rewired-clone rejection and mismatched-grid controls pass. This is synthetic component evidence, not full boat/startup qualification.

The replacement library passes the unchanged strict archive audit with exactly one mask-header change and no unexplained/added/removed entries. Fresh MATLAB loading preserves its R2020b metadata, lock, full script and file bytes. Earlier failed runs exposed fixture workspace/count-query problems; those are preserved and were fixed without weakening force or count assertions.

PR #8 remains an explicit dependency. Both this topic and #9 touch the library archive, so integration must reconcile that binary and rerun the exact resulting head. No merge or production-pin change occurred. Actual-case force/motion equivalence remains open separately.


Comment by amosborne, 2026-10-26T12:28:56Z:

Superseded for the Topwater integration by the user's decision to pursue WEC-Sim's native spectral excitation instead of imported-elevation convolution.

The replacement work is tracked in [#12](https://github.com/amosborne/WEC-Sim/issues/12) (preserve supplied directional phases) and [#13](https://github.com/amosborne/WEC-Sim/issues/13) (consistent directional force weights), followed by exact component/coefficient handoff and complete force/startup/motion verification.

The original defect reports, commits, tests and numerical evidence remain preserved. Closing this item does not disprove the generic history/aggregation/performance findings. Its branch is retained; it is no longer a proposed change for the active Topwater route. No production pin, reference-branch merge or upstream submission is included in this transition.

## Issue 11

Reduce memory and repeated work when calculating imported-wave forces

Scope: our WEC-Sim fork and the imported-elevation path used by `topwater solve --solver=wecsim`. This is separate from the kernel-settings qualification in amosborne/topwater#35.

## Problem

To calculate wave forces, WEC-Sim looks up water height for every output time and every time offset in the force kernel. On uniform clocks, many lookups represent the same source time, but `bodyClass.userDefinedExcitation` builds the entire output-time-by-offset table.

For the unchanged canonical case, that means 3,001 output times and 32,001 kernel offsets: **96,035,001 values per table**. Each double table occupies about 0.716 GiB. Query times, interpolated elevations and one force product can nominally occupy 2.15 GiB before interpolation and MATLAB overhead.

## Measured reproducer

A fresh Windows/R2025b process called the actual held WEC method once, using one unchanged physical recording, all 3,001 output times, the admitted 320 s kernel and complete past/future elevation history.

| Measurement | Before |
|---|---:|
| Native force call | 4.742 s |
| Whole-process peak working set | 3.73 GiB |
| Physical elevation disagreement on old common support | 4.86e-17 m |

The six-channel native result is retained as a numerical replay reference. A simple extrapolation to 736 recordings is about 58 minutes per force pass; the full loop has not been timed. The measured runtime uses held revision `9952634`; the fix must preserve the newer history/group semantics from #8/#10, including their noninteger lag origin and warning/ramp behavior.

Local evidence: `artifacts/workflow/native-force-component-preflight-attempt1/probe/`. A portable synthetic before/after example must accompany the fix PR.

## Proposed change

Use standard zero-padded FFT linear convolution inside WEC-Sim for compatible uniform output grids. Keep the same imported-height interpolation, source-clock ramp, direction-specific kernel, lag origin and discrete integration weights. Use bounded row-by-row blocks for unsupported or ambiguous grids and boundary cases.

The force equation and supplied wave/IRF data remain the reference. This proposal does not calculate applied wave forces in Topwater or change kernel settings.

## Acceptance

- Preserve all eleven history and four group/native-bus controls, including legacy calls, source-clock ramps, incomplete records and modified-clone rejection.
- Add independent scalar witnesses for asymmetric/noninteger lags, past/future impulses, offset clocks, fallback selection and endpoint coverage.
- Replay the retained actual 3,001-by-six force history with maximum and time-weighted RMS differences no greater than 1e-10 of their reference scales; retain explicit zero handling without a denominator floor.
- Measure the same-input runtime and whole-process peak memory before and after. Demonstrate that the intended fast path ran.
- Provide a minimal portable example and state remaining coverage limits. Full-case force/startup/application and motion compatibility remain separate gates.

Target: our `topwater/reference`, with #8/#10 declared prerequisites. No production pin change is part of this issue.

**Upstream candidate: yes**, if numerical/API equivalence and fallback behavior are demonstrated. Avoiding these repeated tables benefits imported-wave workloads generally.


Comment by amosborne, 2026-10-26T12:28:58Z:

Superseded for the Topwater integration by the user's decision to pursue WEC-Sim's native spectral excitation instead of imported-elevation convolution.

The replacement work is tracked in [#12](https://github.com/amosborne/WEC-Sim/issues/12) (preserve supplied directional phases) and [#13](https://github.com/amosborne/WEC-Sim/issues/13) (consistent directional force weights), followed by exact component/coefficient handoff and complete force/startup/motion verification.

The original defect reports, commits, tests and numerical evidence remain preserved. Closing this item does not disprove the generic history/aggregation/performance findings. Its branch is retained; it is no longer a proposed change for the active Topwater route. No production pin, reference-branch merge or upstream submission is included in this transition.

## Issue 13

Keep directional wave-force amplitudes consistent with wave elevation

`spectrumImportFullDir`, in our WEC-Sim fork based on v7.1.0.

The specified wave amplitudes should determine both the water elevation and the wave-driven force. Changing the spacing of directional bins must not introduce an extra force scaling relative to those amplitudes.

Source review indicates inconsistent weighting: `waveClass.waveElevFullDir` includes each directional bin's width when converting the spectrum to wave amplitudes, while `irregExcFullDirF` omits that width. The active rigid-body block supplies the amplitude array directly, without compensating for it. This appears to make the force inconsistent with the incident wave represented by the same input.

### Minimal example and expected behavior

Use two directions and three nonuniform frequencies with prescribed integrated amplitudes. Compare native excitation with an independent harmonic force sum using the **actual prepared phases**, so the separate phase-import problem does not obscure this test. Repeat with directions 0/90 degrees and 0/60 degrees, encoding the same prescribed amplitudes using the corresponding native bin widths.

The predicted discrepancy depends on bin width. Before/after numerical results have **not been measured yet**; they will be added before the fix is ready for review. A separate nonuniform-direction control should exercise the existing width rule.

### Proposed change and verification

Include directional widths exactly once in native FullDir excitation, using the same rule as elevation. Derive widths from the direction vector already supplied to the force function, preserving its signature and the Simulink library. Do not change frequency weights, phases, hydrodynamic coefficients, force signs or mean-drift conventions.

Keep this fix independent of phase import. Test the scalar harmonic calculation and actual native force block before complete-packet or motion comparisons. Low-frequency support and coefficient-interpolation changes belong to separate integration work.

**Upstream PR candidate:** yes, if the minimal example confirms the source finding. Consistent wave and force normalization is a generic requirement for full directional spectra. Any upstream submission requires separate user authorization.


Comment by amosborne, 2026-10-26T12:49:27Z:

Folded into [#12](https://github.com/amosborne/WEC-Sim/issues/12) at the user's request. Phase preservation and consistent directional force weighting will be delivered as one coherent change on one branch and one PR.

This issue and its evidence are retained until the final cleanup; the combined issue carries the measured examples and remaining verification.


## Retained evidence

Historical investigations and exact producer/input identities remain in the linked expert notes and original artifacts. In particular:

- `docs/experts/wecsim_expert/spectral-force-route-review-2026-10-26.md` explains why the spectral route superseded imported elevation.
- `docs/experts/wecsim_expert/spectral-handoff-implementation-2026-10-26.md` records the earlier expanded route; its optional extensions were subsequently removed.
- `artifacts/workflow/native-convolution-source-v1/` and the original convolution attempts preserve the failed alternative without changing their acceptance policies.
- `artifacts/workflow/wecsim-cleanup-preservation1/github-before.json` retains the full API snapshot of these issues/comments, PRs/reviews and branch identities. SHA256 `9da76d3c388c7fb274e1fedbc36b5788f9b11cabaf15f299b53166513febe6ca`.

## Local source recovery

One Git bundle preserves all rebuild-store refs, including detached worktree heads and unpublished trials: `artifacts/workflow/wecsim-cleanup-preservation1/wecsim-trials.bundle`, SHA256 `d27eaf012274f080babccd1ffd1552013bae4aa9e387fff83affe98d4603c920`. The source store is shallow at releasea23d299: preserve `shallow-boundaries.txt` beside the bundle. A plain mirror clone was correctly rejected because that boundary metadata was absent; an independent bare restore using the original boundary recovered all24refs, all four unpublished trial commits and6b, and passed connectivity checks without network or alternates. Six exact dirty/untracked source snapshots and binary patches are retained beside it, with hashes in `recovery-manifest.json`. The active953 Topwater worktree, clean995 legacy checkout and official release remain separately recoverable. At capture time, no local worktree had been removed.

Independent cleanup review also found24ignored Simulink compile files in the two retired `tests/slprj` trees (1,490,350bytes). Exact snapshots and per-file hashes are now preserved in `ignored-snapshots.json` beside the bundle. At capture time, no local worktree had been removed.

Final pre-deletion API snapshot: `artifacts/workflow/wecsim-cleanup-execution1/before.json`, SHA256 `b99b4ad30653a77ab6217f6407ffa4d3e53667ca1611de608aba3f2eeed4cb7e`. Substantive issue text/comments match the original capture; only an encoded BOM prefix was normalized for display. Both raw captures remain preserved.

## Completed retirement

After validated adoption at7e1b854, GitHub natively archived PR8/10; obsolete issues1/2/11/13 were deleted with their text preserved above. Nine retired remote heads (including the old runtime branch) were removed with atomic expected-head leases. Only main and topwater/reference remain. The exact995 runtime is recoverable through archive/topwater-legacy-9952634. Fork Actions are disabled to prevent inherited docs deployment recreating gh-pages.

Eleven audited local trial worktrees and13 trial/legacy local branch refs were retired. Active953 and6b validation provenance, the shared Git store and official release checkout remain. Snapshot/hash checks preceded removal; numerical study directories were left intact. Full operation receipts are in artifacts/workflow/wecsim-cleanup-execution1.
