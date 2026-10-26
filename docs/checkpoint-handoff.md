# Checkpoint handoff

Updated 2026-10-26 UTC. Latest user steering controls scope. Read AGENTS.md,
coordination and shared-contract; historical assignments are not authorization.

## Active goal and live execution

Finish the visualization fix and simple-boat validation of the corrected native
WEC-Sim spectral path. Working branch remains `baseline/fluid-body`.

**No numerical job is live.** `native-spectral-full-force-attempt1` is terminal
success, coordinator `/root`, tool session57465 terminal. Exact Jobempty/driver0/
watchdog0/frozen/hostCIM/idle checks passed and lease released. Check the execution
ledger and actual process before another claim; do not restart a terminal job.
Manifest: `artifacts/workflow/native-spectral-full-force-attempt1/manifest.json`.
Output: `artifacts/workflow/native-spectral-full-allocation1/comparison/`.
Frozen input bundle db680661/config8427140c binds7721 files. Allocation is600s
native child,720s whole,8GiB continuous Job commit cap plus sampled working-set
limit (sampling gaps are explicit). No automatic budget increase or retry.
All producers/inputs/dependencies remain frozen until exact terminal proof.

This run compares all3001 original0..60s/.02 six-channel native forces against
the offline original spectral oracle, plus startup gates. Both per-channel
maximum error <=1e-10 reference peak and trapezoidal RMS error <=1e-10 reference
RMS; exact-zero channels require exact zero. No phase/time/gain fitting.
Measured fullforce/startup agreement is now independently accepted: all six
maximum/RMS gates passed, worst relative errors3.414e-15/4.360e-15; phase error
5.412e-15m. Controlled whole497.594s, peak commit3584798720B. Acceptance receipt
`native-spectral-full-force-review-v1/acceptance.json` SHAf460dc37 binds report
137ecc5a plus runtime/native outputs. Actual signal bus and motion remain unproved.
The old full-batch redundant isvalid cleanup warning is explicitly triaged as
secondary after successful verified cleanup; no rerun solely for that warning.

## Completed and verified

- Visualization checkpoint6ad9241 is committed/pushed. Fixed gray dashed mean
  waterline, no label; full gray hull mesh, black opening outline, solid blue
  moving contact. README explains fixed-mean hydrodynamics/exterior continuation.
  Separation is a qualitative geometry warning, not an accuracy bound.
- All50 renderer tests passed. The final public video has751 decoded H264 frames,
  800x800,12.5fps,60.08s; command136.093s. Video SHA9e654250c924749188aa184847349d5e98204369519dc41e464d4ef49b1191cf.
  Physics NPZ/PNG/protected logs unchanged. User forbids agent image review and
  previews during full rendering; significant-change previews are for USER only.
- Topwater spectral integration b1df338 is committed/pushed:38 Python,14 native
  MATLAB,9 official H5 handoff controls passed at WEC candidatee24c1a7. Raw hydro,
  radiation CIC/Ainf and original packet remain unchanged; excitation uses native
  harmonic evaluation with exact phases/weights/component excitation nodes.
- Resource probe `native-spectral-resource-attempt5` passed, exact processes/hash/
  idle checks complete and lease released before the current full run. It used
  full671-frequency/485-direction bank and full3001 wave clock;17 calls2.703s,
  whole108.672s, preparation92.132s. Firstcall.375s, steadymean.113905s/max.137577s.
  Cache sentinel only first call (direction sum4106.91). Commit peak3489710080B;
  sampledWS3948793856B with19 disclosed gaps. Resource acceptance only:
  `native-spectral-resource-review-v1/acceptance.json`, SHA3bf0c645.
- Full allocation rationale: force projection342s mean/413s max-steady plus~106s
  observed overhead;600/720s limits are bounded scheduling choices, not guarantees.
  `native-spectral-full-authorization1/qualification.json` SHA1bacfed9 binds probe,
  independent review and explicit root budget. Full wrapper cf0c8877 is reviewed.

## Clock correction and preserved failures

The separate clock-only diagnostic passed:3001 increasing samples0..60,340 entries
with at most1 local ULP difference (7.105427357601002e-15s); first difference30.06s.
Actual simulation clock exactly matches MATLAB colon; prepared clock exactly
matches integer multiplication. First17 queries are byte-identical.
NMR accepted exact own-constructor/count/endpoints/finite/monotonic checks plus
indexed1ULP bound for this private qualification. Force query times/values/gates
are unchanged. Root `run_workflow_case.m` now saves actual prepared query times
and retains native wave times separately (source SHAc8b250e6, committed after
the accepted full comparison). Clock/callback manifest72ab1c21 is reviewed.
Old private v1 files, original full comparator6e6f0a6e and all four failed probes
are preserved. Unused allocation5 was preparation only, never numerically run.
Monitor v5 handles same-handle rundown under one fixed1s grace; actual signal
is required, expiry precedes acceptance, primary failures terminate owned Job.
Receipt publication handles bounded Windows sharing locks; six cheap controls pass.
Use FILE_SHARE_DELETE when reading a live atomic receipt.

## Remaining acceptance and ownership

1. Fullforce/startup run is accepted and terminal. Next allocate actual bus only.
2. Actual native signal-bus v2 source manifest750df5b3 is independently approved.
   Requires full-force PASS plus independent receipt before allocation. All3001
   outputs must satisfy unchanged1e-10peak/exactzero gates; actual Clock values
   and output timestamps retained, constructor/1ULP proof, no resampling.
   Stdlib preparation already copied accepted input bytes to
   `native-spectral-actual-bus-allocation1/case`, config at its parent. Prep
   freeze196cee3a binds7744 originals/harness/review inputs. Recovery expert is
   authoring `native-model-watchdog-v1` to freeze new+original case and run hidden
   bus/motion batches using the approved observer; source review before launch.
3. Motion v2 source under `native-spectral-motion-validation-v1`, manifestb05351f9,
   is independently approved. Requires accepted forces AND bus before execution.
   Root later creates isolated e24 + ONLY qualified PR9 Body_Elements.slx overlay
   SHA7b0c063c8c4db2dced7e42257002989d69cdbb7a789a21f0a044ca26c793d828.
   Verify outer-model once-only force application, initial rest/zero radiation
   history and all six motions: max world-position vector error1mm and relative
   rotation angle.01deg, no fitting. Reference trajectory/forces never drive model.
4. Canonical physics files, dependency pin/reference branch and merges remain held.
   Do not mark goal complete from resource/small-test evidence.

Root is the sole numerical owner. WEC expert owns producer/MATLAB source; recovery
expert owns private watchers; NMR independently reviews, authors run no jobs.
Only private source/evidence work may occur in parallel. Recovery owns model
watchdog; WEC expert is coordinating staged motion execution and its frozen
private overlay. Preserve all existing worktrees and changes.

## Fork/reference holds and scope

Own-fork issue12/draft PR14 now point to pushed exact candidate
`e24c1a73372baf37aebe37a8b6bf0cd0d2368b0c`; source clean. Reference stays
`25b58bda96ab0272cbfb4b0755539d5bb8ca227e`; production pin9952634 remains HELD.
Old own-fork issues1/2/11 and PR8/10 are closed superseded, issue13 folded into12;
original branches/evidence retained. PR9 at67e7c93 inertia correction independently
qualified, not merged. No upstream issue/PR or merge without explicit permission.
Do not implement convolution fallback, FFT route, streaming or native prototype.
Topwater-native remains a nonworking prototype. Issues31-34 remain separate physical
qualification gaps; this compatibility check does not certify general accuracy.

Accepted operator: native-fullband-producer-attempt1/coefficients/operator-manifest.json
SHAf1f788f9, nativeMAT0effe8d1/netCDFb4a1dc07. Reuse it; no new BEM campaign.
Long historical checkpoints and failed-study detail are retained in
[the pre-fullforce archive](experts/history/checkpoint-handoff-2026-10-26-pre-fullforce.md).
