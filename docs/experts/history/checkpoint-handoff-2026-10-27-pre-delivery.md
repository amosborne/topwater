# Historical checkpoint: superseded by current handoff

# Checkpoint handoff

Updated 2026-10-27 UTC. Latest user steering controls scope. Read AGENTS.md,
coordination and shared-contract; historical assignments are not authorization.

## Active goal and live execution

Finish the visualization fix and simple-boat validation of the corrected native
WEC-Sim spectral path. Working branch remains `baseline/fluid-body`.

**CURRENT LIVE JOB:** `native-spectral-private-assembly-attempt1`, coordinator
/root, tool session **65774**, lease held. Guard config/source/output under
`native-spectral-private-assembly1/`; external freeze54a56b3e binds9227 files.
Candidate `native-spectral-private-candidate1`, evidence `native-spectral-private-evidence1`.
Independent source approval: delivery manifest7b9b9aa2/driver97464850/providerfe12e956
(Cap separately approved physical equivalence)/freezer36d230d5. Generic Python
guard4118fd03 reuses exact tested v4observer; all4 real containment controls pass.
Root one-stage limits3000s whole/12GiB commit and sampled working set, no automatic
extension/retry. Reuse accepted saved motion; unchanged production total field,
no new body integration or coefficient campaign. No producer/input/dependency
edits while live. Use SAME session65774 and live receipt with FILE_SHARE_DELETE.
Last exact live receipt464.7s running, driver33768, no error, peak commit949682176B.
Radiation960/960 (96.3s), convolution and diffraction671/671 (77.8s) complete;
new total field synthesized. Geometry/diagnostics56/751 frames done at41.5s.
Do not restart. After terminal verify exact ownedPIDs/freeze/idle/release LAST.
Then separate ordinary `uv run topwater render` on the PRIVATE candidate folder,
with source/runtime/6completedfiles frozen and bounded guard; no preview/imageview.
Then use delivery --verify-final with its newly built finalfreeze (all8candidatefiles
and assemblyreport), full decode751frames, independent final audit/logtriage.
No new body solve, field policy, guard rewrite or coefficient campaign is needed.
WEC initial logtriage lives in native-spectral-private-log-triage-v1/initial-notes.md;
request finaltriage only after complete logs. Existing restoring183%moment/132%force
and338unavailable warnings match old counts; mesh/irregular limitations persist.
Higher radiationwarningcounts5418/4506 reflect empty vs prefilled source cache,
not new warningfamily or accuracy waiver. All author agents stopped.

Coupled motion `native-spectral-motion-attempt1`, session95832 is TERMINAL SUCCESS.
All7903 frozen files unchanged; all observed owned PIDs absent via hostCIM;
independent idle readback passed and lease released LAST before current job.
Whole1382.281s, peak commit6590115840B; sampledWS6428798976B with13 disclosedgaps.
Independent acceptanced934f305 binds reporta9cadd6e: outer-force/rest/CIC pass;
max position5.975797863270701e-05m (<.001), SO3.0015698538181314137deg (<.01).
Root checkpointd9b2f1e is pushed; issue30 comment6431310948 records acceptance.

Actual bus attempt2/session39311 is TERMINAL SUCCESS, processes/frozen inputs
verified and lease released before motion staging. All3001x6 actual native bus
samples satisfy unchanged gates; actual Clock/query/output times are identical.
Reportdcb734f4, independent BUS-only acceptance547becf2; whole386.656s.

Prior bus attempt1/session31435 is TERMINAL FAILURE. Caseallocation1/watch1
preserved; exact Jobempty/frozen/hostCIM/idle checks complete and lease released.
The bus attempt stopped at209.11s on optional image metadata for NEW owned
PID27728 (membership/creation verified, wait258/exit0, no previously known image).
No force-array comparison exists from attempt1. V4 now permits fixed same-handle
1s grace for new ownership/creation-verified rundown candidates, actual signal/
membership/expiry proof and mandatory native-engine executable binding retained.
Its NNcontrols passed digest6/observer9/real-runtime5, independently accepted.
Metadata watch2 freeze failed beforeMATLAB on cross-API EXE permission semantics;
v4 compares file type across APIs while retaining FULL same-API checks.

Preceding `native-spectral-full-force-attempt1` is terminal success (session57465),
exact Jobempty/driver0/watchdog0/frozen/hostCIM/idle checks complete, lease released.
Its output is `native-spectral-full-allocation1/comparison`, freezedb680661.

The completed force study compares all3001 original0..60s/.02 native forces against
the offline original spectral oracle, plus startup gates. Both per-channel
maximum error <=1e-10 reference peak and trapezoidal RMS error <=1e-10 reference
RMS; exact-zero channels require exact zero. No phase/time/gain fitting.
Measured fullforce/startup agreement is now independently accepted: all six
maximum/RMS gates passed, worst relative errors3.414e-15/4.360e-15; phase error
5.412e-15m. Controlled whole497.594s, peak commit3584798720B. Acceptance receipt
`native-spectral-full-force-review-v1/acceptance.json` SHAf460dc37 binds report
137ecc5a plus runtime/native outputs. Actual signal bus and coupled motion are independently accepted (see current evidence above).
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

1. Forces, native bus and coupled motion are independently accepted and terminal.
   Current private assembly65774 must finish exact field/export/identity/terminal gates.
2. After terminal process/hash verification, independently review actual outer
   once-only force, initial rest/zero CIC history/full replay and all6 world poses.
3. Motion v2 certificate manifestb05351f9 and staged validator v3 manifest2f37980c
   are independently source-approved. Requires accepted forces AND bus before
   preparation/execution. Prepare-only imports science under a separate root
   lease but starts no MATLAB; run-staged verifies the final frozen case and
   rewrites no inputs. Raw CRLF bytes remain frozen; only canonical LF conversion
   proves the exact single logging-hook recipe.
   Root created detached `wecsim-spectral-inertia-validation` at e24 with ONLY
   qualified source/lib/WEC-Sim/WECSim_Lib_Body_Elements.slx overlay
   SHA7b0c063c8c4db2dced7e42257002989d69cdbb7a789a21f0a044ca26c793d828.
   Verify outer-model once-only force application, initial rest/zero radiation
   history and all six motions: max world-position vector error1mm and relative
   rotation angle.01deg, no fitting. Reference trajectory/forces never drive model.
4. After accepted motion, deliver a PRIVATE candidate through the SAME production
   pose conversion, restoring diagnostics, total-surface reconstruction, Result
   export/result.png and public render CLI. Reuse saved new motion; do not rerun
   body integration. Exact8files; provenance/delta reports outside the folder.
   The prior public video proves renderer behavior only, not corrected-path export.
   Necessary same-production field RHS solves are authorized; no new coefficient,
   mesh-refinement or convergence campaign. Old completed-field arrays cannot be
   reused with new motion. WEC owns single surface_context.py b2881929 and driver;
   Cap independently reviewed context. Default centroid BEM quadrature is unchanged.
5. Canonical physics files, dependency pin/reference branch and merges remain held.
   Do not mark goal complete from resource/small-test evidence.

Root is the sole numerical owner. WEC expert owns producer/MATLAB source; recovery
expert owns private watchers; NMR independently reviews, authors run no jobs.
Only private source/evidence work may occur in parallel. Reviewed model-watchdog
v4 supports BUS/staged-Motion and within-verification hash deduplication.
The tiny v2 fixture exposed Windows path-stat versus handle-fstat ctime semantics:
all fields agree except stable ctime values, confirmed in
`digest-boundary-check-attempt1/stat-observation.json`. V3 retains separate full
path-before/after and handle-before/after checks, common-field cross identity
(file TYPE, not cross-API permissions) and fresh pre/post scopes. No physical
guards/tolerances change. Do not edit live V4 inputs. WEC expert owns
approved staged motion v3 source. Preserve all existing worktrees and changes.

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
