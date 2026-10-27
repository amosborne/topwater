# Checkpoint handoff

Updated 2026-10-27 UTC. Latest user steering controls scope. Follow AGENTS.md,
coordination and the shared physical contract. Historical assignments are not
permission to resume work.

## Current goal and status

Goal: finish the visualization fix and simple-boat validation of the corrected
WEC-Sim spectral path. Active implementation branch: baseline/fluid-body.

**NO LIVE NUMERICAL JOB.** Final private eight-file audit is terminal success.
Root verified all owned process identities absent through host CIM, frozen inputs,
separate idle-ledger readback and lease release LAST for each allocation. No
numerical job or producer edit is queued. Final independent end-to-end acceptance
76f5616a and warning triage147c92ad are complete. This checkpoint records the completed scoped delivery; verify its
commit is pushed before marking the active goal complete. Issue30 final update:
https://github.com/amosborne/topwater/issues/30#issuecomment-6431729288

Pushed checkpoints:6ad9241 visualization, b1df338 spectral integration,
d348e4c timestamp/full-force qualification, d9b2f1e accepted coupled motion.
Own-fork issue12 and draft PR14 now correctly state boat force/bus/motion accepted,
with the separately qualified PR9 inertia overlay explicitly identified.

## Actual delivery and evidence

Private candidate: artifacts/workflow/native-spectral-private-candidate1/
contains EXACTLY boat.npz, sea.npz, result.npz, result.png, video.mp4, sea.log,
result.log and video.log. Reports/cache/provenance are outside that folder under
native-spectral-private-evidence1/. Public canonical physics remains protected.

- Native full-force/startup: acceptance f460dc37. All3001x6 original samples pass
  unchanged1e-10 peak/RMS gates; exact-zero channels exact. Worst relative maximum
  3.414e-15 and RMS4.360e-15. No fitting, phase/gain change or external driving.
- Actual native bus: acceptance547becf2. All3001x6 samples pass; actual query and
  output times equal the prepared integer-multiply clock exactly.
- Actual coupled motion: acceptance d934f305, reporta9cadd6e. Maximum world-position
  vector difference0.059758mm (<1mm), relative SO3 angle0.001569854deg (<.01deg).
  Outer once-only application/rest/zero initial radiation/full CIC replay pass.
  Negative velocity history is inferred from fresh/native-zero initialization;
  persistent negative velocity arrays were not directly measured.
- New total-surface/Result/PNG assembly: reporta6c096b4, source97464850/contextfe12e956,
  operatorf1f788f9/production surface240198ee, external freeze54a56b3e (9227files).
  Saved new motion was reused; no new body integration or coefficient campaign.
  Empty source bank required5760 radiation and671 combined-diffraction RHS solves
  with SAME accepted17-panel/default-centroid BEM operator. No old field reused.
- Full ordinary public render CLI:751 frames,800x800,12.5fps,60.08s encoded duration;
  physical snapshots0..60s. Video SHAe5fce0e4b89114268c127a5dea57b4e51a108117b5d889576ddb33dbf902c2c4.
  video.log ENDsuccess/warnings0. All9240 frozen render inputs unchanged.
- Final audit: native-spectral-private-evidence1/final-report.json, final freeze
  84c42f0f (9237files). Standard finite Result/schema,3001 motion/751 surface frames,
  exact accepted new poses/0..60 endpoints/input identities/log END/PNG signature
  and full751-frame decode passed. Root terminal/CIM/frozen/idle/release complete.

Generic Python guard4118fd03 reuses tested v4observer53acfac6; four real containment
controls passed before delivery. No further guard changes/retries are needed.

## Visualization and output requirements

Current renderer141046c5 passed50 display controls (mean-waterline-preview-attempt3).
Current plotf39a7572 passed23 focused controls (seven-tick-plot-tests-attempt1).
Do not attribute the earlier286-test workflow run to these later source hashes.

Water surface is visible up to the moving contact; full gray hull mesh, black
opening outline and solid blue waterline retain their layer order. Fixed dashed
mean-waterline reference has NO label. It is computed once at initial pose/flat
z0. Separation is a qualitative geometry departure, not an accuracy bound.
Result.png uses common translation/rotation row scales and seven nice symmetric
ticks, zero centered, edge ticks at plot limits. README explains fixed-mean
hydrodynamics, radiation/restoring feedback and explicit exterior-trace continuation.
User forbids previews or agent image inspection during full rendering. No images
were viewed in this delivery; automated complete video decoding is allowed.

## Field deltas, warnings and timing

New total-minus-incident projected-area/time RMS:5.98768mm over full domain,
2.24251% of incident RMS; one-beam RMS7.99367mm/2.96028%; maximum65.0892mm.
These are field contributions, not numerical error bounds. Direct old/new field
subtraction is explicitly not comparable because saved physical XY grids differ;
no regridding or fitted comparison was performed.

Known physical warnings persist: calm restoring departures183%moment/132%force,
338 unavailable flat-reference samples, represented wet-set departure29.15%,
coarse mesh/irregular-frequency and provisional causal field projection. Issues31-34
remain unresolved; compatibility does not certify physical accuracy. Higher new
radiation-warning counts reflect empty vs prefilled source-cache accounting.
Final WEC triage is under native-spectral-private-log-triage-v1 (report147c92ad).
NMR final independent acceptance is private-delivery-review-v1/acceptance.json,
SHA76f5616aac1a515ba26d8ee49b7f9277ba938b07e175bb96d96b44df94326215.
It binds all8actualfiles/finalreport624ed7bb/source/test/force/bus/motion/triage evidence.
Versioned evidence: experts/wecsim_expert/spectral-end-to-end-qualification-2026-10-27.md.

Measured motion validation1382.281s; surface/export968.172s; render CLI153.797s
(controlled155.531s); final audit29.219s. Geometry505.6s is mostly wet-set area335.4s
and delta metrics143.6s; waterline intersections11.5s. Detailed accounting and
MATLAB timing gaps: docs/workflow-timing-2026-10-27.md. No speed campaign started.

## Holds and owners

Production dependency pin9952634 remains HELD. Reference25b58bda stays unchanged.
WEC topic exact e24c1a7 is pushed; PR14 remains draft. Motion private detached
wecsim-spectral-inertia-validation overlays ONLY PR9 library7b0c063c. PR9 exact
67e7c93 remains unmerged. No upstream issue/PR, merge or pin promotion is authorized.
Topwater-native remains a nonworking prototype. No FFT/convolution fallback,
streaming, new coefficient/convergence/mesh campaign or another body run is queued.

Root owns integration/checkpoints/numerical execution. WEC owns fork/adapter,
Cap approved final provider equivalence, NMR independently reviews physical gates.
All scientific authors run no jobs. Preserve every worktree and failed study.
Historical detail: experts/history/checkpoint-handoff-2026-10-27-pre-delivery.md
and checkpoint-handoff-2026-10-26-pre-fullforce.md.
