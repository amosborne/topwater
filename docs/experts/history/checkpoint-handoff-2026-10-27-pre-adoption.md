# Historical checkpoint; superseded by the active adoption goal

# Checkpoint handoff

Updated 2026-10-27 UTC. Latest user steering controls scope. Follow AGENTS.md,
coordination and the shared physical contract. Historical assignments are not
permission to resume work.

## Current goal and status

The visualization fix and corrected-path simple-boat validation goal is complete
at pushed checkpoint8703b02. Active implementation branch: baseline/fluid-body.

ARCHIVAL CORRECTION (verified2026-10-27): GitHub DOES support native individual
PR archival. Topwater `is:pr is:archived` returns17 PRs13-29; merged11/12 are not
archived. WEC-Sim returns0 archived PRs. Earlier assistant claims that only closed
states/archive tags existed were WRONG. Official feature:
https://docs.github.com/en/communities/moderating-comments-and-conversations/archive-pull-requests
User means this actual feature for retired WEC PR8/10, not mere closure/export.
Cleanup plan corrected accordingly; no WEC archival/deletions performed yet.

LATEST USER REQUEST: both PR9 and PR14 are now merged. Make a WEC-Sim fork
cleanup PLAN FIRST: PRs to archive, issues to delete, branches to delete. No cleanup
branch/issue deletion is authorized yet. Draft plan:
experts/wecsim_expert/fork-cleanup-plan-2026-10-27.md. Main/ref kept; six topic+
dev/gh-pages proposed for deletion; old995branch held pending validatedpinmigration.
Issues1/2/11/13 proposeddeleteafteronehistoricalexport; keep3/4/5/12. PR8/10mustbe
nativelyARCHIVED(notjustclosed); keepmerged6/7/9/14. No upstream changes or newrelease rebase.

Actual WEC reference now6b774228e1f6b31dee473e626204b0ef65e5ec38. PR9 user squash
merged4cf628e, PR14 user squash merged into6b774228. Their finaltopicheads were
b2ceac6 and1556710. PR9 reduced8files(oneproductionproperty), PR14 reduced10files
(two productionfiles). PR14 maintextgeneric/noTopwateranalysis; all10native
regressions passed, independentreceipt4ffbf7f8. EarlierWEC customAPI/BEMIOwriter
extensions removed. Issue3andPR9 ownstatuscomments archivedthenremoved asuser
requested; authoritativebodiescurrent. No play-by-play comments going forward.

NO LIVE NUMERICAL JOB. Latest standard-adapter-python-attempt2 is TERMINAL PASS,
all38focusedPython tests2.99s. Actualnewsourceimportpath asserted. Rootverifiedall
ownedPIDs absent viaCIM, frozenfilesmatch, idle/readback andlease releaseLAST.
Attempt1 had26pass/12SETUPERRORS from inaccessibleglobalpytest-tempdirectory;
noimplementationchange, attempt2 usednew workspacebasetemp andpassed.

PRESERVE ongoing adapter work (prioritypausedforcleanupplanning, notcancelled):
Topwater isolated artifacts/worktrees/standard-wecsim-input,
branchwork/standard-wecsim-input, committed9539192bb733663584483f33ef2fb3de273d6fcf.
Nonauthorsourceapproved e1998261; native21validation/handoff tests NOT YET RUN,
canonicalforce/fullinitializer/bus/motionnewroutequalification NOT YET RUN.
No further numerical job is queued under the plan-first request.
Rootbranchbaseline/fluid-body source stilloldexpandedpath; adaptercommit hasNOT
been cherry-picked, productionpin995unchanged. Do notdeleteitsoldruntimebranch.

Existing documented publicstruct input was proved by tiny13.75s syntheticgate,
receiptb787454d (CIC/Ainf/mass/restoringexact; wave1.11e-16/force8.88e-16). It is
not genuineboat/fullinitializationproof. New spectral adapter usesofficialBEMIO
readCAPYTAINE/radiationIRF/rebuildHydroStruct andSTANDARDnested excitation fields;
noH5spectral/customview/accessor/newAPI/fExtinjection/fakeunusedA/B. LegacyH5/IRF
studies remain underoriginalpolicies. Sourceauthors stopped; preserveallworktrees.

Isolated WEC worktrees: wecsim-standard-input1556710 andwecsim-inertia-cleanupb2ceac6;
old e24/67e sources preserved. Inertiacleanupbyte-equivalence proof confirmsALL163
retained source/test/library files equal oldqualified67ebytes; onlyunusedauthoring
tools+READMEremoved. Native3+sample8evidencevalidforidenticalexecutablebytes;
no newnativeexecutionatb2wasclaimed. Oldthree authoringtools stillversioned67e.

Prior visualization/simpleboat expanded-route goal completed8703b02. Its evidence
below remains historical and does not establish the newly adapted standard route.

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
