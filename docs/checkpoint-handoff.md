# Checkpoint handoff

Updated 2026-10-27 UTC. Follow AGENTS.md and coordination. Latest user steering
controls scope; historical plans and diagnostic gates do not authorize campaigns.

## Active goal and user acceptance

Resolve Topwater #31. User CONFIRMED **1 cm maximum numerical water-elevation
error INCLUDING reference uncertainty, against the chosen LINEARIZED model**.
This is reconstruction accuracy, not a nonlinear-model or experimental comparison.
Freeze the accepted WEC-Sim motion/forcing, mean hull and chosen hydrodynamic data.
Compare the same exterior locations and clocks through the full60s. Preserve
incident/diffraction/radiation; never adjust motion or phases to obtain agreement.

The earlier23% figures are transfer-function diagnostics, NOT centimetre errors
in the saved water. The actual-trajectory seven-probe screen is now complete,
but its reference remains unqualified. Do not claim certified error or completion
from method disagreement, proxy residuals or a few points.
The user explicitly rejects pursuing perfection for its own sake.

Priority: #31, then measured end-to-end speed, then #32-34. Upstream contributions
are deferred. Active branch remains baseline/fluid-body; no main pushes,
assistant merges, bookkeeping PRs or upstream contact/issues/PRs.
No WEC-Sim rerun, body-source correction, mesh campaign or native dynamics is
currently allocated. Root owns serial numerical execution and integration.

## Latest implementation and immediate work

Full-H offline implementation source235989db and focused tests544ec0fa/19595990
passed independent source reviewa75ccbde27caeddd96e5ad97deccd72793aa57ddc3c96cbf787ec8614a03d3fb
at artifacts/workflow/issue31-offline-source-review1/acceptance.json.
Actual focused suite:42 PASSED in17.81s, test_offline_radiation.py and
test_surface_reconstruction.py. Evidence: issue31-offline-tests-attempt3/terminal.json
and issue31-offline-tests-allocation3/child.log +guard-receipt.json. Earlier tests
attempt1/2 failed sandbox temporary-directory/rename permissions, not numerical
assertions; attempt3 used workspace OS/pytest temp and host execution with the
same source/settings. Guard0/driver0/Jobempty/frozen; hostCIM exactowned absent;
finish74554 terminal0, separateidle/readback/release-last. Ledger idle/null and
lease absent. Source freeze lifted. No full canonical solve/render rerun.

Implementation preserves raw startup, full complex H and all native samples;
logs policy/startup/pre-response per spatial chunk. ONE production readiness flag
remains FALSE until whole-domain qualification. Root updated docs/cli/solve.md
to describe pending settings/failures/limits, not claim final accuracy.

NEXT CAP source-only assignment: concrete whole5929x751 qualification helper/plan
using accepted saved motion, unchanged incident/diffraction/trace/geometry checks;
include a SHORT field+FFT runtime/memory pilot before any large allocation.
Sourcebank/chunking cost must be measured. No WEC rerun or gate flip. Prospective
per-point reference controls and uncertainty, matched clocks/basis, independent
frequency-sum spot checks and common incident/diffraction identity are required.
Root executes serially only after review.

Input-tail helper/plan are source-ready at issue31-input-tail-screen1/, author
hashesb47b4529/33c39bef; NMR review pending. Use final accepted hashes for execution.
Decision screen uses full-H tail envelope PLUS existing per-probe Uref_est <=10mm;
no default0.5mm assumption, no automatic saved-a semantics charge for an unchanged
full-H reference. Passing is necessary only, not whole-field or tail qualification.
Sparse sample input values use direct DFT, not complex interpolation with an
arbitrary time-origin phase. Preparation source is under input-tail-allocation1/.

## Completed study evidence and current boundaries

Both actual-height screen and tail extension are terminal0/guard0/Jobempty/frozen.
HostCIM, independent idle readback and release-last completed; latest finish
session36942 exited0 (earlier screen finish11288 also exited0).
Ledger idle/owner null and lease absent. No numerical process is allocated.
Revalidate current liveness before the next job; a checkpoint alone is not proof.

Actual report: artifacts/workflow/issue31-height-screen1/results/report.json,
with height-evidence.npz alongside. ALLsix accepted NEW motions were used.
Cosine24 versus full-complex reference-candidate disagreement reaches42.6603mm
at44.68s on an exterior continuation trace; the direct near exterior probe
reaches31.9271mm and far probe3.46673mm. These are method differences, NOT
certified errors against the linearized model. The conservative observed
sensitivity sum reaches2.91712mm, dominated by future-extension/cut controls.
Missing response above24rad/s remains unknown; reference_qualified=false and
full_surface_acceptance=false. Fourteen controls took0.844s (guard2.094s).

Independent actual-results review accepts this screen and its limits:
artifacts/workflow/issue31-height-review1/acceptance.json, SHA9d55d233d3ae1534a35023dfc222cf6d3fee095e42b4d24319ca0abf1cf62664.
Next: the smallest justified extra frequency-tail/reference check, then a reconstruction correction and whole-domain
qualification. Do not restart body/sourceK/mesh or primitive precision campaigns.
The user explicitly separates linear-model numerical accuracy from the harder
nonlinear/experimental model-error question.

The PRIVATE24-to48rad/s same-operator extension is now complete:
artifacts/workflow/issue31-tail-extension1/results/report.json, transfer-extension.npz
and height-evidence.npz. Child26.422s; exit0/guard0/Jobempty/frozen; hostCIM all
owned identities absent; finish session36942 terminal0 with independent idle
readback/release-last. Ledger idle/owner null, lease absent. Recomputed20/22/24
overlap has EXACTzero discrepancy, and old0..24 samples are unchanged. The maximum
actual48-minus24 elevation change is0.0591576mm. The36-to48 absolute input-weighted
envelope is at most0.0386197mm; Uref_est reaches2.86777mm. These are seven-probe
engineering convergence evidence, not full-domain acceptance or an infinite-tail
bound. Independent actual-results reviewf0904700ca9c9988db34e1eda173ecec5bcf4bf086e3ae42fc1d342ae6e828e6
accepts this conditional engineering evidence; source review58ae4b8f passed.
Review receipt: artifacts/workflow/issue31-tail-review1/acceptance.json.
Root accepts48band as a CONDITIONAL comparison reference with an explicit0.5mm
remaining-tail engineering allowance, pending analytic feedthrough separation and
no unresolved high-frequency growth/resonance. This is an assumption, not a
proven infinite-tail bound. Existing sparse Heave32/64 evidence does not prove
all-mode tail absence. Final qualification must state assumptions/uncertainty
and cover the full declared region, not just these seven probes.

Next bounded CAP source-only work: directly test a full-complex signed-lag causal
candidate against saved-motion reference, considering analytically derived
instantaneous acceleration before inversion. Future-extension sensitivity is NOT
a bound on actual negative-lag contribution; do not assume truncation passes.
The finite-low discrete operator endpoint differs from native symbolic-infinity
(which switches Green split). NMR provisionally accepts the algebraic endpoint
derivation under explicit conditions, with the jump term exactly once and
invertibility/conditioning checks. No fitted feedthrough or sourceK correction.
This is private implementation/testing, not a selected production correction.
Engineering convergence and explicit tail assumptions are sufficient; do not
demand an infinite-frequency theorem. No numerical job currently allocated.

The private causal candidate was EXECUTED and REJECTED. Its source/plan remain
under artifacts/workflow/issue31-causal-screen1/; attempt1 failed cache-directory
setup, attempt2 failed symbolic DOF use; both exact sources are preserved there.
Only setup fixes were made before attempt3, source7353191b/plan059acd34,
review803c76ac. Actual attempt3 report/arrays:
artifacts/workflow/issue31-causal-attempt3/results/. Child1.078s; exit0/guard0/
driver0/Jobempty/frozen; hostCIM allowned absent; finish1577 terminal0, independent
idle readback/release-last. Ledger idle/owner null, lease absent.

Delay/fold/direct numerical controls PASS, but candidate-vs-regular-reference
max194.933mm at43.74s/probe4; old cosine max42.593mm. Candidate band36-to48
sensitivity11.125mm and dt sensitivity1.397mm. Full-H versus saved-A*a reference
change is only0.286675mm. Endpoint condition4920.88, residual4.44e-16. This rejects
the destructive causal-only projection, not a sign/scaling/setup failure.
Independent actual review: artifacts/workflow/issue31-causal-review1/acceptance.json,
SHA8698b3ca8e2c5eac18097424cd9f4a32fd36f688295a7698fb0676aa2106153d.

CURRENT AUTHORIZATION: CAP implements intact full-complex offline FFT(v) surface
reconstruction in _capytaine_surface.py and focused field tests. No fitted
endpoint, runtime asymptotic decomposition, causal cut, body/sourceK/mesh change,
motion rerun or forcing change. Preserve incident/diffraction, exterior trace
policy and physical/geometry failure gates. Reuse sourcebank and spatial chunks;
generic inputs, no canonical whitelist. Explicit zero prehistory and smooth
future continuation; preserve RAW t0 by default, record startup/pre-response and
validate initial memory against derived A*a. No new initial-value override API
or automatic zeroing of physical direct response. If a gate conflicts, report
it rather than silently relax it. Any later qualified correction is budgeted.
No tests or
numerics until root's serial allocation; independent source/results review first.
Stop authoring at source-ready, not production acceptance. Whole5929x751 total
elevation qualification, full-run controls, regressions, logs/docs and outputs
remain required. A universal causal-only kernel is not an extra user requirement.

NMR owns a bounded read-only practical tail-qualification proposal. Sparse
regular-H32/48/64/96 maxima .466/.248/.275/.835 show no-growth is NOT established.
The0.5mm tail assumption therefore remains conditional, not accepted as an
unconditional bound. Seek actual input-weighted tail evidence proportional to
10mm, not primitive/transfer perfection or an automatic frequency/mesh campaign.
Offline full-complex synthesis does not itself prove exact physical rest-start
causality; future sensitivities are selected-continuation tests, not a bound on
all possible anti-causal influence. Do not claim reference or goal complete yet.

Next cheap reference task: TOP/receipt_reliability authors source-only
artifacts/workflow/issue31-input-tail-screen1/ helper/plan, using existing sparse
H48/64/96 and actual NEW input spectrum. NMR independently defines/reviews the
weighting and explicit envelope assumptions. No new BEM or run delegated. First
measure actual input weighting, then decide whether ONE justified48-to96 response
extension is needed. Do not interpolate sparse samples as known tail or reuse
unsupported0.5mm blindly. Root allocates all execution only after source review.

The preceding bandwidth-control report is
artifacts/workflow/issue31-green-bandwidth1/results/report.json:
known-causal point flux at actual deep/shallow centroids, no body solve. The
shallow omitted-band discrepancy reaches1.31% of its flux-transfer scale; the
same-band cosine/KK discrepancy is about0.116%. This establishes a finite-band
contribution, NOT the cause of the entire23% boat-mode residual or an actual
water-height error. Independent review0e8722d0 accepts only that diagnostic.

Completed PRIVATE seven-probe actual-height screen, ALLsix NEW saved motions:
Use existing full-complex transfer cache, not another primitive-precision study.
Reproduce current cosine24 reconstruction with NEWv; compare full-complex
frequency reconstruction with explicit zero prehistory, smooth future-extension,
padding/time-integration and supported-band controls. First screen can reject
or expose reference uncertainty; it cannot qualify the full published region.
Reference uncertainty allocation2mm is internal planning within10mm,
not a request for machine precision or an extra user physical requirement.
Unknown frequency-tail error must remain explicit, not silently assumed zero.

Owners: capytaine_expert authors the bounded helper/implementation; independent
numerical_methods_researcher reviews reference design, source and actual results.
Root selects parameters/leases/acceptance. No next numerical launch until its
bounded next plan and prospective controls have independent review. The screen's FFT
periods256/512/1024s must be truly distinct; dt.02/.01, smooth future continuations
10/20s, bands12/18/24. Preserve all original motion knots and compare all3001
input times. Interior future cuts20/40s test earlier-prefix independence in
addition to post60 extension/prestart checks. Unknown >24rad/s support remains
explicit; numerical sensitivities alone are not a certified tail bound.

## Exact saved inputs and coordinate scope

- Accepted NEW motion: artifacts/workflow/reference-adoption-case3/native/motion.mat,
  SHAaf127598fd6bcf3b8bd9bbe727e10f8109df8f9a73f7d87411aa18f69f9c31df.
  Native time3001x1; q/velocity/acceleration3001x6. Mode order is fixed mean native
  Surge/Sway/Heave/Roll/Pitch/Yaw, not Euler-rate/world-velocity substitution.
- Accepted snapshot: same campaign motion-comparison.npz, SHAd588619a. Current
  completed-field cache q/v/a payloads match it exactly.
- Current field cache: artifacts/workflow/reference-adoption-private-evidence1/
  field-work/surface-recovery/completed-field.npz, SHA2000203da3645e3861a667cfd0fe0ee287bc352a3cc8b8a2f3367ba30bc63ec7.
  It contains all NEW motion, poses/coordinate transforms,751 surface times,
  xy_native_m5929x2, total/incident elevation751x5929. It does NOT retain a
  radiation/diffraction source bank or separate radiation/diffraction surfaces.
- Existing transfer cache: artifacts/workflow/surface-response-attempt1/field/
  response.npz, SHA110380140cbab1dbcd2e6f72a6e5efc7ae0c0a6598af3a5fe62959065e7e8492.
  Complex Hv960x7x6 at .025..24rad/s; old radiation response arrays must NOT be
  reused with the new trajectory. Recompute from NEWv. Original/current mean
  hulls differ only by recorded roundoff up to1.4844e-18m, not byte identity.
- The first three probes are direct [6,0],[0,1.25],[3.5,0]. Four others are
  exterior traces for two mean-interior targets; report xy AND target_xy.
  Seven probes do not establish a maximum over5929x751. Final scope must cover
  all declared valid published exterior nodes/times or supply a valid bound;
  continuation regions and their shared policy must be reported explicitly.

## Preserved diagnostic findings

Source and canonical outputs remain unchanged by #31; all new kernels are
PRIVATE studies, not accepted production corrections. No library fork changed.
Detailed original gates/results remain in
[the archived checkpoint](experts/history/checkpoint-handoff-2026-10-27-transfer-screening.md)
and [Capytaine diagnostic record](experts/capytaine_expert/issue31-diagnostic-plan-2026-10-27.md).
Do not delete/relabel failed studies when adopting the user height requirement.

Mixed analytical-Rankine/numerical-wave integration and a separate native
precision/approximation path have been identified. Neither a Green split nor GL2
alone is qualified. Native point precision retains both E1 and angular arithmetic
limitations; do not blame E1 casting alone.
The double-precision fixed-point witness passes its bounded checks (review52ddd89f).
The adaptive SciPy Heave study passes ALL field/reference controls but fails the
finite-band cosine recovery proxy (reviewb63d3a9d, report98036466,158.562s).
That is not a sourceK-causality proof or an actual height error. Its precise
frequency/derivative gates were diagnostic engineering choices, not the user's
end requirement. Do not launch more sourceK or whole-model investigations merely
to perfect those proxies; first measure actual-height error and reference limits.

Previous reference-adoption and cleanup goal is COMPLETE at55e44a5. Its evidence
below remains the immutable starting point, not authorization to repeat its runs.

Merged fork reference and DEFAULT runtime pin are now
6b774228e1f6b31dee473e626204b0ef65e5ec38 over official v7.1.0/a23d299.
The user merged PR6/7/9/14; keep these and closed fix issues4/5/3/12.
The standard953 adapter is integrated in the active main workspace. Four physical
handoff files are exact953 after CRLF normalization; backend differs only by the
pin. Held-pin regression now explicitly selects995. Default checkout
artifacts/dependencies/topwater-wec-sim is clean trackedtopwater/reference6b.
No upstream submission or assistant merge occurred.

## Execution and accepted adoption

NO LIVE NUMERICAL JOB. All latest producer/encoder/owned MATLAB PIDs are terminal,
frozen-input checks and separate idle/readback/release-last passed; lease absent.
Last numerical job: issue31-bandwidth-attempt1, terminal0 (diagnostic only). Recorded status
is a checkpoint, not current liveness proof; revalidate before any future job.
No repeat body integration is needed to investigate #31.

Independent accepted receipts under artifacts/workflow:
- reference-adoption-motion-review-v1/acceptance.json286f3a82: NEW6b+953 full3001x6
  forces/startup/phase/actualonce-only application/rest/negative-history/fullCIC/
  clocks/motion. Max difference59.758um/.00156985deg (limits1mm/.01deg).
  Reportcase3ec4d8d8d, freeze756cc7cf. Whole1212.407s; simulation745.50s,
  nativeforce305.10s, init49.09s/model51.77s. No physical-accuracy claim.
- reference-adoption-private-review-v1/acceptance.json4622b301: fresh TOTALfield,
  standard Result/PNG/logs and actual ordinaryCLI render, exact8files,
  3001motion/751surface/751FULLdecodedvideo. Finalreport2104ca1f, assembly1f16fa47,
  sourceb078f2fe, fieldfreeze53be5d65/renderd8ecf8b8/final8f28fef0.
  Whole field987s, CLI render126.906s, audit26.032s. No extra body integration or
  coefficient campaign. Fresh-field source raw953CRLF528e/LF240 both verified.
- reference-adoption-main-review-v1/acceptance.json3712ab33: actual MAIN38tests
  pass2.69s plus DEFAULT6b/nooverride/owned hiddenR2025b actualpreflight, all three
  licenses, no loaded models, server exit0/heldnativehandle/shared-name gone.
  Actualbackend712f0547; preflightfreeze172fb5a6; helper32.563s. The receipt binds
  tested working bytes; the completion audit verifies committed source identity.
- Earlier38Python+21native handoff0cffbbc3, generic10spectrum4ffbf7f8, tiny3public
  casesb787454d remain supporting selected-path evidence.

Guard-v2 d168ff67 measuredcc4ff62f passed9realWindows controls. Live optional
publication defers only Win5/32/33 without retry sleeps; resource controls and
mandatory initial/final receipts stay active/fail closed. Keep oldv1/v4 and
failed attempts unchanged. Case1 sandboxprefs failed beforephysics; case2 guard
publication interrupted motion; revised memory-control test validates actual
128MiB denial64MiB/success256MiB without asserting undocumented peak semantics.

## Canonical output and source state

The reviewed eight private files are PROMOTED byte-for-byte to examples/simple-boat:
boat.npz/sea.npz/result.npz/result.png/video.mp4/sea.log/result.log/video.log.
Promotion proof reference-adoption-promotion1/report.json binds private462/main3712.
Boat4168e1a7/sea9ba25ba3 unchanged; newresult2b615b56/PNG275ac4c5/videoe5fce0e4.
Original8files343MB retained and hash-verified in reference-adoption-public-rollback1
against the original field freeze. Public logs faithfully retain private generation
paths and explicitly resumed accepted motion; do not claim another public solveCLI.
Actual qualification was native body + fresh production field + ordinaryCLIrender
+ actual adopted-main defaultpreflight. No agent image views/previews occurred.

All three logsENDsuccess; sea/video0WARNING/result34 knownrows (including wrapped
continuations): restoring6/tabulation-info1/projection1/mesh17/irregular8/wetset1.
Native batch one endpoint-coverage heuristic warning plus two INFO messages;
all consumed headings supplied, force replay exact within fixed tolerances.
Genuine15deg BEM angular accuracy remains unqualified. T-I RMS5.98768mm=2.24251%,
onebeam7.99367mm=2.96028%, maximum65.0892mm at44.56s. These are field disturbance,
not error bounds. Old/newXY differ; no indexwise subtraction/regridding.
See docs/cli/warning-triage.md for current complete dispositions/physical31-34.

Current code uses existing bodyClass(hydroStruct), official readCAPYTAINE /
radiationIRF/rebuildHydroStruct. WEC computes wave forces and motion; no custom
component-view/accessor/BEMIO writer extension. Spectral H5 skipped, consumedCIC/
Ainf preserved, unused rawA/B/T omitted. Native Topwater remains an incomplete
prototype. No remeshing/incident-only fallback/nonlinear-physics claim.

## Completed fork cleanup

Adoption checkpoint `7e1b854` is pushed. Technical adoption and cleanup passed
independent current-state audit3968adb4. Final publication and remote verification
are recorded in `reference-adoption-completion-audit1/report.json`; consult that
receipt before claiming whole-goal completion.
See [completion record](experts/wecsim_expert/fork-cleanup-plan-2026-10-27.md).

Only remote branches `main`0753b2e (unchanged upstream mirror) and
`topwater/reference`6b remain. Nine obsolete heads were removed with expected-SHA
leases. Verified annotated tag `archive/topwater-legacy-9952634` preserves exact995;
the old runtime branch is removed and official v7.1.0 tag is unchanged.
PR8/10 were natively archived; obsolete issues1/2/11/13 were deleted after
preservation. Merged PR6/7/9/14 and closed fix issues3/4/5/12 remain.

Actions is disabled ONLY in our WEC-Sim fork, preventing inherited documentation
deployment from recreating gh-pages. Topwater Actions remains enabled. Review
documentation deployment before deliberately re-enabling fork CI; current checks
were run locally. Upstream and mirror source are unchanged.

Eleven obsolete local worktrees and thirteen local trial/legacy heads were removed.
Keep default6b, stocka23, the clean detacheda23 common parent, linked6b validation
and qualified953 Topwater source. Their shared Git store remains intact. Other
Topwater worktrees and original failed numerical studies are preserved.

Recovery evidence: `artifacts/workflow/wecsim-cleanup-preservation1` contains the
184.6MB bundle d27eaf01, original shallow boundary, six dirty/untracked source
snapshots and24 ignored Simulink files. Independent offline recovery WITH the
shallow boundary restored24refs/all four unpublished trials and passed connectivity.
A plain clone without that boundary failed and is not the recovery method.
Historical issue text and redirected links are retained in the committed
retired-workflow record. Execution evidence is in `wecsim-cleanup-execution1`;
independent audit is `reference-adoption-completion-review-v1/audit.json`.
Final publication is tracked by `reference-adoption-completion-audit1/report.json`.
Fork cleanup is complete; #31 work is controlled by the current scope above.

## Persistent limits and visuals

Physical issues31-34 remain open; compatibility is not physical accuracy. Only
#31 is currently authorized; do not expand into the deferred investigations. Keep warnings visible/triaged.
No Topwater-specific WEC capability without a demonstrated standard-interface
blocker and strong engineering rationale. Old expanded e24 route acceptance76f5616a
is historical only. See experts/history/checkpoint-handoff-2026-10-27-pre-adoption.md.
Visualization141046c5 (50 controls), plotf39a7572 (23 controls) remain unchanged:
water/fullgraymesh/blackoutline/solidbluecontact plus unlabeled fixedgraymeanoutline;
seven nice symmetric y ticks, centeredzero/common translation/rotation row limits.
