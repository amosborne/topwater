# Checkpoint handoff

Updated 2026-10-27 UTC. Follow AGENTS.md, coordination and shared contract.
Latest user instructions control scope; historical assignments are not permission.

## Active goal and authorization

Current goal: resolve Topwater issue #31, the causal radiation-surface
amplitude/phase discrepancy in `topwater solve --solver=wecsim`.
User priority: #31 first, then measured end-to-end performance, then #32-34.
Upstream contributions are deferred. No upstream contact/issues/PRs or merges.
Active implementation branch remains `baseline/fluid-body`.

Start with bounded saved-transfer and analytic controls, one mode and near/far
points. Establish cause and prospective acceptance gates before correcting the
method. Preserve accepted boat motion/forcing, incident/diffraction physics,
physical inputs and old studies. Numerical work is serial under the coordinator's
lease. Authors do not independently approve consequential physical corrections.
Full reference-convergence/applicability campaigns (#32-34), native dynamics and
optimization are not current assignments; targeted discriminators needed for #31
may be scheduled explicitly. No numerical job is currently allocated.

First saved-data discriminator completed in4.594s under guard120s/2GiB;
4385frozen files, terminal0/Jobempty/hostCIM/frozen/readback/release-last passed.
Report: artifacts/workflow/issue31-stored-diagnostic1/results/report.json.
Near/far Heave complex recovery remains about31.3%/3.48% at cutoff24, lag60;
lag100 barely changes it. Matching-endpoint analytic controls are far closer.
Independent diagnostic review c6dd7158 accepted the interpretation only.
This is diagnosis, not production qualification. Do not fit a replacement to these residuals.

The fixed17facet Heave Green-split comparison completed in0.5s (guard3.329s),
report `artifacts/workflow/issue31-green-diagnostic1/results3/report.json`
(9da95c7d), accepted diagnostic review a28c7cf7. Same-source field evaluation
changes remove most of the suspicious low-split high-frequency growth, but both
splits retain substantial free-surface boundary residuals. Neither is a selected
correction. Two earlier helper/API failures are preserved; original helpersv1/v2
and all attempts remain. Successful helper SHAee80509d avoids cached LU wrappers
and normalizes endpoint sources before numeric gradient stacking. No Capytaine
library or production code was changed.
The same17facet GL2/centroid comparison completed in0.672s (guard2.547s),
report `artifacts/workflow/issue31-quadrature-diagnostic1/results/report.json`
(2b9e3ce4), accepted diagnostic review85f85a18. GL2 reduces some raw residuals
but does not qualify either split. Potential finite-difference witnesses expose
an additional potential/gradient inconsistency in the high split; low reported
boundary residuals alone are not proof of correctness.
The extracted-Rankine mixed-integration identity predicts the observed low-split
residual proportional to omega^2 and high-split residual independent of frequency.
Kernel identity diagnostic completed in0.265s (guard1.891s), report
`artifacts/workflow/issue31-kernel-diagnostic1/results/report.json` (185a1cbd),
accepted independent review6e4ac9de. The extracted-Rankine identities predict
raw residuals within1.3e-7absolute on fixedpanel0; directtheta1001/2001 evaluation
also identifies a separate table-interpolation inconsistency. This is a cause
finding, not a qualified correction or physical error bound.

Current author assignment: private FIRST candidate only, no production change.
Use common full-Green panel quadrature for exterior potential AND gradient by
subtracting analytical Rankine terms and adding the same numerical quadrature
terms. Keep original bodyK/source solves and production cosine inversion fixed.
Dense .025..24rad/s Heave, near[0,1.25]/far[6,0], six original heldouts; directtheta
1001 training versus2001 reference, productionlag60 gate (lag100 diagnostic only).
Use per-point S=maxheldout|Hreference|; complex1%S at every heldout, amplitude2%
and phase1degree where response>=5%S; S=0 uses1e-12seconds absolute and undefined
phase. Reference theta convergence target0.1%S. No fitting or source/model tuning.
Source/API review and explicit root allocation precede execution. If oldsourceK
blocks the causal gate, stop and diagnose source-operator consistency rather than
adjust the kernel to hide it. Six-mode and full saved-motion qualification follow
only an informative first candidate. Canonical outputs remain unchanged.

All owned processes are terminal; frozen-input/Jobempty/hostCIM/idle/readback/
release-last passed. No numerical allocation is currently live.

Coordinator adopts prospective candidate gates in
experts/numerical_methods_researcher/issue31-method-gates-2026-10-27.md BEFORE
any candidate validation: per-point/per-mode complex1%S, amplitude2%/phase1degree
where response>=5%S, explicit zero handling, allsix/heldout/nearfar/waterline,
causal analytic bursts with1mm and1%peak criteria and independent references.
These are engineering reconstruction tolerances, not physical accuracy bounds;
they deliberately distinguish the current material discrepancy while preserving
small-response safeguards. Exact transient error requires an independent causal
reference. Actual BEM transient sensitivities cannot substitute for that proof.

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
Last numerical job: issue31-kernel-attempt1, terminal0. Recorded status
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
