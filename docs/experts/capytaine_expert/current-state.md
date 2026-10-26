# Current state: Capytaine expert

## Current reusable workflow and path audit — 2026-10-25

The public canonical workflow is delivered at checkpoint `bda901d`.
`artifacts/workflow/canonical-pipeline-attempt4/` records terminal success:
60 s WEC-Sim solve in2393.743s and render in226.932s,3001 motion samples,
751 surface frames and751 decoded video frames. Root verified the exact process
tree, frozen identities and final artifacts, then released the execution lease.
The final shared-axis `result.png` was backfilled from the same completed motion
without rerunning physics. Current scheduling belongs to the handoff/ledger.

Production preserves the producer's total-water grid and original connectivity;
it does not cut or remesh the water at the boat. Hull-interior fluid is invalid
implicitly, and the physical waterline is computed separately. The result
includes incident, diffracted and radiated elevation from the saved motion;
reconstruction does not change body forcing or trajectory.

This is a provisional linear reference, not physical accuracy qualification.
The selected field policy remains cutoff24rad/s, deltaomega<=.025rad/s,
full-history causal projection, zero direct terms and explicit exterior
disturbance traces. Near-field full-complex recovery has not passed:42.7mm
maximum/12.5mm RMS signed-lag sensitivity is not an error bound. Canonical
field/prepared-load differences reach47.5506N and90.0598Nm by component;
minimal-mesh/irregular-frequency and restoring/wet-set applicability concerns
remain open in issues31-34. See [actual warning triage](../../cli/warning-triage.md)
and [phase investigation](../numerical_methods_researcher/causal-phase-issue31-investigation-2026-10-25.md).
Native Topwater dynamics remain held; no speed work is allocated.

Latest user scope is rigorous correctness review of the actual/intended
Capytaine/BEMIO/WEC-Sim workflow, including supported inputs whose symmetry can
hide a defect. Never-reached library branches are outside scope. Read-only
Capytaine audit verified installed3.0.0 and current producer SHA256
`2858b847efea11a97c50bd9d12c744bf3069723e32e5fa401309c4e2ee0f352a`.
Native displacement/velocity/pressure/free-surface signs and one-time FK addition
are consistent; no new native Capytaine sign/factor defect was established.
Exact facet K/CB/volume and full inertia are the interface authority; legacy
native hydrostatic sidecars are overwritten deliberately before H5 use.

The audit identified a finite-data boundary hole: producer checks cover A/B/Ainf
and dense excitation H through3.4rad/s, while the imported three-heading
diffraction extends through12rad/s. BEMIO's existing maximum-residual assertions
are not explicit finite gates. BEMIO ownership is preparing a narrow guard and
single-NaN pre/post witness; this note claims no correction or executed test.
Global surface-policy acceptance is permission for provisional method use,
not universal certification; misleading full-complex gate wording and broader
coverage require reconciliation without arbitrary boat/case whitelists.

Bounded next controls are explicit finite/roundtrip checks with a single-NaN
injection and finite control; the existing stored-data phase discriminator;
and knot/between-knot/low-frequency field-load comparisons. These are proposals,
not new numerical allocations or convergence claims. Every fix PR requires a
minimal incorrect-before/corrected-after example. Root owns scheduling and
integration. This reviewer authored no physical correction and ran no jobs.

## Historical failure-recovery checkpoint

The first canonical workflow reached201/751 export frames and exited2 after
1989.285s. A KH.dat cleanup failure masked its primary geometry/metric error;
that masked attempt alone does not establish flooding or topology failure.
Completed success is recorded above; the failed evidence remains preserved.

Recovery caches completed fields/full motion privately in OS-temp before export
and retains frame/context evidence on failure. Geometry-only replay requires
no body integration or BEM. Physical/contact/startup gates were unchanged.
Independent review approved physical-hull snapshot compaction and safe invalid-
time diagnostics. Root's recovery/interface suite passed209tests/1optional skip
in13.05s under `artifacts/workflow/recovery-interface-tests-attempt1`.
At that checkpoint pipeline attempt2 was still pending; it is historical status,
not an instruction to retry or restore clipped water export.
See [versioned source/evidence note](surface-reconstruction-source-2026-10-25.md).

## Preserved minimum-reference evidence

Minimum17-panel producer and first whole10s/501sample native WEC-Sim motion are
independently reviewed and coordinator-accepted at dt0.02. Final render3 video
and plots passed independent literal encoded-frame audit. Authorization and execution
ownership belong to the [handoff](../../checkpoint-handoff.md) and shared ledger;
no numerical work is allocated to this role. Topwater solver/refinement remain held.

Physical canonical input remains unchanged, SHA256
ea18844904e03bbfd72abc3396240c97b1994b3584e36c8149f8ae5f60a99e4f.
Original unfinished exporter and all failed/rejected/native centroid outputs
are preserved; history and identities are in
[minimum producer evidence](minimum-producer-2026-10-24.md).

Accepted producer output: artifacts/reference/minimum-20261024/native-gl2.
Native Capytaine3.0.0 performed1086solves on17triangles in69.9801s; all14 original
identity/finite/metadata checks passed. Separate supported nativeGL2 hydrostatics
reused all frequency arrays unchanged with zero new BEM solves. All8 independent
analytic checks and external executed review passed; maxK discrepancy4.37e-11.
CG/CB/fullM/K/Ainf identities retained; body_name importer alias explicit.

Official BEMIO/H5 and actual runtime CG/CB/volume/fullK/Ainf readback pass.
Single phase conversion residual1.455e-11. Consumer FileSolid products-order
correction independently reviewed/executed PASS; fullM+A roundtrip~1.6e-12,
original asymmetric residual retained. Proper R2025b save migration/graphs/
ports/scripts, cold-start preparation and worldCG/STL/nativeRzRyRx/worldomega
mapping reviewed. [Joint interface review](consumer-interface-review-2026-10-24.md).
Consumer owns BEMIO/startup/engine corrections; no consumer physical edit authored
by this reviewer. Startup/dispatch corrections have separate independent reviews.

Accepted motion: artifacts/reference/minimum-20261024/case-attempt1. Actualq0/v0/
radiation0 and signed prepared/applied501x6 forcing pass. Independent all501Solid
NewtonEuler balance<=1.71e-13N/2.14e-14Nm and exact forcebus sum establish internal
consistency. Initial apparent(M+A)*a mismatch is exactly stock1e-7TransportDelay
zero initial output plus forward extrapolation at dt0.02, recovered from stored
outputs<=1.18e-12; no frame/inertia/force-sign correction is justified.
[Native motion review](native-motion-review-2026-10-24.md).

Temporal and mesh motion accuracy remain unestablished.17panel frequency
resolution/irregular-frequency warnings remain limits. No dt study/refinement
allocated. Accepted video: case-attempt1/render-attempt3/reference-isometric.mp4,
literalFFmpegframes0/125/250 reviewed at0/5/10s; h2641280x900,251frames25fps,
10.04s endpoint hold documented. Nativepose/telemetry/centeredfullhull/equalmetre
isometric labels pass, exact sharedpacket and all5renderinputs unchanged.

Latest user steering: future higher-roll case must use the project's own
sea-state generator and generator dominant heading, preserving this baseline.
Bespoke period/heading RAO optimization is canceled; no such script or numerical
study was launched. A bounded a-priori static geometry operating-envelope study
is complete: samefullhull/nativeGL2, roll+/-12deg with freeheave/trim and
explicitrelative-immersion/trim poses, wet-set turnover/rim/contact flags, exact
buoyancy/righting vsfullmeanK. Attempt2 terminal0:15.218s/951geometryevals, all9
postcheckgroups and independent executed review PASS;49/49free roots,175static
poses (11rim-immersion flags).5/10deg RELATIVE calm-water heel gives righting
secant deviations12.33%/22.33%, tangent24.03%/39.32%, wet-set turnover15.87%/28.82%.
Stable/no-opening through12deg does not rescue the fixedmeanK assumption. The
exploratory5%/10% labels yield+/-1.5deg; not adopted finalcase gates or automatic
GLOBALroll ceiling. Wave-following requires actual generated-case relative-water
geometry evidence, not a blanket rescue. First failed abs-list attempt and BOM
claim-ledger protocol deviation preserved honestly; guarded retry used separate
verified claim/readback and exactchildwait, separate idle-underlease verification
before final release. No further numericalallocation, BEM or trajectory adapter.
[Geometry study](geometry-envelope-2026-10-24.md).

Next producer completed after independent/joint source review and sole allocation:
2838new native calls in67.119s terminal0, same17panels, .05:.05:12radiation/directAinf, realheading0
officialBEMIOcompatibility plus3x384exact untouched default-offset0 FK/diff/H pairs.
Paired MAT schema keeps native exp(-iomega*t) phase and full6DOF, no Cartesian
missing-data fabrication. Preliminary SCREEN-only A/B interpolation is explicit;
H is direct. All9 postchecks pass; independent FK max2.91e-11, phase reconstruction
5.82e-11, spatial frame closure2.22e-16. Sources frozen through exact terminal;
separate idle-under-lease readback passed and lease released. Executed independent
producer/interface review PASS; mesh/kernel/temporal/relative-wave applicability unassessed.
No further numerical allocation. [Executed producer](generated-bank-native-export-2026-10-24.md).

Stored native/CIC comparison exposes an unresolved radiation issue: heave
A(12)-Ainf=-8230.20kg, with substantial B(12); April particular-history difference
0.39567deg roll/0.05393m heave. This does not attribute error to cutoff, mesh or
irregular frequencies. One combined diagnostic source is independently reviewed
PASS (ff75f6bc):66 radiation scout calls across17/68 unchanged-surface panels plus
384 exact Jan05 existingoffset90 pairs, all raw coupling/phase/frame preserved.
No272 cascade, coefficient repair or canonical heading adoption. One allocated
450-call run completed56.9945s terminal0; all geometry/reuse/pair/hash checks pass,
lease released with separate idle-under-lease proof. Raw68-panel omega12 has
negative heave B=-48117.68kg/s; direct Ainf and highfinite probes show strong mesh/
asymptotic mismatch. Independent executed integrity review PASS; numerical
adequacy fails. Diagnostic only, no coefficient replacement.
No further numerical launch allocation.
[Diagnostic proposal](radiation-heading-scout-2026-10-24.md).

Read-only native infinity API audit finds no exporter misuse: omega=inf is
supported, k is symbolic infinity, free_surface remains0, and native symbolic
normalization produces Ainf. Finite default low_freq GF decomposition and native
infinite high_freq branch differ; coarse finite/direct discrepancy remains unresolved.
Same68 hull plus native interior lid24-call control completed55.103s terminal0,
independent executed integrity PASS;83 panels cover77.27% waterplane. Ainf unchanged
to2.39e-12; B12heave/pitch flip positive, but sway/yaw remain negative. Passive
operator NOT obtained; partial-lid sensitivity does not uniquely attribute error.
No coefficients replaced, no further numerical allocation. [API audit/control](infinite-frequency-lid-control-2026-10-24.md).

Earlier held-release Jan90 native run is REJECTED far outside intended operating envelope;
WEC owner/researcher diagnose the stored actual dynamics. Geometry audit must not
guess cause. Source-only actual-wave geometry tool has independent source PASS
(runner5ea25c14/core7dda063b) with analytic wet/rim bounds, actual pose schema, explicit budget/coverage
unknowns and unexecuted fixtures. Its rejected motion remains excluded/unbound.
New seeded Jan90 actual3001-sample motion has provisional global roll5.127deg;
fresh audit config binds ONLY that motion (SHAd69652e3), rawbus/runtimeprovenance,
unchanged hull/bank/frame and nonzero seeded q0/v0. Algorithms unchanged. Bound
manifest232cee24 was checkpointed/allocated with READYstatus only. One actual
geometry job ended terminal0 at300s with INCOMPLETE coverage123/3001 through2.44s;
2878 later samples UNKNOWN. Fixtures allPASS;122wet/rimtargets met and rims dry,
lastrowcoarsewet/rimunknown. No full-case geometry/linearity acceptance. Old narrow
precision and original source preserved. ONE reportingprecision/batched clipping
recovery has source-only authorization below; no recovery numerical grant/job.
[Attempt1 evidence/recovery proposal](actual-wave-geometry-attempt1-2026-10-24.md).
Recovery SOURCEONLY: new batched module/
runner/config prepared, original scalar oracle and all inputs unchanged. Prospective
1%area/10mmrim reporting bounds,600s/1CPU/10Mtri/depth8; no physicalgate change.
AST only; fixtures unexecuted, independent batch source/math PASS; reporting-only
linear JSONL/checkpoint and strict pinned-prefix overlap change also reviewed PASS
at runner61c3598f/core3b04e8f2/configa411d7aa. No recovery
job allocated; renderer preview owns next numerical slot.
[Batched source proposal](actual-wave-geometry-batch-recovery-2026-10-24.md).
[Geometry source proposal](actual-wave-geometry-postprocess-proposal-2026-10-24.md).
Producer [README](../../../reference/capytaine/README.md).
