# Current state: independent numerical-methods review

## Current issue 31 review, 2026-10-27

The user requires **at most 10 mm numerical surface-height error against the
chosen linearized model**. Nonlinear and real-world accuracy are separate.
The coordinator owns execution and acceptance; this reviewer has authored no
physical correction and run no numerical job. Current authorization and process
state belong in the [handoff](../../checkpoint-handoff.md) and execution ledger.

The accepted WEC-Sim motion and forcing remain unchanged. The current review
concerns surface postprocessing. Prior transfer-response and kernel diagnostics
are useful explanations, but they do not replace the user's height requirement.

The all-six-mode saved-motion height screen uses the new accepted 0--60 s
trajectory at all 3,001 native times and seven cached transfer probes. Compared
with the full-complex response candidate, the current cosine reconstruction
differs by **42.66 mm** at an exterior continuation trace, **31.93 mm** at a
direct near-hull probe and **3.47 mm** at the far probe. These are actual
trajectory-dependent method differences, not transfer-normalized percentages.
The [independent height-screen review](../../../artifacts/workflow/issue31-height-review1/acceptance.json),
SHA `9d55d233`, binds the input, source, numerical report and execution evidence.

Exact FFT periods, time-step and input-interpolation controls, pre-response, and
smooth future continuations after 20, 40 and 60 s were checked. Estimated
reference sensitivity is 0.29--2.92 mm, dominated by future dependence. It is an
engineering sensitivity estimate, not a certified error bound. The 2 mm reserve
is planning guidance, not another user requirement. Missing response above
24 rad/s and incomplete spatial coverage prevent a qualified full-surface claim.
Four probes are explicit exterior traces for mean-interior continuation targets.

A known-causal point-Green control shows that a finite-band phase discrepancy
can occur without any body-source operator: the shallow source's 24 rad/s
projection differs from its full response by up to 1.31% of its reference peak.
This does not establish the cause of the entire body-radiation discrepancy.
See the [bandwidth-control review](../../../artifacts/workflow/issue31-bandwidth-review1/acceptance.json).
Separate studies demonstrated mixed Rankine/wave quadrature and native
evaluation-precision limitations; no production correction is accepted from
those studies. See the [adaptive field review](../../../artifacts/workflow/issue31-scipy-review1/acceptance.json).

The same-operator 24--48 rad/s extension is complete and independently reviewed:
[tail review](../../../artifacts/workflow/issue31-tail-review1/acceptance.json).
Overlap at 20, 22 and 24 rad/s is exactly identical for every mode and probe;
all old samples were preserved. Extending the actual-motion response changes
height by at most **0.0592 mm**. The absolute, uncancelled 36--48 band envelope
is at most **0.0387 mm**. This supports selected-case bandwidth convergence;
it does not bound the uncomputed response above 48 rad/s. The existing observed
reference sensitivity sum reaches 2.868 mm and remains dominated by future
prefix dependence. A proposed 0.5 mm remaining-tail allocation is an explicit
engineering assumption, not a measured error bound or accepted requirement.

Next, compare a cached full-complex signed-lag causal candidate directly with
actual saved-motion heights. Separate analytically derived instantaneous
acceleration before inversion; do not fold negative lag forward or infer a
causal-error bound from the future-continuation controls. The finite-low discrete
operator endpoint differs from native symbolic infinity, which changes Green
split. Its proposed feedthrough is a discrete-model endpoint, not added physical
fluid dynamics. No physical correction or full reference qualification has been
accepted. Seven probes remain a screen; final acceptance needs the declared
valid full surface and all saved times, with common incident/diffraction terms
preserved. No further BEM or mesh sweep is recommended by this review.

## Historical provisional policy, 2026-10-25

Updated 2026-10-25: coordinator accepted a provisional computed linear total
surface for the selected stronger-buoy 60 s motion on the unchanged 17-panel
mean hull: radiation cutoff 24 rad/s, spacing 0.025 rad/s, full-history causal
kernel, saved modal velocities, zero prehistory/direct acceleration term, exact
sideband diffraction, and explicit exterior-trace continuation. Body motion and
forcing remain unchanged. [Independent review](surface-reconstruction-review-2026-10-25.md).
[Response evidence](../../../artifacts/workflow/surface-response-attempt1/field/report.json)
has spacing sensitivity <=3.26012e-6 m and cutoff sensitivity <=0.000753226 m;
[phase/projection evidence](../../../artifacts/workflow/surface-projection-attempt1/report.json)
has raw-versus-causal differences up to 0.0426650 m (RMS up to 0.0125110 m),
with a noncausal raw contribution up to 0.0772652 m. These are material numerical
sensitivities, not physical-error bounds or convergence/accuracy certification.
Mean-hull continuation, exact-field versus prepared-load discrepancy and sampled
startup limits remain explicit. Complete artifact/initial-state gates are still
required; no broader campaign or native dynamics is authorized by this review.
Reviewer inspected source/reports only and authored no physical correction or
numerical job. The older notes below retain their original historical scope.

Latest user steering supersedes old queued follow-ups: finalcase60s, rest/zero
velocity andzero radiation memory, viscousOFF, physics-based incident-band policy.
Old unpruned/history-seeded case remains evidence, not the requested finalstartup.
FullCap radiation+diffraction/FK+diffraction applies; projectlambda2beam filter is
alreadyenforced. No Cap/WEC radiationbug or newlyqualifiedcutoff established.
See currenthandoff; no source role implies numerical execution permission.

Updated2026-10-24. Authorization/currentnumericalowner belong to the
[handoff](../../checkpoint-handoff.md) and execution ledger. This role reviews;
no implementation, physical correction or numerical simulation was authored.

REF-001 minimum native motion and centered video are independently audited.
Historical reviews: [hydrostatics](hydrostatics-review-2026-10-24.md),
[startup](startup-review-2026-10-24.md),
[mask dispatch](excitation-mask-review-2026-10-24.md) and
[incident packet](spectral-packet-proposal-2026-10-24.md).

REF-002 targets ONE genuine buoy-generator case with appreciable all6DOF,
roll about5-10deg, and an a-priori geometry-based linear assumption. A canonical
suite and nonlinear trajectory adapter are not requested. Topwater solver held.
The original held-release Jan90 run is REJECTED; the unchanged-operator0.5deg
calm control decays. Confirmed coupled negative radiation-power directions remain;
no operator-passivity or numerical-accuracy acceptance follows from bounded motion.

The seeded Jan90 native60s case (motiond69652e3) passes independent complete3001
initialstate/excitation/CIC/world-wrench/Ar/grid checks. Peaks are0.368/1.549/1.520m
and5.127/1.723/2.200deg. Roll>=5deg occurs ONLY at final60sendpoint; no sustained
statistics claim. Native acceleration-delay startup remains unchanged/nonzero
residual reported. No exactnativeperiodicity or nonlinear fluid initialization.
[Actual case and review evidence](../wecsim_expert/seeded-native-case-proposal-2026-10-24.md).

Static geometry study is independently accepted as geometry, not nonlinear
trajectory proof: [review](geometry-operating-envelope-review-2026-10-24.md).
At relative5degheel restoring departs12.33%, tangent24.03%, wet-set15.87%;2.5deg
freeheel gives5.96/12.83/8.53%. Exploratorylabels are not adopted physical gates.
Scalar actual-wave audit covers ONLY123/3001poses: valid/incomplete, early maximum
wetturnover at least9.7087%, resolvedrimdry. Full60geometry remains unqualified.
Batched recovery source PASS61c3598f/core3b04e8f2, fixtures prepared/unexecuted;
1%area/10mmrim are reportingprecision only. All3001times or explicitUNKNOWN;
oldmatchedintervals must overlap, originalscalar source/results preserved.

Renderer preview2 independently literal/provenance PASS: source82083d9b,
10m equal-unit XYtracking/fixedglobalz, actual0/30/60 images andall3001plot,
all751displayframeclipguards/source-data/phase checks. Preflightba7d7f9f;
rootreceipt119b3893. Visualization acceptance is separate from physical geometry,
linearity, mesh/time accuracy and passivity. No encoder or further job inferred.
