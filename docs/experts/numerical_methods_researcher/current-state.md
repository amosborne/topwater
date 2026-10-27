# Current state: independent numerical-methods review

Final integration addendum: the promoted policy passed50 selected integration
tests in24.20s with no skips; independent receipt review accepted exact coverage
and unchanged numerical module AST apart from the readiness Boolean.

## Accepted selected-case reconstruction, 2026-10-27

The simple-boat case passes the user's **10 mm maximum elevation difference
including numerical reference uncertainty**, against the chosen linearized
model. All **5,929 targets at 751 times** were checked. The largest combined
estimate is **5.7844 mm**; the actual runtime/reference difference is at most
**0.3412 mm**, with maximum time-weighted RMS **0.08965 mm**. These are numerical
reconstruction quantities, not nonlinear or experimental error estimates.

The [selected-case acceptance note](issue31-selected-case-acceptance-2026-10-27.md)
links the exact reports and hashes. The independent root-executed saved-array
audit checked every combined value, all 47 chunks, common incident/diffraction
identities, unchanged motion and the complete exported water geometry. All 751
existing physical geometry/contact/startup frames passed. This reviewer authored
the independent validator, not the physical correction, and launched no solver
or numerical job.

The reference remains conditional on the explicitly accepted, location-specific
regular-response envelope: twice the maximum calibrated from 36--48 and 64/96
rad/s, with known-32/prior-probe checks, bounds uncomputed frequencies. This is
an engineering assumption excluding unseen larger resonances or faster growth,
not a theorem or a mesh/model-resolution qualification. The known Hermite input
tail is enclosed independently; no uniform 0.5 mm allowance is used.

The fixed-mean wetted-area warning remains: worst symmetric difference is
18.27% at 44.64 s. Retain it as a model-applicability limitation requiring later
evaluation; it does not invalidate reconstruction accuracy against that same
linearized model. Broader model questions in issues #32--34 remain separate.

Canonical Result, motion plot, logs and the regenerated video are published.
The [video log](../../../examples/simple-boat/video.log) records decoded
verification of 751 frames, 12.5 fps, 800 by 800 pixels and the complete 60 s
source interval. Readiness was enabled after cleanup by changing only its
Boolean; the coordinator verified numerical-function ASTs were unchanged.
Final integration checks pass50 tests with no skips. Checkpoint push and issue
closure remain coordinator follow-up, recorded in the handoff and GitHub issue.

## Prior investigation summaries

The following are historical intermediate findings, superseded for current
acceptance by the note above. Failed studies and their original policies remain
preserved; statements that qualification was pending describe those checkpoints.

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

The private full-complex positive-lag candidate was tested and **rejected**:
[causal screen review](../../../artifacts/workflow/issue31-causal-review1/acceptance.json).
It removes a 194.93 mm saved-motion contribution despite passing exact delay,
fold and direct-restoration controls. Its band sensitivity is 11.13 mm. Merely
keeping positive lag is therefore unsuitable here. The analytical discrete
acceleration endpoint is numerically solvable; replacing spectral acceleration
with saved acceleration changes the full reference by only 0.287 mm. Sparse
high-frequency data do not establish asymptotic decay or no-growth above 48.

Recommend preserving full-complex response for offline reconstruction, with
explicit startup, future continuation and numerical uncertainty. A universal
runtime causal kernel is not an additional user requirement. However, the
chosen linear equations are causal: frequency-domain synthesis alone does not
prove the rest-start physical solution, and the roughly 3 mm future controls
are selected-case sensitivities rather than a universal anti-causal error bound.
The remaining-tail allocation stays conditional. No correction or unconditional
reference truth is accepted. Next implementation and total-field qualification
belong to the coordinator and solver owner; no further kernel tuning, BEM or
mesh sweep is recommended by this review.

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
