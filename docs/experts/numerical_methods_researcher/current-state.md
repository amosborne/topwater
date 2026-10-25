# Current state: independent numerical-methods review

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
