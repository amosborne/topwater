# Current state: independent numerical-methods review

## Current: approved chart and graticule video delivered, 2026-10-28

Final functional delivery is independently accepted. `result.png` contains the
six motion histories and a full-height original-Sea polar spectrum, supplied
explicitly through `save_solve_outputs(..., sea=sea)`. Its horizontal colorbar
shows original-carrier bin RMS elevation, `sqrt(0.5*sum(a_i**2))`, with linear
cividis colors over 0--0.12 m. Zero bins remain visible. Default incoming heading
points up; the dotted line shows actual initial heading. Frequency labels point
down, with no angular spokes, global title/footer or polar captions. The PNG is
3150 by 1170 pixels. No spectrum data is stored in Result or drawn in the video.

The user-approved video uses a full-bleed origin-centred square, metre edge
graduations, smaller inside header and a Z(m)/0 gauge. The gauge shares the
water's fixed coolwarm normalization, rounded outward to +/-1.6 m for this case,
and world-Y metres per pixel. All 1626 frames decoded at 800 by 800, 12.5 fps
and 130.08 s. Static transparent instruments are cached once and composited
above the changing scene before telemetry; real-buffer regression tests cover
orientation, transparency, draw order and changing poses.

The final seven-module cohort passed **187 tests in 17.22 s**. Delivery took
492.061 s including tests, PNG publication, ordinary CLI video rendering and
final decode/identity checks. Frozen source identities remained unchanged;
terminal and resource guard exited 0 with an empty job. Root confirmed all
recorded processes absent, finish session 97094 exit 0, and independent
idle/null/no-lease readback.

The final receipt is
[`actual-delivery-acceptance.json`](../../../artifacts/workflow/approved-layout-delivery-review1/actual-delivery-acceptance.json),
SHA256 `e967f076671cb97aee3364b4171a1bd37fba808ea9557f327a49b191aa8c6aa2`.
The eleven-array archive remains byte-for-byte the fresh native-unpadded Result,
SHA256 `aa050fa5280a4238ed85ef88dedec701da0469f7dcefd05bcdf8aa197e484e11`.
Motion, full water field, original inputs and physical warning history remain
unchanged by the final display work. No body, BEM, field or calm analysis was
rerun for this delivery.

The active physical reconstruction remains exact native 10 ms velocity FFT,
without oversampling, padding or invented future, using the exact mean-waterline
disturbance trace and pack-only export. Saved calm observations still show about
1.1 mm at the first/last frames but about 14 cm over complete rest intervals,
mostly present in incident water. They are not a wraparound proof or numerical
error certificate. The new 130 s numerical-height accuracy remains unqualified;
the historical 60 s conditional acceptance retains its original scope. Private
palette/log comparisons and superseded insets/styles/previews remain preserved.
No further numerical-method work is allocated; final metadata, issue publication
and Git checkpointing are coordinator-owned.

## Historical completed milestone: simplified surface and superseded video inset

The 130 s simple-boat case has been reconstructed afresh from the unchanged
13001-sample native WEC-Sim trajectory. No new body integration or MATLAB call
was made. The genuine rebuilt operator matched all 14 retained numeric native
arrays and all netCDF variables. Portable evidence is
[`simplified-surface-and-inset-run-2026-10-28.json`](../simplified-surface-and-inset-run-2026-10-28.json).

The current radiation policy uses the exact saved 10 ms velocities, with no
oversampling, zero padding or invented future. The 13001-point DFT has period
130.01 s. Full complex six-mode response, its phase convention, 0--48 rad/s
cutoff and real/imaginary transfer interpolation are retained. Incident and
mean-hull diffraction are still added separately. The producer source is
`4c60f9e370fb827fb036431bdc94325b27667612275aa6daaaf11503e6699506`;
the new completed-field cache is
`3814bf5cc928206f24e48a6356f359964c4dd71e390055273b61a068301f9ab4`.

Interior disturbance targets now use the exact nearest mean-waterline boundary
trace; actual exterior queries and incident evaluation remain at their true
locations. Strict interior Green queries remain prohibited. The small boundary
probe found finite values approached from outside. Its original 0.1% screen
failed and remains failed evidence; finite-offset spatial variation supports
engineering adoption of the exact trace, not a new accuracy certificate.

Normal solve and replay only validate and pack the full uncut water graph.
Moving waterlines belong to render. Opening-clearance, wet-set comparisons,
disturbance statistics and clipped total-startup analysis are removed from the
normal path by user instruction. Historical analysis remains explicitly optional.
The exact-trace/native-FFT cohort passed 57 tests in 4.34 s.

The user explicitly rejected a second padded reconstruction comparison. That
proposal and the earlier probe-policy proposals are superseded. Instead, a
separate saved-field diagnostic reports physical total and incident elevation
on exterior public-grid vertices during both complete rest intervals. For the
new unpadded field, maximum total height is 1.066 mm at startup and 1.059 mm at
the final frame, but reaches 138.844 mm in [0,10] s and 140.897 mm in [120,130] s.
Incident maxima are nearly the same. These measurements do not identify
wraparound causation or certify numerical error; quiet endpoints do not imply
quiet complete rest windows. The 1 cm operational calm criterion is distinct
from a 1 cm numerical-error allowance against the chosen linearized model.

The superseded delivered video included a fixed original-sea polar spectrum inset.
Three original carrier arrays were added to Result, including all 384 carriers;
all 11 prior NPY members and the solve-log prefix remained byte-identical.
WVHT displays 1.66 m. Ordinary render produced and decoded all 1626 frames at
800 by 800, 12.5 fps and 130.08 s encoded duration. Motion PNG and original
boat/sea/sea.log remain unchanged. The integrated data/render/CLI/surface cohort
passed 177 tests in 10.83 s; this overlaps earlier coverage and is not added to
it as a distinct-test count. No preview or image inspection was performed.

Independent actual delivery acceptance is
[`actual-delivery-acceptance.json`](../../../artifacts/workflow/sea-inset-delivery-review1/actual-delivery-acceptance.json),
SHA256 `e59baae834e3e24843a0d01e1a97d44759c4d2e6a51a2a45b301a86373ad9fa8`.
Root subsequently confirmed terminal success, frozen inputs, empty guard job,
all recorded processes absent, finish exit 0 and independent idle/null/no-lease
readback in
[`execution-closure.json`](../../../artifacts/workflow/sea-inset-delivery-review1/execution-closure.json).
No further physics or numerical work is allocated. The new 130 s unpadded
surface's numerical-height accuracy remains unqualified; the inset adds no
physical-result claim.

The historical 60 s selected-case acceptance remains valid only under its
original inputs, policies and conditional response-envelope assumptions:
5.784425 mm maximum combined difference and estimated reference uncertainty.
It does not transfer automatically to this 130 s input. Preserve that evidence
and all failed studies below under their original policies.

## Historical: fixed 10 ms cutoff-search implementation accepted

modulewave_packets.py a9b14b34..., tests12cd8243..., doc e1a41108... preappend
reviewed with10msdriver e01aa496. 21testspass0.10s; threeAPIcalls14.05/26.33/
61.82ms, cutoff.276923/.353846/1.253846Hz andcertifiedbounds<10mm. Solealltime
envelopegate; no dwellfit. Realharmonicincrementalsum,oneFFT+onefinalIFFT,
oddclock handling, norescale/noIFTscan. Guardadds0/0/1bins; firstcertifiednot
continuousglobalminimum. Rootjobfullycleaned. Evidencepacket-cutoff-search-
2026-10-27.json and docs/wave-packet-search.md. No nextjob/productionintegration.

## Historical: corrected real-Fourier outcome, 2026-10-27

Usercorrected Ncounting:DC+negativeconjugatesfree,Npositive. Priorcomplex
projection wrongforrequestedtask. Correctedsource/actualaccepted;9casesrootrun
6.785s,fullycleaned. Nooptimizer/rescaling. N40realpeak4.49/12.78/56.89mm;
quadraturevsFFTcoefs<=3.34e-14,ifft/direct<=9.26e-14,imag<=5.62e-16.
Plots wave-real-quintic-n-* and real-quintic-packets-2026-10-27.json arecurrent.
Noanalyticmagnitude-as-envelope claim, physicalqualitypass or nextallocation.

## Earlier complex-projection study: superseded and incorrect for this task

## Historical outcome: optimizer-free quintic examples, 2026-10-27

User requested N10/20/40 for0.01/.1/1Hz, phases0/45/90atMIDPOINT65, window
rest10/ramp5/dwell100/ramp5/rest10. Directquintic Fourier source+actualreview
accepted consistency; root ran9cases in6.8867s and fullycleanedserialslot.
Nooptimizer/LP/greedy method used. Source/plots/numbers in
[example note](../../research/pulse40-three-examples.md) and
[evidence](quintic-packets-2026-10-27.json). Positivecomplex-template projection
discards DC/negative content, not ordinaryrealwaveFFT. N40worstscore29.11/7.63/
5.82percent; lowf0loss ismethod-specific, notphysicalfundamental. Noqualitypass,
productionchange or nextjob. Preserve oldstudies below underoriginalpolicies.

## Historical intermediate proposals (superseded)


Current authorized example: rest10/ramp5/dwell100/ramp5/rest10,N40,phi0,
basefrequencies0.01/0.1/1Hz. Source pulse40-three-design1/demo.py independently
accepted abbba7b5014c9d4a65654270b0242a56de806ded3d171d8a19b726106b9df948,
plan.jsoncf24e4919d433e4b2cf65366dc3a334ca4418db9d58cf6b03f7ba284b29dd02e.
Finite positive dictionary greedyselection then angular minimaxLP; score includes
dwellcomplexerror,quietamp,rampovershootabove1. Measured continuous bound;
no1%gate/globaloptimality/physicalclaim. Rootsoleexecution pendingguardready;
no researcher numericaljob. Current state below is historical discussion.

Current accepted parameterization:f0,phi0 atdwellstart,R,D,Ncomponentsmaximum;
construction chooses frequencies/coefficients. User asks objective definition.
Read-only recommendation:minimaxnormalizedcomplex-envelopeerror forunitdwell
andzerocalm overfinitehorizon; freeboundedramps, documentedfrequency/normguards,
continuousvalidation. Notadopted/noimplementation/noallocation. Selectable
frequencies make general problemnonconvex; report bestvalidatedcandidate not
unproven globaloptimum. Prior bandwidth/binwidth input proposals superseded.

LATEST: each original wave component gets its own deterministic packet(f0,D),
not a mandatory shared envelope. Separate sea carrier grid0.01–1/.001Hz from
packet offset bandwidth/spacing. D100s shared; phases at commondwellstart;
preserve directions. Dcontrols spectral width scale, notminphysicalfrequency;
ramps control neededoffsetbandwidth; spacing controls uniformgridrecurrence;
ripple/overshoot need explicit constraints. New construction policy not chosen.
Priorrootjob completed restricted55s/12s in8.922s, allterminal/idle/releasepassed;
its edge-clearance premise is superseded, not a general packet requirement.
No new numerical job, production change or canonical replacement allocated.

Latest candidate is 0.01–1 Hz at 0.001 Hz spacing, revisiting the narrower
band below. Read-only review: 991 inclusive centres; exact harmonic grid repeats
at 1000 s; slots need not all carry energy. Final sidebands and boat-specific
short-wave eligibility still need checking. No implementation/numerical job.

Earlier user requirement fixes all final nonzero incoming components, including
packet sidebands, to0.05–0.5Hz. Runtime/canonical case remain unchanged. Earlier
minimax6s/K30 example reaches0.0104167Hz and is historical, not compliant with
the new input requirement. Generator metadata shows mostly0.01Hz source bands,
8stratified carriers per band (typical0.00125Hz strata), packet grid0.005rad/s
(~0.000796Hz). Retain that packet spacing as a starting investigation value;
source sampling density is not a convergence result. No new job is allocated.

User now permits discarding polynomial/smoothness constraints for minimal
spectral-signature research. DWELL semantics and specified carrier phases at
common dwell start remain the design directions. Top3 primary-source research
is in [the note](../../research/minimal-pulse-spectra-2026-10-27.md): minimax /
equiripple (first candidate), DPSS concentration, sparse frequency selection.
This role independently reviewed the separate researcher's note; no blocker.
Rectangle baseline jumps must lie inside transition gaps, not dwell endpoints.
No production policy or runtime savings is qualified or implemented.

User-requested minimax illustration is now executed and independently reviewed:
[comparison](../../research/minimax-packet-comparison-2026-10-27.md). SameK30
dictionary as quintic projection; R8/12s and10sDWELL. Sampled worst errors
0.42208/0.057321% vs Q1.28731/0.376701%; inverse4.44e-16m, work2.641s.
Both dense refinement flags remainFALSE (small residuals explicitly reported),
so this is an approximate grid minimax candidate, not continuum optimality or
physical acceptance. Between-node allowances retained. Source by this role,
reviewed by separate researcher, execution by root, all cleanup/release passed.

Root completed the latest phase-at-DWELL illustration: specifiedphi0 both gives
both1m crests at14.5s, unchanged ramps8/12s and10sDWELL. Six discrete amplitude /
phase panels use tau=t-14.5 atx0. New waveforms are intentional, old cases kept.
Global/dwell inverse errors4.44e-16m; selected direct sum2.00e-15m. Root also ran
actualsea period stats and cheap common-second ramp bandwidth screens. See
[the construction/evidence note](../../wave-packet-timing.md). All are bounded
temporal diagnostics, not propagated fields or whole-hull physical acceptance.
Sources reviewed by TOP; root alone executed; all serial cleanup/release passed.
Canonical inputs/results and physical solver code are unchanged. No next job
is allocated. Prior #31 acceptance below remains scoped to its original case.

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
