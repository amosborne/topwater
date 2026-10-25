# Current state: WEC-Sim / BEMIO expert

Updated2026-10-24. Scope and authorization: [handoff](../../checkpoint-handoff.md).
No active numerical allocation; latest reference/video pipeline is complete.
No automatic follow-up is authorized. Coordinator owns
requirements, acceptance, integration and serial execution allocation.

## Current priority

Shared Matplotlib2Dtopdown proof is COMPLETE: ROOTpeakposterPASS/USERlinked.
[Watch5s topdown](../../../artifacts/reference/stronger-buoy-20261024/matplotlib-clip-attempt1.mp4).
SHA954343c52b6bc5690fcf421bc0ac6804c2a4643b2571c26ffc34682f44e87eed;
800square/H264/63frames12.5fps/5.04s, savedtimes50.00..54.96. All17tests/full63decode
PASS. Render11.25s, wholejob17.797s; MatplotlibAgg blit/directstream, no GPUneeded.
[Final2Devidence](telemetry-2d-readiness-2026-10-24.md).

Geometry/time/connectivity+dtypeBYTEIDENTICAL tooriginalscene1079ae90. NewNPZ
2475d3be addsprovidedroll_deg/pitch_deg/heave_m fromEXACTacceptednativeqrows,
no interpolation/rezero/waveprep/headings/physics. GenericTopwaterCLI is
`uv run topwater render INPUT.npz OUTPUT.mp4`: two positionals, wholeNPZ/fixed800,
MP4only/defaultMatplotlibcore/no--extra. Renderer never infers WECpose.
ReferenceQAposters/proofs areexplicit, separate frompublicdefaultoutput.

No furtherWEC code/jobs. Executionowner visualization finalizesledger/release+
uvlock; ROOT handleshandoff/integration/checkpoint. Retain all3Dmovies/source/
caches pending laterUSERreview/conditionalcleanup authorization; no cleanupnow.

## Previous accepted PyVista preview

Reusable off-the-shelf PyVista/VTK preview is COMPLETE and ROOTposteraccepted.
[Watch5s square clip](../../../artifacts/reference/stronger-buoy-20261024/pyvista-clip-attempt3.mp4).
SHAc80e14001a9e3818291e428cea71378b3e158b44831d084afcc34b67d1ec2b0f;
800square/H264/63frames12.5fps/5.04s, physicaltimes50.00..54.96. Full63decodePASS.
ActualIntelArcPro/cachedrender31.515s; initialreusabledata preparation100.75s.
Executionowner visualization finalizes exactterminal/lease via sharedledger.
No further WEC code/jobs; ROOT handles checkpoint/integration.

ProducerWECadapter selectsacceptedposes+sameincidentfield into neutralSceneSequence;
genericTopwater library owns standardPyVista/VTK drawing and publicCLI:
`uv run topwater render INPUT.npz OUTPUT.mp4`, BOTHrequiredpositionals, WHOLE NPZ,
fixed800square/cadencefrominput, onlynamedMP4defaultoutput. No WECknowledge/hull
constants in library; convexmeshrequirement validated, differentbox smoke usedsameCLI.
PreparedsceneSHA1079ae90 reused forcachedlighting; no wave/physics/mesh changes.
[Workflow/evidence](pyvista-readiness-2026-10-24.md). OptionalVTKshadow/waterline
artifacts were rejected/preserved; finalflagsFalse/standardcontinuouspalette.
Visualconvexclosure/prescribedwater is not flooding or nonlinearaccuracy proof.

## Previous completed stronger presentation

USER requested clearer existing stronger60s motion; FINAL PRESENTATION COMPLETE.
[Watch artifact](../../../artifacts/reference/stronger-buoy-20261024/presentation-encode-attempt2/reference-presentation.mp4).
MP4SHA7a9debb31a440b6967fd3340f2fcd6f89072699dd7ecc2e36136d3bca6079aa8,
H2641280x900/751frames12.5fps/60.08s. All751decode andROOTliteral0/30/51.44/60
PASS; userlinked. No activejob; exact14472/renderer/ffmpegterminal, all24frozen
inputs/depsmatch, FINALrelease8bbf9c. ROOT owns final checkpoint/integration.

24deg obliquebroadside/10m/equalmetres, noaxes/grid/wiremesh; compactHUD,
solidinside/outside andactualcached-contact waterline. Distantwater320x225 display
pass, hull/contact/water-neighborhood1280. Finitepatch edges and coarsening remain
visible limitations; no physicalfield/motion amplification or solver changes.
Actualencode1064s after selectedpreview5.469s. All source/cache/input/QA evidence:
[Presentation receipt](presentation-readiness-2026-10-24.md).

The actual stronger record46086/2023-12-23T00:40Z Hm01.65747/Tm027.23657,
seed42/384unscaledcarriers, bow182.067deg, spatialtrain tc40/plateau15/support25.
Native original17H5/CIC10/dt.02/all6free/qv0/memory0/ramp0/visc0; actual60s
implementation/energy independentlyPASS. Peaks[.562,.837,.759]m and
[6.618565,2.677850,1.325248]deg; netradiation32971.040J, temporaryreturn4074.896J,
cumulativemin0; incident46582.193J. These are trajectory diagnostics, not full
fluid-energy closure/passivity or nonlinear/mesh/time accuracy qualification.
[Case evidence](stronger-buoy-readiness-2026-10-24.md),
[original stronger movie](stronger-video-readiness-2026-10-24.md).

The cached localwaterplane proxy has maxrelativeheel3.10deg/trim2.88deg/
relheave.09m; this does not certify combinedwavewetting/flooding/linearity.
Geometry applicability remains approximate. No furthercases/physics/mesh/bug
studies or Topwater solver work are authorized. Preserve originalstronger and
mildvideos (mildpeakroll1.481deg), rejected tiny-cap inputs and superseded
stationary-history prototype under their original interpretations.

## Latest executed official controls

TWO unchanged official v7.1.0 supplied baselines completed on maintained fork
f962103a348d63101ce669c2611f7d7c109755c3 (officialbasea23d299), clean checkout.
Standard RM3 regular400s/4001samples: PASS31.257s, no suppliedgolden. Official
regularCIC200s/2001samples/CIC20: PASS28.280s; fullgolden float/spar/relativeheave
errors<=4.16e-13m vs official1e-5m tolerance. Independent implementation/integrity
checks pass. Native initialv/Frad0 and fullgrids/forcesums pass; all134/138frozen
identities unchanged. Their official100s ramp/PTO1.2e6Ns/m/constraints are retained;
these are openforced systems, not our freeboat or localcalm-arrival startup.
[Source, exact output hashes and executed receipts](official-examples-readiness-2026-10-24.md).
Both exact processes terminal, separate idle proofs/final releases complete.

## Preserved physical evidence and limits

REF001 minimum10s/501native samples is COMPLETE/frozen with actual-rest/calmtail,
fullmass/frame/force checks and centered equal-scale finalvideo. Internal native
consistency accepted; mesh/time accuracy unestablished.
[Minimum reference](minimum-reference-2026-10-24.md).
Latest stored10s energy check: cumulative resisting CIC work nevernegative at
501samples, total+1.253275845J; instantaneouspower reaches-.0345015W. Temporary
energy return is not itself nonpassivity. Added-mass mass-split work is separate;
this trajectory does not prove coefficient passivity. Original rawSHA80a16af1?
and motion78c9a354? match independent prior review.

REF002 old full384component buoy banks, unchanged17panel exports, officialH5/CIC40
and harmonic-only screens remain immutable historical evidence.
[Consumer record](generated-bank-consumer-execution-2026-10-24.md).
Old Jan90 held-release case is REJECTED (107degroll/free-mode drift); full60s saved.
[Case](jan90-native-case-proposal-2026-10-24.md).
Small calm0.5degroll control decayed to.00004258deg in60s; this does not establish
all-state stability/passivity. [Control](free-decay-control-proposal-2026-10-24.md).
Nonzero particular q/v+CIC-history prototype and its native fixtures passed
implementation checks but are SUPERSEDED for final requested zero-memory/localcalm
startup. Its5.127deg peak occurred only at final60s sample, not sustainedstatistics.
[Preserved prototype](seeded-native-case-proposal-2026-10-24.md).
No further history adapter or seeded trajectory work is active.

Full native17panel B has negative harmonic-power directions; Cap owns direct
pressure-work/Kochin witnesses and producer interpretation. Consumer signs,
convolution consumption and nativewrench checks do not establish an implementation
bug or repair unreliable BEM coefficients. No coefficient symmetrization/damping
fit/blanket accuracy claim is accepted. Geometry/mesh uncertainty remains.

Topwater solver/nonlinear adapter/custom integrator work remains HELD. Only the
existing sea generator is available for later separately scoped final input work.
Execution ledger owns liveness; engine9168 is historical identity, not liveness
proof. Reconcile authoritative outside-sandbox IPC and exact processes before
numerical launch; restricted process visibility never justifies engine replacement.
