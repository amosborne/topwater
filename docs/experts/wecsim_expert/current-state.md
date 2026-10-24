# Current state: WEC-Sim / BEMIO expert

Updated2026-10-24. Scope and authorization: [handoff](../../checkpoint-handoff.md).
No active numerical allocation; source/readiness work only. Coordinator owns
requirements, acceptance, integration and serial execution allocation.

## Current priority

Prepare ONE mild single incoming physical packet diagnostic, full60s, before
returning to final buoy-sea selection. Water must be initially locally calm,
body q/v/radiationmemory0, all6free, viscousOFF; no stationary-sea held release,
nonzero history seed, hidden ramp, restraints, force fitting or cropped startup.
A finite-band packet has noncompact tails: quantify original REF001 numerical
calm gates over full hullregion and all6direct/retained initial forces BEFORE
any body run. No new incident/body numerical job granted yet.

[New mild packet source/readiness](mild-packet-diagnostic-readiness-2026-10-24.md):
prospective .015m/.25Hz/45deg/cos4 halfbandwidth.5rad/s/focus40s, unchanged reviewed
REF001 dispersive evaluator with only one nonzero slot. Full[-20,80] incident
padding yields[0,60] excitation. Existing17panel/CIC40 operator is UNQUALIFIED;
this simple diagnostic is not the final collected-buoy sea or a passivity remedy.
Source/freeze independentlyPASS; ONE incident preflight completed309.922s and
was REJECTED: fullregion eta2.069um/speed4.113um/s andarrival/bandpass, but initial
heave retained.129757N/spectral.147764N exceed fixed.1N cap. Otherfiveforceaxespass.
No body simulated, retiming/rerun/capchange. Completephysicalpacket/report preserved;
allinputs unchanged, exactprocessesterminal/separateidleproof/finalrelease complete.
Bodydraft SOURCEPASS but UNRUN/no bodymanifest orgrant. Coordinator now considers
SOURCEONLY extension of ACCEPTED originalREF001threepacketfield to60s, preserving
alloldcoefficients/initialcalmbounds/F0 andprovingfirst10excitationprefixidentity.
[Exact extension](reference-extension-readiness-2026-10-24.md) incident now independently
EXECUTEDPASS in.468s: allfield/waveprefix identities andfirst501forces EXACT0error;
originalfullplanformcalm/F0 inherited. Exactprocessesterminal/leasefree. BodySOURCEPASS
latest07e48e6f; originalH5/CIC10/fullM/nativechain frozen in newbodymanifest0766458d.
Native60sbody notrun/granted; independentmanifestreview pending. No newpacketfit.

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
