# Checkpoint handoff

Updated2026-10-24. Latest user steering controls scope. Read
[AGENTS.md](../AGENTS.md), [coordination](experts/coordination.md),
[physical contract](experts/shared-contract.md) and the relevant owner state.

## Current authorization and immediate priority

LATEST USER REGROUP (supersedes queued diagnostic work): stop chasing unproven
bugs. Select a genuinely LESS energetic collected-buoy sea using the existing
generator, return to locally calm water/boat at rest/zero memory, physically
incoming waves with NO temporal ramp, run the complete60s and PRODUCE THE VIDEO.
No official PTO-zero run, further mesh study or energy-bug investigation is
allocated. Preserve those source-ready proposals as UNRUN.

WEC owns arrival/consumer/native60s/video; Topwater input owner
`/root/calm_buoy_input` owns existing-generator observation/export ONLY;
researcher owns bounded independent physical-arrival/source review. Initial
selection proposal is a real milder observation (roughly0.3-0.6m Hm0, not a scaled
2.27m realization), seed42, bow aligned with dominant arriving energy. No phase
or amplitude fitting. Incident packetization broadens the stationary spectrum:
report actual finite-train energy/directional changes and support explicitly,
with source provenance and review; a stationary
global realization is not a calm incoming front. Source/readiness first, then
serial incident/native/video jobs. No new numerical grant currently.

WEC-Sim uses open-water frequency-domain hydrodynamics and body-only time
integration, not a finite simulated ocean patch. Define a physically consistent
wide incoming incident field; camera zoom is independent. Do not claim the
prescribed incident display is a computed scattered/radiated fluid domain.

Deliver ONE credible buoy-generator-derived native WEC-Sim case with appreciable
all6DOF and roll about5-10deg. User explicitly confirms **60s is the full validation
length**. Preserve the complete startup; no cropped window or canonical suite.

LATEST steering supersedes the unpruned stationary-history prototype:

- Prune incident frequencies that the actual physics/numerical procedure cannot
  support accurately. Establish the theoretical basis first; no arbitrary cutoff
  or tuning to roll amplitude. Source remains collected-buoy data. Retained energy,
  directional spread and phases must not be rescaled/fitted. Report removed energy.
- Startup is TRUE locally calm water with incoming waves initially upstream,
  boat at rest, zero velocity and
  **zero radiation memory**. The prior stationary q/v/history seed is a different
  initial-value problem and is NOT the final requested startup. Specify the
  incident field and initial pose consistently. Body rest in an established sea
  does not satisfy this requirement. No hidden force ramp; do not reintroduce a
  steady-history seed to remove physical free-mode transients.
- **Viscous damping stays off.** No added damping to suppress numerical problems.
  A later viscous model would need independent justification and matching physics
  in both WEC-Sim and Topwater, not reference-only tuning.
- No nonlinear trajectory model/adapter or Topwater solver implementation/tests/
  simulation/optimization is authorized. Existing-generator reference work is
  authorized. Assume a linear range a priori from geometry/restoring curves;
  geometry diagnostics do not prove nonlinear dynamic accuracy.

Underlying chain: Capytaine3.0.0 full linear radiation+diffraction, excitation
FK+diffraction; maintained WEC-Simv7.1.0 fork pinnedf962103a. It is NOT an FK-only
long-wave approximation. Official theory assumes inviscid/incompressible/
irrotational flow, wave amplitude small vs wavelength and body motion small vs
body dimensions. Project's wavelength>=2beam support envelope is retained, not
misrepresented as a universal diffraction-theory limit:9ftbeam givesfmax~.53346Hz
in deep water, tighter than existing.6Hz limit. Need separate physical envelope,
mesh/time/numerical coverage and radiation-kernel integration band.

No WEC-Sim or Capytaine radiation implementation bug is established. Raw coarse
native B already has negative cycle-averaged harmonic-power directions in the
full6DOF matrix; instantaneous energy return or negative offdiagonal/component
work can be physically correct. Verify sign/conjugate variables and a direct
native pressure-work/radiated-energy witness before assigning cause. Mesh/
irregular-frequency/formulation sensitivity is evidence of unreliable numerical
results, not automatically a library bug. WEC-Sim importing/CIC/wrench checks pass.

Latest user request: test assumptions against WEC-Sim's own provided examples.
COMPLETED two serial native runs on pinned fork, official physical inputs unchanged:
standard RM3 regular400s and provided regularCIC200s regression. Both PASS;
CIC full-history float/spar/relative-heave errors <=4.157e-13m versus official
1e-5m tolerance. Initial velocities/radiation forces zero. CIC instantaneous
resisting power can be negative while cumulative work stays nonnegative here.
This supports the exercised implementation, not our boat BEM accuracy/passivity.
Their native100s wave ramp/PTO/constraints are retained, not our final startup.
[Executed receipts](experts/wecsim_expert/official-examples-readiness-2026-10-24.md).

Producer direct pressure-work witness also executed:17-panel coupled2.2rad/s
outgoing mean power -447.744W per squared motion scale; same physical coupled
motion on68panels gives +1362.465W. Positive17-panel heave control passes.
Mesh sensitivity is established; no Capytaine bug or qualified replacement/cutoff.
Original10s coefficients contain the same offending matrix; actual10s cumulative
CIC resisting work remains nonnegative (+1.2533J final). No absolute Kochin power
normalization claim. [Witness](experts/capytaine_expert/native-energy-witness-proposal-2026-10-24.md).

NEXT bounded diagnostic: one mild physical incoming packet60s,
true locally calm startup, all6free, zero velocity/memory and viscous off,
using original packet machinery. New singlepacket incident attempt completed
309.922s and was REJECTED: eta2.069um<15um and speed4.113um/s<23.562um/s, but
initial heave retained/direct .129757/.147764N exceeds unchanged.1N cap.
All results preserved in `artifacts/reference/mild-packet-20261024/incident-attempt1`.
Exact producer32260/guard30676/session65212 terminal, separate idle verification
and final release complete. No body run, retiming or retry.

Selected next SOURCE/READINESS ONLY: extend physically UNCHANGED accepted REF001
threepacket input to60s. Keep ORIGINAL H5/CIC10/dt/native algorithm; only duration
and required future incident padding change. Prove exact spectral-array identity,
inherit accepted full-hull calm bounds, preserve old incident rows and require
first501 six-force samples to match frozen10s input. New incident-extension
preflight and then body each need separate source/freeze review/allocation.
Extension incident independently EXECUTED PASS:0.468s, exact8spectral arrays,
all3old wave byte prefixes, all3001new six-force samples independently reproduced,
first501 forces bit-identical to old prepared/native-expected input. Inherited
eta3.324um/u25.496um/s full-hull bounds and initial6force gates remain valid.
Report `artifacts/reference/reference-extension-20261024/incident-attempt1/diagnostic.json`,
SHAab486f2aace3b2c21820858c791f3fac3231160c74291531de2a09d4cb842340.
Exact27948/15616 terminal0, idle/separate verification/final release complete.
Original body extension EXECUTED PASS32.182s: full3001state/force/CIC/wrench/Ar
checks, first10s q/v repeatability <=2.78e-17/7.81e-18. Peaks translations
[.005964,.005844,.010048]m and rotations[.0223,.0402,.0269]deg. Same small input,
not final strong-motion case. Exact4532/24680 terminal, all142frozen match,
idle/separate readback/final release done. Bodyoutput
`artifacts/reference/reference-extension-20261024/case-attempt1`.
No encoder/automatic retry/new mesh job.

Latest user is skeptical that mesh density causes the radiation issue. Agree
root cause is UNESTABLISHED: sign change is discretization sensitivity only.
Cap read-only reassessment completed: quadrature/collocation/formulation/GF/
conditioning alternatives remain unisolated. No further refinement job or
coefficient repair is allocated. Official examples
were wave/PTO-driven, not energy-neutral; no complete wave/body/PTO energy budget
has been closed. Positive cumulative CIC work is trajectory evidence only.

Previous user requested progressively modifying the passing official WEC example;
the latest regroup above SUPERSEDES that diagnostic path.
WEC owns bounded SOURCE/PLAN/first-step readiness: native official regularCIC,
one controlled change at a time; first proposed step PTO damping1.2e6->0 with
official H5/heave constraints/ramp100/CIC20/dt.1/200s retained. Specify actual
work/mechanical-energy accounting and decision branches, review before any job.
No step was executed. STOPPED/UNRUN; no automatic suite or full6 claim from heave.

Stored-data audit changes the diagnosis: ORIGINAL60s cumulative CIC work never
negative, final+7.390385J (10s prefix+1.253276J). REJECTEDJan90held-release also
never negative cumulatively, final+237893.65J; at first5degroll34.64s already
+18042.95J. Thus negative harmonic coefficients are a warning, NOT demonstrated
net energy generation or established cause of runaway in either actual trajectory.
Researcher is completing energy/split context from stored data; no new solver job.
This is a diagnostic with an unqualified operator,
not the final buoy case. Do NOT finish or encode the old unpruned
seeded prototype as the canonical validation case.

## Preserved evidence and prototypes

REF-001 minimum17panel10s incoming-wave reference COMPLETE, independently audited:
`artifacts/reference/minimum-20261024/case-attempt1`. Final centered video:
`render-attempt3/reference-isometric.mp4` (samefolder), SHAa79a75cd88fa03c8b8d46d833cd445645e58a906995a826be9e5674ef345664f.
[Workflow](../reference/wecsim/README.md),
[baseline record](experts/wecsim_expert/minimum-reference-2026-10-24.md).
Internal consistency accepted; mesh/time accuracy still unestablished.

Neutral24ftx9ft2500kg boat unchanged:
`reference/capytaine/canonical-boat.json`, SHAea18844904e03bbfd72abc3396240c97b1994b3584e36c8149f8ae5f60a99e4f.
Native17panel BEM centroid coefficients; exactGL2 hydrostatics independently checked.
Full dry inertia/CG/frame/product mapping preserved. Four confirmed maintainedfork
corrections (import history, component sum, inertia products boundary, actual
single-total import dispatch) are historical exercised-path fixes, not radiation
bug evidence. [Fork pin](experts/wecsim-fork.json).

Three full generator prototypes: seed42/384components, Hm0~2.27-2.28m,
`artifacts/reference/generated-sea-candidates-20261024/manifest.json`. Original
whole-band/no-pruning policy remains attached to these old outputs. No scaling.
Generatorread-onlyreviewconfirms current48bands.030..495Hz alreadybelow.53346;
rawmatching-source libraryhas42,575completeobservations, notonly9catalogexamples.
Prunedsubset mustretainoriginal8*band+stratumIDs andvalues, reportenergy/period/
directionchanges andparents; no rescale/re-realize/fittedphase. Oldwhole-spectrum
loaderclosurecannotbequietlyreusedwithunchangedprovenance for a subset.
[Generator evidence](experts/topwater_expert/sea-generator-review-2026-10-24.md).
Full matrices/native pair H:
`generated-bank-native-20261024/attempt1` underartifacts/reference; official H5/
40sIRF in`generated-bank-consumer-20261024/bemio-attempt1`.17/68 scout and partiallid
controls expose strong pitch mesh sensitivity/negative damping, no qualified
replacement. NativeinfinityAPI supported; finite/infiniteGF decomposition differs.
[Producer current state](experts/capytaine_expert/current-state.md).

Old rest/held Jan90 sea case is REJECTED:107degroll/large drift. A calm0.5degroll
control with the same operator decays over60s; this does not establish passivity.
Stationary-history seed prototype then yields peaks[.368,1.549,1.520m,
5.127,1.723,2.200deg], all3001state/excitation/CIC/worldwrench checks independently
PASS. Roll>=5deg only at FINAL60.00sample, not sustainedstatistics. It is now
SUPERSEDED for final startup and incident-band policy, not falsified mathematically.
Output`generated-bank-consumer-20261024/case-jan90-seeded-attempt1`; motion
SHAd69652e3ca838afd2df6ddc249762d58321d833324038c51798399a36d055a95.
[Seeded prototype record](experts/wecsim_expert/seeded-native-case-proposal-2026-10-24.md).

Its initializer fixtures have native/source/executed independent PASS; original
Ar1e-7delay/defaultzero/extrap unchanged, nonzero startup residual reported.
Do not confuse physical radiation memory with this numerical mass-split delay.
Preserve those tools; no history adapter in the final zero-memory startup.

Geometry scalar audit validINCOMPLETE123/3001poses underoriginal.1%/1mmpolicy.
Batch1%/10mm audit was INTENTIONALLY STOPPED under new pruning steering, not
completed/numericalfailure: checkpoint1641poses,1645parseableJSONLrows, no finalfull
JSON/postcheck. Exactproducer34632 stopped after authoritative identity match;
guard37272 collected4294967295, toolsession97782exit-1. Stop record:
`artifacts/reference/actual-wave-geometry-seeded-20261024/batch-attempt1-steering-stop-record.json`.
No full60geometry/linearity acceptance. Oldearly wetturnover>=9.7087%, resolvedrim
dry; results are oldunprunedbank only. Source/fixtures/bounds remain preserved.

Preview2 independently literal/provenance PASS,10mXYtracking/fixedz/equalunits,
0/30/60PNG+all3001plot. Renderer source82083d9b, readyartifactcontrols861b1282,
rootvisualizationreceipt6cdb258e andpreviewreceipt119b3893; scope ONLYprovisional
visualization. NOencoder everallocated/run. Source model/data remain untouched;
new pruned input will need new binding/appropriate component-count handling.

## Execution and integration

[Execution ledger](experts/execution-state.json) owns liveness. All jobs terminal;
latest original60s body extension terminal/PASS, inputs verified,
idle-under-lease recorded, separate readback and final lease release completed.
No numerical owner/lease currently. SharedMATLAB9168/launcher37656 are historical
identities, not proof of liveness; authoritative reconciliation before reuse.
One serial numerical owner, frozenactualinputs/code/deps, claim?SEPARATEverify?run;
terminal?idleUNDERlease?SEPARATEreadback?FINALrelease. UTF8-sig ledger readers.

Active branch baseline/fluid-body, latest pushed checkpoint75199f6. Preserve otherworktrees,
all failed/interrupted studies and old original policies. Official examples and
true calm-arrival clarification are the current integration checkpoint.
Do not pushmain/createbookkeepingPR/merge/upstreammessages. Coordinator owns
scope/acceptance/integration; Cap nativeproducer; WEC dynamics/startup; Topwater
existinggenerator only; researcher independentmethods. Renderer agent is stopped.
