# Checkpoint handoff

Updated2026-10-24. Latest user steering controls scope. Read
[AGENTS.md](../AGENTS.md), [coordination](experts/coordination.md),
[physical contract](experts/shared-contract.md) and the relevant owner state.

## Current authorization and immediate priority

Deliver ONE credible buoy-generator-derived native WEC-Sim case with appreciable
all6DOF and roll about5-10deg. User explicitly confirms **60s is the full validation
length**. Preserve the complete startup; no cropped window or canonical suite.

LATEST steering supersedes the unpruned stationary-history prototype:

- Prune incident frequencies that the actual physics/numerical procedure cannot
  support accurately. Establish the theoretical basis first; no arbitrary cutoff
  or tuning to roll amplitude. Source remains collected-buoy data. Retained energy,
  directional spread and phases must not be rescaled/fitted. Report removed energy.
- Startup is the intended rest/zero-prior-body-motion problem: zero velocity and
  **zero radiation memory**. The prior stationary q/v/history seed is a different
  initial-value problem and is NOT the final requested startup. Specify the
  incident field and initial pose consistently; do not quietly reintroduce a
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

NEXT assignments are bounded read-only/source proposals for physics-based band
policy, direct native energy/sign verification and correct zero-memory startup.
Producer energy witness SOURCE independently PASS60c618b7, controlsc7a3605d,
proposal8689cc83. Exactly3native solves planned:17coupled2.2,17Heave1positive,
68SAMEfixedcoupled2.2. No absoluteKochinpower normalization claim, no fitted
prefactor, no cutoff/bug conclusion from3points. Checkpoint then explicit grant.
No new numerical job is allocated yet. Do NOT finish or encode the old unpruned
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
intentional stoppedbatch idleunderlease/separate7aed12/readback+finalrelease done.
No numerical owner/lease currently. SharedMATLAB9168/launcher37656 are historical
identities, not proof of liveness; authoritative reconciliation before reuse.
One serial numerical owner, frozenactualinputs/code/deps, claim?SEPARATEverify?run;
terminal?idleUNDERlease?SEPARATEreadback?FINALrelease. UTF8-sig ledger readers.

Active branch baseline/fluid-body, latest pushed4c90db7. Preserve otherworktrees,
all failed/interrupted studies and old original policies. Issue30 updated through
e5dc02a; latest pruning/zero-memory steering still needs checkpoint/issue refresh.
Do not pushmain/createbookkeepingPR/merge/upstreammessages. Coordinator owns
scope/acceptance/integration; Cap nativeproducer; WEC dynamics/startup; Topwater
existinggenerator only; researcher independentmethods. Renderer agent is stopped.
