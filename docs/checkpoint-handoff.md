# Checkpoint handoff

Updated2026-10-24. Latest user steering controls scope. Read
[AGENTS.md](../AGENTS.md), [coordination](experts/coordination.md),
[physical contract](experts/shared-contract.md) and the relevant owner state.

## Current authorization and immediate priority

LATEST USER REGROUP: stop unproven bug investigations. Deliver ONE genuinely
less energetic collected-buoy sea, locally calm water/boat rest/zero memory,
physical incoming waves with NO temporal ramp, full60s and centered VIDEO.
The video is the immediate remaining deliverable. No PTO-zero/density/energy-bug
study is allocated. Those source-ready proposals remain STOPPED/UNRUN.
Topwater solver implementation/tests/optimization remain HELD.

## Current mild-buoy result: body complete, video pending

Actual observed record46086/2023-03-19T08:40Z, source Hm0.5351635m/Tm027.42693s,
seed42/384unchanged generator carriers. Bow201.308143deg toward dominant arrival.
ExportPASS1.110s; variance/replay/frame closure verified. No amplitude/phase fit.
Parent output: `artifacts/reference/calm-buoy-20261024/export-attempt1`.

Required excitation coverage COMPLETE/PASS69.308s:816FK+diffraction solves,
same17panels, periodic headings. Original radiation H5/CIC10 unchanged.
Actual output `excitation-attempt1` under the same artifact root. Independent
input/output/FK/seam/readback checks PASS. No mesh/convergence qualification.

Physical train is an INITIAL SPATIAL window per carrier: W(s/cg+30), plateau
+/-15s/support+/-25s/quintic transition, then exact signed-kappa dispersion.
Native ramp0, no temporal force window/history seed. Parent amplitudes/phases
unchanged. Finite trains broaden the stationary spectrum; actual finite60s
statistics/support/pruned/reversed sideband measures are reported separately.
Belowomega.1, exactconstant-pressure H0=K[:,Heave] plus a reviewed quadraticomega
bridge supplies the long-wave limit. No new native radiation solve/repair.
IMMUTABLE accepted physical data: `train-attempt2/train.npz` and `body-force.mat`.
The original train-attempt2 report remains REJECTED under its old0.1N cap.

NEW explicit engineering calm policy: eps=.001, ||F0||<=epsmg=24.525N,
||M0||<=epsmgL=179.405Nm (L7.3152m), plus each component<=epsSi in SAME native
basis. Existing eta/u gates unchanged. This is not universal accuracy proof;
old0.1N rejections remain preserved. Independent NEW assessment PASS20.984s,
`assessment-attempt1/diagnostic.json`, SHA367a08f33af6fec4aca2cf410d2d8492a294dcaa7ba2e02d6ccefd595113f1e2.
Fullhull24eta/72velocity boxes completely covered: eta<=21.034um,
speed<=2.290736mm/s, force norm12.7286N(~.052%weight), moment1.1623Nm.

Native60s BODY COMPLETE/PASS27.776s, independently audited all3001samples:
`artifacts/reference/calm-buoy-20261024/case-attempt1`.
Original H5/CIC10/dt.02/all6free, q/v/memory0, ramp0/viscous0, no prerun/crop.
Applied/prepared6forces and fullwrench sum EXACT0 error; CIC/rotatedSolidNE/
nativeAr/fullinertia/frame/finite checks PASS. Peaks native surge/sway/heave
[.238006,.274488,.240298]m; roll/pitch/yaw[1.481170,.662514,.616524]deg.
Native modelworld is meanboat-aligned; station-global R201.308/t=[0,0,0] is
embedding metadata applied ONCE to both field and body for display/export.
No M/K/H/force double rotation or extra CG phase.

Motion SHA f822fde31bae1f5fd4e0691d6fd9467143786af6631cc981c591df63cd8daaf4;
raw native6b8a9adc5224d6bbb92541af3c6350771b99ee638af98e7c30d4999d31739e57;
summaryaab25bc9a17af2900269555debdc0732baf6a76effe87090af4b835895b539bf.
Exact29068/39568 terminal0/1fa9ae, all146frozen identities match, idle/separate
readback/final release dd4b47. No numerical owner unless ledger says otherwise.

ROOT genuinely accepted visualization after independent body review:
`case-attempt1/render-motion-acceptance.json`, SHA7c52776c9342d3ed69b9089d40d543f357c623586d446f4e396dd846852a6767,
status provisional_model_trajectory. No mesh/time/nonlinear accuracy proof.

## Immediate next: video

Video owner `/root/calm_video` owns only new renderer/control/binding helpers;
WECbody owner `/root/wecsim_reference` stopped after body evidence.
`render_calm_buoy_reference.py` SOURCE independent PASS5076954c; exactsame
train2/profileAPI,12m close/equalmetres/fixedz/meanheading camera,751nativeframes
at12.5fps including0/60. Native-to-global mapping and CG applied once.
Pure binding COMPLETE: artifact controls0edd1468, preview proposal61d7eaf49dc6f8dfae915fb74177ab1f5704033e07a860a3187f7eea857ac508
in `render-preview-attempt1-source/proposal.json` under artifact root.
Actual binding/freeze and PREVIEW independently PASS28.609s. All751surfaces
cached,0/30/60PNGs and completeall6motionplot literally reviewed byROOT and
independent researcher. All751hullposes reproduce rawdata exactly, no clipping;
surfacecache-direct error2.203um below56.757um displaybound. Stable preview
`render-preview-attempt1`, preflightbaf5301d/cachef6af6151. Exact27940/29344
terminal/released. GenuineROOT preview receipt `render-preview-review.json`,
SHA86e43d7047f32f0c04a3cbd2b7577d066f4a0136c74c59a9659cee92c1381c09.
NEXT pure encoder binding/freeze review, then SEPARATE encoder allocation(max900)
reusing exact accepted cache. No encoder allocated yet. Never finish old seeded
prototype video. Body/packet/assessment/renderer/cache must stay unchanged.

WEC-Sim uses open-water frequency-domain hydrodynamics/body-only integration,
not a finite ocean patch. The physical incident field extends beyond camera;
render is prescribed incident surface, not computed scattered/radiated water.

Key evidence: [mild input](experts/topwater_expert/calm-buoy-export-readiness-2026-10-24.md),
[excitation coverage](experts/capytaine_expert/calm-excitation-coverage-2026-10-24.md),
[arrival/body](experts/wecsim_expert/calm-buoy-arrival-readiness-2026-10-24.md),
[video](experts/wecsim_expert/calm-video-readiness-2026-10-24.md).

## Interpretation corrections preserved

No WEC-Sim/Capytaine radiation bug is established.17-to68 sign change showed
spatial-discretization sensitivity, not panel density as identified cause.
Both official RM3 examples pass; providedCIC golden<=4.157e-13m vs1e-5m tolerance.
Original unchanged incoming mild case also passed60s with its10s prefix identical.
Actual cumulative CIC work never negative in either original60s(+7.390J) or old
Jan90established-sea case(+237.894kJ). The latter's107degroll/11.5m sway made it
unsuitable for linear validation, not a demonstrated numerical failure.
[Actual-work review](experts/numerical_methods_researcher/actual-radiation-work-review-2026-10-24.md).
No complete wave/body/PTO/fluid-storage energy balance was closed. Radiation loss
is expected; instantaneous stored-energy return is not itself energy creation.

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
