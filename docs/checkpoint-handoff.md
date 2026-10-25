# Checkpoint handoff

Updated2026-10-24. Latest user steering controls scope. Read
[AGENTS.md](../AGENTS.md), [coordination](experts/coordination.md),
[physical contract](experts/shared-contract.md) and the relevant owner state.

## Current authorization and immediate priority

LATEST USER: render FULL60s withaccepted2Dstyle, whileparallel READONLY search
findslowhangingperformancewins. USERexplicitlyapprovesdeletingALL3Drendering
code/deps ('2D is good'). No deletionofphysical3Dgeometry/solvers oroldartifacts.
Owners holdcodebriefly forROOTcheckpointcurrentrefinement, then clean3Dpaths
BEFOREfreezingnewfullproducer. Genericowner `/root/visualization` owns2Dlibrary/
PyVista-backend/dependencycleanup andserialfullrender; WECowner ownsobsolete
3Drender-onlyhelpers cleanup/shareddatautilities+newFULL751preparedinput.
UseORIGINALfull12m751cache+exactnativeHUD, MINIMALsame-fieldborderonlyifneeded
forfixedcamera/alltimecoverage; central65x65BITEXACT, no96mfullfield/wavephysics.
Proposalfirst, max~600s/1CPU/noMATLAB, full751decode+peakPNG+phasewalltimings.
No fulljobgrantedyet. `/root/startup_review` ownsparallelREADONLYcode/docsperf
analysis, no numericaljobs/edits/benchmarks. DoNOTmutatefrozenproducerduringjob.

Requested minimal top-down Matplotlib preview and standard library workflow are
COMPLETE atcheckpoint662a148 (pushed/issue30updated). LATEST USER authorizes
refinement: monospacefixed-widthHUD, opaqueoffwhitehullinterior, mutedwaterpalette,
darkbluewaterline andFULLCLI/render/encoderphase timing withspeedrecommendations.
ROOT GRANTED onecached90sjob to `/root/visualization`, exactmanifest
`matplotlib-clip-attempt2-source/proposal.json` SHA49de450cdc57026323345a8f76e6b773c5e18e9a3406fd75d4d36d5cd94de421.
Twofocusedheader/footprinttests, same2475telemetryNPZ, cleanPUBLICCLI63, peakPNG/
full63decode. No augmentation/waves/deps/physics. Freeze/claim/separateverify/
exactterminal/idleunderlease/separateverify/finalrelease asusual. Sendstilland
timingspromptly; no extraencoderbenchmarks. No full60s encode allocated.
Refinementattempt2 COMPLETE: externalCLI12.205s/render10.625s, geometry3.326,
Aggdraw4.283, pipewait.110, encoderflush2.187. ROOTacceptedopaquehull/monofixed
header/mutedwaterstyle. Tests2/full63decodePASS; terminal/frozen/release2abc98.
LATEST USER askslargerheaderacrossmost800pxwidth. ROOT GRANTED FONT-ONLYattempt3,
manifest`matplotlib-clip-attempt3-source/proposal.json`
SHA9a77d47fe82c07b54763bd51b6958854c7187eb2290a5ca80075b1349f1b4e87,
sameowner/same2475input/style/timers,14.5pt(~766pxheader). Onecached90sCLI63+
peakPNG/full63decode; no tests/deps/waves/augmentation. Serialleaseprocedure.
Topwater solver/physics implementation remains HELD. USER hasnowapproved2Dand
3Drendercode/dependencyremoval; oldsource/evidence remainrecoverableinhistory.

Watch `artifacts/reference/stronger-buoy-20261024/matplotlib-clip-attempt1.mp4`.
SHA256 954343c52b6bc5690fcf421bc0ac6804c2a4643b2571c26ffc34682f44e87eed.
Still: `matplotlib-qa-attempt1/poster.png` under the same artifact root.
Source50.00-54.96s, H264800x800/63frames/12.5fps/5.04s. ROOT viewed peak51.44s
poster and sent still/video links to USER. Full63-frame decode and17 focused
tests passed. Same cached hull/water/time arrays BITEXACT; accepted native q
supplies roll_deg/pitch_deg/heave_m, no inference or re-zeroing.

Style: fixed top-down map, continuous symmetric blue-white-red surface elevation
about mean0 (not seabed depth), no water grid, light-gray hull mesh lines, black
mesh-boundary upper rim, actual blue hull-water intersection with white halo.
Video header shows TIME/ROLL/PITCH/HEAVE. Convex meshes are an explicit validated
requirement; supplied mesh/data determine shape and resolution.

Public command: `uv run topwater render data/scene.npz videos/preview.mp4`.
BOTH arguments are required positional paths; examples prefer relative paths,
absolute paths also work. Always render the whole NPZ at its uniform timestamps.
Fixed800x800; no input/output, start/end or size flags. Optional JSON settings
control presentation. Matplotlib and encoder are CORE deps, so no --extra for2D.
Only named MP4 persists from normal rendering; scratch is cleaned and timing/
diagnostics return on stdout. Temporary reference QA is allowed while tuning.

SceneSequence contains named arrays time_s, hull_points[T,N,3], hull_faces[F,3],
water_points[T,M,3], water_faces[K,3], metadata, optional neutral telemetry arrays
serialized as telemetry__roll_deg/pitch_deg/heave_m. Header requires supplied
telemetry. Metres/right-handed-z-up; fixed triangular connectivity per sequence,
uniform timestamps. Producers own adapters; renderer has no WEC-Sim/boat constants.
Augmented input is `scene-telemetry-attempt1.npz` under the artifact root.

One Agg figure/reused artists, real background blitting, direct canvas-buffer
streaming to bundled FFmpeg. Culling draws8064 of107648 water faces. Actual
render11.25s (setup.765/validation.094/frameavg.120s), wholetests+augmentation+
CLI+QA17.797s. No wave preparation or GPU; triangle/triangle geometry computes
the actual waterline separately from standard Matplotlib drawing. No custom
rasterizer. Original prepared field required100.75s once; styles reuse it.

Initial implementation owner `/root/visualization` is topwater_expert; WEC owner
`/root/wecsim_reference` supplies telemetry/data adapters. Dedicated
visualization_expert role is defined with inherited settings, but this client's
native spawn returned unknownagent_type. It is configured, NOT loaded here;
recheck after client refresh and transfer at a stable checkpoint. No separate
CLI coordinator was started. Execution ledger owns exact terminal/frozen/idle/
release evidence: session81549, release5ec2c5. All numerical processes terminal.

Previous shared-PyVista checkpoint8450c2c was pushed and issue30 updated. Its
accepted oblique MP4 and failed/rejected attempts remain preserved. No further
3D development is authorized while awaiting 2D review; no deletion yet.
See [library evidence](experts/topwater_expert/visualization-workflow-2026-10-24.md),
[adapter evidence](experts/wecsim_expert/telemetry-2d-readiness-2026-10-24.md),
[visualization role](experts/visualization_expert/current-state.md).


PREVIOUS USER REQUEST COMPLETE: render the new video and provide where to watch it, using
the agreed readability changes. Remove axes and surface
grid, dark-to-light blue filled water, lower camera, shaded flat hull faces,
clear inside/outside and waterline, slight hull transparency if readable.
Same frozen stronger-case motion/packet/cache and actual metres/time; no new
sea/solver/physics/geometry correction or amplitude exaggeration. WEC owner
`/root/wecsim_reference` owns new presentation renderer/controls/source-readiness.
NEW VIDEO COMPLETE. Watch local
`artifacts/reference/stronger-buoy-20261024/presentation-encode-attempt2/reference-presentation.mp4`.
SHA2567a9debb31a440b6967fd3340f2fcd6f89072699dd7ecc2e36136d3bca6079aa8.
H2641280x900/12.5fps/751frames/60.08s (full60s simulation plus endpoint frame).
Selected24deg lower camera/10m equal-metre span, no axes/grid, fixed dark-to-light
blue height palette, flat shaded ivory exterior/warm interior, actual cached
waterline and slight underwater visibility. Camera followsXY only, fixedheading
and verticalreference. Physical inputs and all751 cached poses/waves unchanged.
Original movies and failed/slower rendering previews are preserved.

Final render1064s/17m44s. Full751decode PASS; ROOT viewed encoded0/30/51.44/60s.
Exact cached-water ray residual<=4.44e-16m/contact<=2.22e-16m. Distantwater uses
320x225 display pass; hull/HUD/contact and water queries near hull retain1280x900.
Finite waterpatch edges remain visible; shading/contact are presentation, not
computed flooding/wetting or a new physical validation. No sea-state, solver or
Topwater physics work is allocated. Execution ledger owns terminal
and release evidence for originalsession14472/manifestcbf276209e7af879a17177b89d0296f9b803cffb751ff8a3601f67f736498467.

LATEST user requested MORE energetic actual buoy sea after the mild clip,
with radiation-energy diagnostics and a video. STRONGER CASE AND VIDEO NOW
COMPLETE. No further case/bug/mesh/Topwater work is allocated. Preserve the
provisional physics/geometry limits below; do not turn amplitude achievement or
positive radiation work into an accuracy/linearity proof.

Selected actual record46086/2023-12-23T00:40Z,row29984,
Hm01.657467948/Tm027.236567149. Neutral exportPASS1.093s,
variance.1717m2/384seed42carriers/defaultbow182.067194deg. Actual measured a/MEM
change; frequencies/phases/IDs arrayexact to mild seed42, no uniform scaling.
Artifact root `artifacts/reference/stronger-buoy-20261024`.

Physical spatial tc40/plateau15/support25/quintic/signed-kappa train completed
2.891s; same nativefixed17fullheading H/H0bridge/originalH5/CIC10/dt.02.
Calm assessmentPASS28.266s, full24eta/72velocity boxes: eta<=.83572mm,
u<=3.686mm/s, normF6.104N<24.525, normM12.016Nm<179.405 and relative axes/arrival.
NativebodyPASS31.349s, all3001raw/state/IC/Fe/CIC/NE/Ar/mass/frame checks
independently PASS; ALLprepared/applied6forces exact0error. No ramp/viscous/PTO/
mooring/history seed/crop; q/v/memory0. Modelworld nativealigned; station-global
R182.067194/t0 metadata applied once, no M/K/H double rotation.
Peaks xyz[.562349,.836958,.758899]m and rpy[6.618565,2.677850,1.325248]deg.
Workaudit independently reproduced: netCIC32971.03955J/mincumulative0,
temporaryreturn4074.89555J, incident46582.19287J. Expected netbody-to-radiation
transfer for this trajectory, not a complete fluid-energy balance/passivity proof.

Strongermotion3eae21f9efa9987d77225b0ebe3f2cac4babec0a482ea01bbbf7801370e6b072;
raw d350bac332d13ad94d3684395beec16ec95cc67c0d78f8423fcc6f6d7782f6da;
summary709f3c034d14bb0b6fee0f2609a9fd4918d5be1cc245c92617839af02f29e0c4.
ROOTmotion visualization receipt a9169af8 and preview receipt8386d482 genuine.
PreviewPASS35.641s, full751pose/cache/frame/clip checks andROOT+independent
literal0/30/60+all6plotPASS. EncoderPASS220.375s, exactsession41948, terminal a7248d/exit0,
FFmpeg35700exit0; guard38024/producer30412/renderer11444/FFmpeg all absent
CIM8aa409. All18inputs/deps unchanged; idle d12912->separatee948bf->release283833.
No remaining numerical job/lease unless ledger says otherwise.

FINALSTRONGER VIDEO `render-encode-attempt1/reference-isometric.mp4`,
SHA42c8ca1c80e961784200d076faff7e9a34a67c782120019e606bbf553ec004e2.
H2641280x900/12.5fps/751frames/60.08s (full native60 endpoint hold).
ROOTfull751-frame decode/metadata/frozenhash check PASS; actualdecoded0/30/60
literally PASSROOT+independent researcher. `render-encode-attempt1/final-review/verification.json`.
Video linked to user; 12m equal-metre close view, sameexactprescribed incident
field, no motion/wave amplification. Openwater physics is not limited to camera.

LINEAR-USE ASSUMPTION IS APPROXIMATE: cached751-time local-waterplane proxy gives
relativeheel3.10deg/trim2.88deg/relativeheave.090m, not6.62deg relativeflatheel.
At near-globalrollpeak relativeheel2.70deg. Existing3deg static curves show7.3%
righting-secant/15.3%tangent departures; combined curvature/wetting/rim/flooding
and nonlinear/dynamic accuracy remain unqualified. See stronger-relative-water-
proxy note. User asked an a-priori geometry-based assumption, not proof from an
unavailable nonlinear model. Do NOT claim strict negligible geometry change.
Minimum-mesh end-to-end workflow and requested motion/video are complete;
mesh/time/stricter physical applicability qualification is a distinct future phase.

Key stronger records: stronger-buoy-export-readiness (Topwater), stronger-buoy-
readiness and stronger-video-readiness (WEC), stronger-relative-water-proxy.
Commands/source identities in `reference/wecsim/stronger-buoy-pipeline.json`.
Minimum-reference completion audit: [requirements and evidence](experts/minimum-reference-completion-2026-10-24.md).
Final repository checkpoint/issue30 refresh follows the completed video audit.

Earlier regroup remains standing:
Deliver ONE genuinely
less energetic collected-buoy sea, locally calm water/boat rest/zero memory,
physical incoming waves with NO temporal ramp, full60s and centered VIDEO.
No PTO-zero/density/energy-bug
study is allocated. Those source-ready proposals remain STOPPED/UNRUN.
Topwater solver implementation/tests/optimization remain HELD.

## Preserved mild-buoy result: body AND video complete

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

## Preserved mild video evidence

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
Encoder COMPLETE206.563s from exact accepted cache. Output
`render-encode-attempt1/reference-isometric.mp4`, SHA205400a829eef3ffb9b1fdf83a359299dc05e4edf15f1f051f3765fb5ca7e8f9.
H2641280x900/12.5fps/751frames/60.08s; native60s endpoint held one frame.
Owner hit usage limit after launch; ROOT recovered without rerun. Authoritative
CIMab8428 exact17992/3976+children absent; renderer complete report follows
FFmpegreturncode0 assertion. Original tool/producer exit status unavailable,
recorded NULL rather than invented0. All18/dependency/output hashes verified,
full751frames decoded, actualdecoded0/30/60PNG literally PASSROOT+researcher.
`render-encode-attempt1/final-review/verification.json` and recovery.json record
evidence. ROOTidle-underleasee53374->separate3c7174->finalrelease1663ac.
Mild video delivered; body/packet/assessment/renderer/cache unchanged.

Mild radiation-work audit independently reproduced from all3001saved samples:
netCIC+3271.532294J,mincumulative0, temporarynegative-powerreturn309.915673J,
incidentwork3492.513729J. No net radiation creation observed; not a complete
fluid-energy balance. `case-attempt1/radiation-work-audit.json`.

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
latest stronger encoder terminal/PASS, inputs verified,
idle-under-lease recorded, separate readback and final lease release completed.
No numerical owner/lease currently. SharedMATLAB9168/launcher37656 are historical
identities, not proof of liveness; authoritative reconciliation before reuse.
One serial numerical owner, frozenactualinputs/code/deps, claim -> SEPARATEverify -> run;
terminal -> idleUNDERlease -> SEPARATEreadback -> FINALrelease. UTF8-sig ledger readers.

Active branch baseline/fluid-body; latest published prior checkpoint69d2431.
Final stronger-case/video checkpoint is recorded in git log and issue30.
Preserve otherworktrees, all failed/interrupted studies and old original policies.
Do not pushmain/createbookkeepingPR/merge/upstreammessages. Coordinator owns
scope/acceptance/integration; Cap nativeproducer; WEC dynamics/startup; Topwater
existinggenerator only; researcher independentmethods. Renderer agent is stopped.
