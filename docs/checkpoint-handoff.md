# Checkpoint handoff

Updated2026-10-25. Latest user steering controls scope. Read
[AGENTS.md](../AGENTS.md), [coordination](experts/coordination.md),
[physical contract](experts/shared-contract.md) and the relevant owner state.

## Current authorization and immediate priority

IMPLEMENTATION AUTHORIZED: latest user request is to implement the canonical
example/workflow plan developed in this conversation. One flat
`examples/simple-boat/` contains `boat.npz`, `sea.npz`, `result.npz`, `video.mp4`.
Public generate-sea/solve/render commands take an optional directory (default
cwd), fixed filenames; solve additionally requires duration and defaults native.
Native dynamics implementation remains excluded; native must fail clearly.
WEC-Sim is the first working solve backend and must export a reconstructed total
free surface, not an incident-only surrogate. Reconstruction reads computed
motion without another body solve or feedback. Measure total-minus-incident
surface RMS/resolved peaks (full patch and near hull), without tuning to smallness.

Boat input is an uncapped convex physical triangular hull with exactly one
inferred opening. Result includes complete exported hull geometry (possibly
remeshed), CG, absolute poses and independently sampled variable-topology total
water geometry. Numeric NPZ only, no general metadata/solver marker/pickle.
Shared code handles equilibrium/rest/zero-memory/upstream spatial train; requested
duration does not alter startup. Fatal errors preserve old outputs; physical
applicability departures warn (10% restoring/wetted-set engineering policy).
Every public command needs its own Markdown page including hardcoded choices,
fatal/warning detection and limitations. Final binaries use Git LFS.

Start with bounded mesh/FOSS probes, then contracts/CLI/docs, WEC reconstruction,
canonical run/video. Use expert ownership and independent adversarial review;
commit/push reviewed checkpoints only to baseline/fluid-body. Existing dirty
changes below are the previously authorized completed renderer cleanup and must
be preserved/checkpointed. Numerical work is serial under the central lease.
On entry no lease existed; authoritative CIM inspection confirmed prior pytest
PIDs40692/16928 absent, no Python/FFmpeg/Gmsh jobs; shared MATLAB32780/9168
remains present and must be queried before reuse. No new numerical launch yet.

## Previous completed visualization checkpoint

FULL60s requested2D video and parallel read-only speed review are COMPLETE.
USER approved2D and deletion of3D rendering work. Active3D render backends,
render-only reference helpers/configs, and PyVista/VTK project dependencies are
removed. Physical3D meshes, intersection math, solvers, native WEC-Sim/fork code,
old evidence and MP4s are preserved. No further job is allocated; Topwater
solver/physics implementation remains HELD. Await next user steering.

Watch `artifacts/reference/stronger-buoy-20261024/matplotlib-full60-attempt1.mp4`.
SHA25615510982f6cef3a246f1129222e0525052215a06cfa2e0073af3c53bf803f23b.
H264800x800/751frames/12.5fps/60.08s, including native60s endpoint framehold.
All751 decoded successfully. Owner viewed0/30/51.44/60stills; ROOT viewed0+51.44.
QA is in `matplotlib-full60-qa-attempt1`; peak frame is frame-0643.png.

Final style: top-down fixed view, continuous coolwarm surface elevation about0,
opaque offwhite hull, gray mesh, black rim, darkblue actualwaterline/whitehalo.
Monospace signed fixed-width14.5pt header spans nearly800px and shows TIME/ROLL/
PITCH/HEAVE. Same original accepted751 hull/time and central65x65 XYZ data are
BITEXACT; native q rows0:4:3000 provide exact HUD. No scaling/heading/physics/IC
change. Calm startup remains as previously accepted.

Full input uses75x77 waternodes: original12m cache plus ONLY1550 true-field border
nodes per frame. Waterhalves6.9375/7.125m cover fixedview half6.100702m atalltimes.
Prepared scene~104MB, SHA884667cb0a43f55ec86303c38fae303a76c58a8869c547b1d3290f70b4f330bc.
No96m fullfield regeneration. Preparation40.657s; publicCLI148.792s, renderer
145.75s, wholeproducer+QA197.911s. Geometry/update64.066s, Aggdraw68.732s,
pipewait.530s, encoderstartup1.937s/finalflush1.703s, hullvalidation5.359s,
setup3.391s, input.157s. Pipewait is not total encoder CPU; work overlaps.

Public command remains `uv run topwater render data/scene.npz videos/preview.mp4`.
Both required positional paths support relative and absolute paths. Wholeinput,
uniform timestamp cadence, fixed800square, namedMP4 only; no --extra/input/output/
start/end/size flags. Optional JSON settings now only framing/timeout. Core
Matplotlib+bundledencoder, realAggblit/directbufferstream/reusedartists. Prepared
SceneSequence uses worldmetres/right-handed-z-up/fixedtriangularconnectivity,
named telemetry arrays and explicit convexmesh validation. No boat constants.

Parallel `/root/startup_review` read code/primary sources only: no edits, tests,
benchmarks or concurrent numericaljobs. Promising later exact changes are
coveragefinder reuse (preserve holes/deformation fallback) and tighter contact
broadphase. Faster encoder preset needs quality review atfixedbitrate; regular
grid image fastpath changes interpolation and needs visualapproval. Not applied
during fullrender. [Review](experts/visualization_expert/render-speed-review-2026-10-25.md).

Exact fulljob originalsession42135/producer40788/guard28968 terminal6f2f47/CIM6c6e36;
frozen identities PASS, idle/separate58a912/FINALreleasebacd1e. Ledger idle/no lease.
Owner `/root/visualization` (topwater_expert) wrote/executed initial sharedlibrary;
`/root/wecsim_reference` owns producer adapters. Dedicated visualization_expert
definition/KB/routing exists and inherits parent settings, but native spawn in
this client returned unknownagent_type; configured NOTloaded here. Recheckafter
clientrefresh and transferatstablecheckpoint; no separateCLIcoordinator.

Previous checkpoints8450c2c/662a148/b17f529 pushed baseline/fluid-body. Current
full-render/cleanup checkpoint follows. See [library evidence](experts/topwater_expert/visualization-workflow-2026-10-24.md)
and [WEC producer state](experts/wecsim_expert/current-state.md).


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
