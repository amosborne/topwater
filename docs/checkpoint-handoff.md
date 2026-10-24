# Checkpoint handoff

Updated 2026-10-24. Latest user instructions control scope; stable routing is in
[AGENTS.md](../AGENTS.md). Read [assignments](experts/assignments.md) and the
[physical contract](experts/shared-contract.md) before assigning numerical work.

## Current authorization and priority

REF-001 minimum-mesh reference is COMPLETE and independently audited.
Latest user steering authorizes REF-002: ONE credible, fresh collected-buoy
sea-state-generator case with roll about5-10deg and appreciable other5DOF.
No scaling the old sea, no bespoke fitted spectrum, no canonical suite.
Default box direction faces dominant incoming energy; export all direction/frame
conventions and chosen initialheading. Internal state/heading selection is allowed.

Assume a linear operating range A PRIORI from wetted/buoyancy/righting curves
against WEC-Sim's mean assumptions. User explicitly does not request a nonlinear
trajectory model; prior weakNL/adapter proposals are superseded. Topwater SOLVER
implementation/tests/simulation/optimization remain HELD. Existing-generator
reference-input export and external-reference work are authorized.

## Completed trusted baseline

- Same canonical boat input SHA256ea18844904e03bbfd72abc3396240c97b1994b3584e36c8149f8ae5f60a99e4f.
  Minimum17-panel Capytaine export1086solves/69.9800752s; exactGL2 hydrostatics
  independently checked with zero extra BEM solves, original coefficients retained.
  Local native data: artifacts/reference/minimum-20261024/native-gl2.
- Official BEMIO H5 and executed metadata/phase/matrix signoff pass. Maintained
  fork pin f962103a348d63101ce669c2611f7d7c109755c3, published user-owned branch.
  Four reviewed corrections: history padding, component sum, inertia boundary
  order, actual rigid import dispatch once. No upstream PR/main push.
- Accepted packet-attempt3 + packet-full-planform preserve fullhull continuous
  eta3.324um/velocity25.496um/s bounds, six initialforce gates, band and59.4%
  arrival. No ramp/clipping. [Independent input review](experts/numerical_methods_researcher/spectral-packet-proposal-2026-10-24.md).
- Native case-attempt1:10s/501samples/dt0.02, actualrest, zero velocity/radiation
  initialmemory, runtimeM+A error1.59e-12 and applied501x6excitation error1.61e-12.
  Peak translations5.222/4.970/5.935mm, roll0.0191deg. All6motions finite.
  [Independent output review](experts/capytaine_expert/native-motion-review-2026-10-24.md)
  explains stock1e-7TransportDelay initialzero/extrapolation and verifies all501
  Solid Newton-Euler balances. No correction warranted; temporal/mesh accuracy
  remains UNESTABLISHED, not a continuum/Topwater qualification.
- FINAL centered video: artifacts/reference/minimum-20261024/case-attempt1/render-attempt3/reference-isometric.mp4.
  Literal encoded0/5/10s review PASS; H2641280x900,251frames25fps,10.04s including
  lasttimestamp10 endpoint hold. Equal metres/orthographiciso/CGonce/native
  RzRyRx, transparent prescribed incident surface and mm/degree telemetry.
  SHA256a79a75cd88fa03c8b8d46d833cd445645e58a906995a826be9e5674ef345664f.
  Six-motions.png samefolder. Draft1/2 unchanged/preserved; user has final link.
  [Workflow/artifacts](../reference/wecsim/README.md),
  [run record](experts/wecsim_expert/minimum-reference-2026-10-24.md).

## Active strong-motion case preparation

Existing topwater.sea_states already generates higher-energy whole observed
spectra; no generator/source change or Hs multiplier is needed. Three INTERNAL
candidate exports are complete, seed42/384components each:
artifacts/reference/generated-sea-candidates-20261024/manifest.json.
Hm0=2.27842/2.27614/2.27033m with distinct periods; wholebandvariance, exactreplay,
physicalsupport, source/binary hashes and global/bodyphase closure pass.
Default incoming bow headings165.963/155.230/163.218deg. Keep observed/MEM
spreading/phases/globalbank immutable; these are not yet selected canonical sea.
[Generator review](experts/topwater_expert/sea-generator-review-2026-10-24.md).

Completed exact-geometry study: geometry-envelope-20261024/attempt2 under the
artifact root,49freeheel roots/175staticposes/951evaluations, all independent
checks PASS. Relative flatwater heel5deg: righting12.33%below meanlinear,
tangent departure24.03%, wetsetturnover15.87%; at10deg22.33/39.32/28.82%.
Stable/no rimflood through12deg freebranch. Exploratory5%secant/10%tangent labels
give only±1.5deg relativeheel; not an adopted globalroll or dynamicaccuracy gate.
Interpret relativewaterlevel/slopes and fullfootprint curvature, not globalheave
alone. Strong globalwavefollowingroll is a hypothesis needing candidate evidence.

Completed native bank producer: generated-bank-native-20261024/attempt1,
2838freshsolves/67.119s, independentFK/sign/schema/phase/framechecks PASS. Exact
384x6H for eachdefaultalignedbank, full240-point .05:.05:12 radiation+directAinf
and actualheading0 compatibility data, same17panel centroids/exactGL2metadata.
No oldheading interpolation or omitted directional energy.

Consumer officialBEMIO40sIRF/H5 readback and harmonic sign/frame checks PASS.
Stationary PARTICULAR-only roll peaks: April7.5509deg, March10.2007deg,
January3.8098deg. April otherpeaks: surge1.311m/sway.725m/heave1.582m,
pitch3.913deg/yaw2.500deg. These do NOT qualify actualrelease/startup/drift or
geometry. Completed stored-data native-vs40sK comparison (compare-attempt1,
terminal0/1.859s) gives April roll-history difference0.3957deg/heave0.05393m;
relative-roll fixed-mean proxies about6deg do not demonstrate wave-following
rescue.40s quadrature/window effects are small compared with the native-band
inconsistency. Signed heave gap has dominant NEGATIVE -3573kg; positive known
B33 and a passive omitted upper tail alone cannot explain that sign. No blanket
high-band extension, repaired coefficients or trusted strong-motion run follows.

Completed one independently reviewed450-call diagnostic:66 radiation calls on
unchanged17/68-panel surfaces plus384 exact January offset90 pairs. Output:
artifacts/reference/radiation-heading-scout-20261024/attempt1, terminal0/56.9945s,
all source/surface/reuse/phase/frame integrity checks PASS; lease released after
separate idle verification. Hydrodynamic adequacy FAIL:68-panel B33(12)=-48118kg/s
and B55(12)=-90212kgm2/s imply negative radiation power. Incident-band pitch A55
changes about10944->24870 at1rad/s and11439->25958 at2rad/s; heave changes smaller.
Sparse17-panel high-frequency A33 stays about2500 rather than directAinf10541.
No coefficient repair or trusted canonical case follows. No272 cascade.

Completed Jan90 stored-data screen, terminal0/1.828s, independent integrity
review PASS, lease released. Native/CIC40 particular globalroll5.1188/5.2209deg;
othernativepeaks surge.339m/sway1.556m/heave1.512m/pitch1.729deg/yaw2.311deg.
Exact rotation relative to FIXEDmeanCG plane gives roll2.105/2.133deg and
pitch2.484/2.508deg; relativeheaveproxy.018/.030m. Sampled footprint departures
initial.183m/peak.307m remain material against.402m draft. Retain ONEJanuary
heading90 candidate for a provisional actualnativeWEC run, not canonicalaccuracy
or geometric acceptance. This explicit beamheading is an internal selected-case
hypothesis, not a rendered-box rotation or change to the global bank.

Actual native Jan90 held-release60s case is COMPLETE but REJECTED for the
requested operating envelope: roll107.53deg/sway11.474m/yaw27.88deg. Full native
3001samples and60columns preserved at generated-bank-consumer-20261024/
case-jan90-attempt1. Sourceb58b5bb6; exactq0/v0/radiation0, H5/M/frame checks PASS,
all3001x6 actually applied excitation errors EXACT0. Preprocessing29.624s,
compile2.978s,dynamics4.352s; guarded total43.762s terminal0, exactclient16928/
guard31416 absent and lease released after separate idle readback. NOT goal
achieved, not a canonicalcase; no croppedwindow/rerun/encoder authorized.

Stored native review: all3001 actualSolid NewtonEuler balances close to4e-10;
force sums exact and documented residualacceleration extrapolation closes to9e-9.
Roll firstcrosses5deg34.64s/10deg34.86s/20deg36.54s; free sway/yaw drift precedes
runaway. Harmonic particularv0 sway-.1938m/s,yaw+.643deg/s differs from heldrelease
v0=0, consistent with large free modes. No posthoc croppedcanonicalwindow.
Full36radiation has confirmed negative harmonicpowerdirections in incidentband
(native scaledminimum-184.99; CIC40-196.28); independent checks PASS. This doesnot
attribute capsize. Effective mass stayspositive; simple frozenfeedback gains<1
throughonset do not explain fullcoupled/timevaryinggrowth. No coefficientrepair.

Calmwater0.5degroll control COMPLETE30.555s/terminal0, independently reviewed
PASS: all3001x6excitationexact0, nativeworldwrenchbalance6.83e-13/Ardelay1.14e-12.
Maxroll initial.5deg; final.0000426deg; successive10swindows decay strongly.
Swaydrift only.01576m/yaw.00347deg. Exactclient36696/guard30692 terminal and
lease released after separateidleproof. This small initialdirection doesnot show
severegrowth; doesnot prove fulloperatorstability/passivity/forcedcaseadequacy.

WEC concrete IC proposal is ready (SOURCE/PHYSICALPROPOSAL ONLY):
experts/wecsim_expert/stationary-history-ic-proposal-2026-10-24.md. SameJan90 bank,
finite40sK particular q/v as deterministic seed, exactworldomega and CG/pivot map;
missing CIC velocity prehistory must be added to all THREE native force-history
rows BEFORE native PCHIP extrapolation (PCHIP is not additive). No prescribed
motion after0/customODE/damping/phasefit. Meanlinear initializer is NOT exact
stationary native nonlinearrigidbody/discrete/delay solution; report residuals.

Seededinitializer fixture attempt3 COMPLETE22.951s/terminal0, independent
EXECUTED PASS. Actual17fastCIC/PCHIP error1.89e-14; independentSciPy7.44e-15;
forcehistoryrows2.55e-15/nativejointq-v-worldomega roundoff. Exact39780/1420
terminal/idle-underlease-separateverify-release complete. Attempts1/2 failed
beforephysics (danglingkernelinputline, theninvalidFromWorkspaceafter-final mode);
alllogs/partialmodels/originalsources preserved. Fixturefinalsource5385c4c7.

Seeded60sJan90 nativecase COMPLETE38.318s/terminal0, independent EXECUTED
implementation PASS. Output generated-bank-consumer-20261024/case-jan90-seeded-attempt1.
Actual3001q0/nonzeroCGv-worldomega/Fe allPASS(excitationerrorEXACT0). Independent
pre0+ACTUALpost0 convolution/threeforcehistoryrows/nativePCHIP<=1.10e-11;Fradbus0.
ActualrotatedSolid NewtonEuler<=1.46e-11N/5.92e-12Nm; unchangedArrecurrence2.88e-10.
Peaks surge.367880m/sway1.549284m/heave1.520261m/roll5.12711deg/pitch1.72296deg/
yaw2.19986deg. Roll>=5deg ONLY final60.00s sample (peakcriterion, not sustained
5-10degstatistics); whole60 motion oscillatory/no oldrunaway. Fullsixplot:
artifacts/reference/researcher-review-20261024/jan90-seeded-six-motions.png.
Finalphysicaldisplacementis wavephase, not itselfseculardrift. ActualminusQ40
finaloffsets[.1292,.05936,.000116m,.09381,-.00948,.2694deg].

Source run8d8b8115/preparef5d07823/adapterd33cdfa5/verifierb783e5f4 andfixturePASS
helpers; same384bank/17operator/60s/dt.02/CIC40. Ardefaultzero/extrap UNCHANGED;
initialactualminusmeanlinear accelerations and nativeinitialwrenchresidual are
explicitlynonzero (max6404.65N/124.45Nm), uncompensated/fullstartupretained.
Not exactnativeperiodicity, not nonpassivityrepair. Exact7044/30896 terminal,
frozenIDs unchanged, separateidle-underleaseb436e0 and finalrelease complete.
MotionSHAd69652e3ca838afd2df6ddc249762d58321d833324038c51798399a36d055a95.

ONEcandidate retained for geometry/visualization ONLY; no finalphysicalaccuracy
acceptance. Geometryattempt1 validINCOMPLETE after300s:123/3001poses0..2.44s,
122area/rimtargetsmet; remaining2878 UNKNOWN. Fixtures/independentprefixintegrity
PASS, meanwetarea13.97955m2; prefixmaxXOR lower/upper9.7087/9.8047%, dryrimlower
.46776m, curvaturedeparturealready>=.23573m. No full60linearityclaim. Exact33716/
3436 terminal/idle/separateverify/releasecomplete. Originalsource/policy preserved.

BATCH recovery SOURCE independently PASS runner61c3598f/core3b04e8f2/configa411d7aa;
432oldscalaroraclefixtures PREPARED butUNEXECUTED. Newreportingprecision1%meanwet
area/10mmrim is numerical, notphysicalgate. Prospective1CPU600s/10Mtri/depth8/
65536active/2048wavebatches, all3001times orUNKNOWN. Oldmatchedprefix/newenclosures
must overlap. Linear JSONL/checkpoint reporting avoids quadratic writes; onefinal
fullJSON includesserializationelapsed. NO recoveryjob grantyet; Capawaitingroot
checkpoint+grant. No NLforces.

Preview1 COMPLETE2.344s(render)/6.225s(guard), exact31792/9188terminal/leasefree;
three literal0/30/60s images+all3001sixplot. All751displayframe clipping/phase/data
checksPASS. Rootviewedall3; source/data faithful, closer10mviewnowrequested.
Rootprovisionalvisualizationreceipt6cdb258e andartifactreadycontrols2811d87b:
case-jan90-seeded-attempt1/render-motion-acceptance.json and render-ready-controls.json.
ReceiptEXPLICIT geometry incomplete/nonpassivity/peakonlyendpoint limits, not
physicalqualification. Source6f3f5626 used12mXYtracking/fixedz[-3,3]; fixedcamera
alreadymeanheading-135deg, notinstantyaw. Preview2 COMPLETE1.797s(render)/6.098s(guard), source82083d9b/readyv2controls861b1282.
Root+independentliteral0/30/60 andprovenance PASS, all751clipchecks, phase2.44e-15.
Source/data unchanged; exact31120/35800terminal/released. Rootpreviewreviewreceipt
render-preview-review.json SHA119b3893 binds preflightba7d7f9f. No encodergrant.
Nextserialqueue:batchgeometryrecovery?encoder ONLYafterexplicitrootallocation.

Allstate/history/body sourcesandrenderer1 are committed/pushed e5dc02a; issue30
updated. Newbatchsource andnewrendererbindings notcheckpointed yet. Rendererowner
/root/generated_video; geometryowner/root/briefing_readiness; bodyownerfinished.
Topwater solver/NLadapter remainHELD; old107degcase REJECTED; no autocrop/rerun.

Completed24call SAME68nativepartial-lid control (55.1033s/terminal0) also remains
diagnostic:83lidfaces/77.27%waterplane, physicalhull unchanged, Ainf unchanged.
Heave/pitch B(12) flippositive but sway/yaw remainnegative; stillnonpassive. No
kernelrepair/furtherBEM/refined-motion grant. NativeinfinityAPI iscorrect; coarse
finite low_freq vs infinity high_freq Green representation remains unresolved.

Fullfootprint actualpose geometry evidence remains needed before applicability
claim.2.5deg freeflatwaterheel gives5.96%secant/12.83%tangent departure and8.53%
wetsetturnover; these are context, not adopted acceptance gates or combined
heave/trim/curvature bounds. No nonlineartrajectory model/adapter requested.

After resolving the result-invalidating growth, choose ONE crediblecase and verify
relativeimmersion/slopes/curvature with the actualnativeWEC trajectory. One REF002 bodyrun failed the operating envelope; no accepted canonicalselection yet. Stationary generatedinput is not the old
calm-start packet or an initialized nonlinearfluid. No custom integrator/NLmodel.

## Execution and preservation

[Execution ledger](experts/execution-state.json) owns current liveness. All
baseline numerical/encoder/decode jobs terminal; generator exactPID1960 terminal
and all producer checks complete, lease released. Direct Python exit code was
null after Windows controller Refresh; do not relabel it exit0. Controller0,
complete stdout/empty stderr/terminal identity and output checks are preserved.
Shared MATLAB9168/launcher37656 is historical identity, not liveness proof.
Follow [claim/recovery](experts/coordination.md); one serial owner, actual files
frozen until exact process terminal. Preserve other worktrees, original exporter,
all failures and interrupted quadrature36762 under original labels/policies.

Coordinator owns integration/handoff; Capytaine geometry/nativeproducer; WEC-Sim
consumer/dynamics/video; Topwater expert generator only; researcher independent
geometry assumptions. Baseline+generatorcheckpoint8d3755b and geometry/nativeproducer/consumer/Jan90
checkpoint3d5cedb and rejectednativecase/lidcheckpoint0c0c9e2 are pushed; issue30
updated through calmcontrolcheckpointd09c354; geometrysourcecheckpointbc8dc4f
is also pushed. NewICproposal/handoff refresh pending next checkpoint. Actualbodycase and24call lidcontrol are terminal; no numerical grant remains.
Geometryaudit and generatedrenderer independently SOURCEPASS, stillUNBOUND/nojobs.
Renderer source6f3f5626/controls7212083b,12.5fps751nativeframes/60.08s,equalmetres
12mXYtracking/fixedz. BindonlyNEWusablemotion, previewbeforeencode. Owner for
newrendererfiles: /root/generated_video; WECinitializer: /root/wecsim_reference. Never push main/create bookkeeping PR/merge/upstreampost.
