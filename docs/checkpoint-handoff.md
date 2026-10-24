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

WEC owns SOURCE/PROPOSAL ONLY for60s/dt.02/40smemory held-release case: prospective
q0 follows deterministic localtangentplane, zeroallvelocity/radiationmemory;
priorholdingload and initializedincident/diffraction assumptions explicit. No
nativebody launch allocated yet. Cap owns source-only one24-call SAME68hull
native numerical-lid control at1,2,12,infinity; no launch yet. Source check confirms
native omega=infinity API is correct; default finite low_freq and infinity
high_freq Green decompositions differ numerically on coarse meshes. No repair.
The lid control could test irregular-frequency mechanism; not fullkernel adequacy.
Fullfootprint actualpose geometry evidence remains needed before applicability
claim.2.5deg freeflatwaterheel gives5.96%secant/12.83%tangent departure and8.53%
wetsetturnover; these are context, not adopted acceptance gates or combined
heave/trim/curvature bounds. No nonlineartrajectory model/adapter requested.

Then choose ONE case/initialheading/explicitduration+initialstate and verify
relativeimmersion/slopes/curvature with the actualnativeWEC trajectory. No actual
REF002 bodyrun or canonicalselection yet. Stationary generatedinput is not the old
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
geometry assumptions. Baseline+generatorcheckpoint8d3755b is pushed and issue30
updated; geometry/nativeproducer/consumer/Jan90 screen form the next substantive
checkpoint; source-only follow-ups remain separately unallocated. Never push main/create bookkeeping PR/merge/upstreampost.
