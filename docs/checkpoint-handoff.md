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

Next bounded numerical job: capytaine_expert geometry operating-envelope study,
fullsamehull exactGL2, heel/freeheave-trim/righting plus wetset/waterplane/contact
curves. Source/manifest independent review underway; do not launch without sole
coordinator allocation. Proposed5%secant/10%tangent labels are exploratory
engineering proposals, not adopted dynamic accuracy gates. Interpret displacements
relative to local water level/slope; large global wave-following heave alone does
not prove immersion/flooding. Then bind WEC producer frequency/heading coverage
and faithful linear generated-bank consumer to select ONE actual strong6DOF case.
WEC owner prepares explicit initialization/duration; stationary generator output
is not the old calm-start packet or an initialized nonlinear fluid field.

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
geometry assumptions. Baseline/follow-up checkpoint commit/push and issue30 update
are integration work next. Never push main/create bookkeeping PR/merge/upstreampost.
