# Owned bounded assignments

Task statuses are `queued` (owned but not authorized), `authorized` (within the
current user request), `running`, or `completed`. Current project authorization
is recorded only in [the handoff](../checkpoint-handoff.md); an execution lease
schedules authorized work and never grants scope. Update status in place.

## Radiation-surface height accuracy, 2026-10-27

User acceptance for #31: maximum numerical elevation error plus reference
uncertainty <=10mm against the chosen linearized model, full60s and declared
valid exterior output region. Saved motion/forcing and hydrodynamic data stay
fixed. Internal transfer/derivative diagnostics do not replace this requirement.

- CAP-ISSUE31-HEIGHT: screen completed; next PRIVATE helper authoring authorized,
  capytaine_expert.
  Seven-probe all-six saved-motion screen: reproduce production cosine response
  and compare full-complex frequency reconstruction under prospective reference
  controls. Next bounded proposal extends the SAME chosen response from24 to48
  rad/s, seven probes/all six modes, with overlap and input-weighted tail checks.
  No body/mesh/source-method correction, motion rerun or production edit allocated;
  extending response support does not change the chosen hydrodynamic model.
  Own helper and concise evidence/plan note; root owns numerical launch.
- NUM-ISSUE31-HEIGHT: independent numerical_methods_researcher review of reference
  design, source/API/physics conventions and actual height/uncertainty results.
  No physical implementation authorship or numerical execution delegated.
- Root owns scope, prospective controls, serial execution, acceptance, integration
  and checkpoint publication. Seven probes cannot establish the whole-region
  maximum; unresolved reference/tail uncertainty must remain visible.

Prior diagnostic studies and their failed gates remain preserved. No correction
has yet been accepted into production. Speed and #32-34 remain later work.

The following assignments are historical and confer no current execution authority.

## Canonical workflow implementation, 2026-10-25

Current user authorization is the canonical NPZ/directory CLI plan in the handoff.
Coordinator owns CLI, dependencies, documentation routing, numerical scheduling,
acceptance and Git checkpoints. No native dynamics implementation is allocated.

LOG-WORKFLOW-001, authorized/source implemented: coordinator owns CLI wiring,
documentation and integration; Topwater expert authored command_log.py and its
tests; independent researcher triages actual warnings. Deliver sea.log,
result.log and video.log with consistent timestamps/levels, failure/publication
status and latest-attempt replacement. Preserve numeric NPZ contracts. Current
228-test suite passes; actual complete case/log/video acceptance remains pending.
User declined streaming solve/render; commands remain independent.

- TW-WORKFLOW-001, authorized/source implementation: topwater_expert owns new
  boat_geometry.py/workflow_data.py, corresponding tests and mesh-tooling review.
  Deliver convex one-opening physical hull validation, NPZ contracts, temporary
  closure adapter and bounded geometry probes. Do not alter generic old inputs.
- WEC-WORKFLOW-001, authorized/source implementation: wecsim_expert owns reusable
  backend/setup, equilibrium/train construction, native WEC runtime and solve
  documentation. Deliver generic input-driven pipeline without frozen artifact
  paths or another motion solve during surface reconstruction.
- CAP-SURFACE-001, authorized/source implementation: capytaine_expert owns reusable
  hydrodynamic producer/retained fields and total transient surface reconstruction,
  with producer tests. Joint interface agreement with WEC; independent review
  before physical reliance. No mesh-refinement campaign.

All experts preserve others' edits and stop before tests/numerical runs until the
coordinator explicitly assigns the sole execution lease. Package/dependency and
producer identities are frozen during every allocated run. Geometry/physics
authors do not independently approve their consequential corrections.

## CAP-001: first canonical producer/interface

Status: completed (minimum export and stored BEMIO interface, 2026-10-24).
Owner capytaine_expert; consumer partner wecsim_expert. Native17-panel export,
independently reviewed GL2 hydrostatics and stored official BEMIO import signoff
are in capytaine_expert/minimum-producer-2026-10-24.md and
consumer-interface-review-2026-10-24.md. Runtime body mapping/startup belongs to
REF-001 and is now completed/audited there. No refinement is authorized by this completion.

Scope: make the minimum-mesh native Capytaine export needed by the first WEC-Sim
case, preserving exact physical surface/full inertia/CG and agreed DOF/heading/
phase/units. Coordinate frequency and transient coverage; retain native files
and exact identity. Exclude Topwater code, mesh campaign, intermediate coefficient
acceptance thresholds, and unauthorized fixes to the consumer.

Deliver compatible native files, small jointly approved producer/consumer
contract and reproducible command/manifest. Acceptance: consumer can read the
physical and numerical data without missing/defaulted metadata, phase or unit
ambiguity; no coefficient threshold is required. Freeze controls before launch.
Stop at reviewed minimum-case export/interface; no automatic refinement.

## WEC-001: maintained fork, independent leads and exercised paths

Status: authorized (2026-10-24 minimum-mesh reference goal). Owner wecsim_expert. Fork/source-map/setup inspection completed
in wecsim_expert/native-setup-2026-10-24.md. Remaining reproductions are authorized
within the minimum-case exercised paths, under a sole numerical execution lease.

Scope: pinned v7.1.0 maintained fork; independently reproduce the four specific
marine-gyro leads (yaw-order/frame and velocity/acceleration mapping, unequal
postprocessing grids, angular-bin energy normalization, direction row/column
shape) using independent expectations. Inspect the canonical producer/BEMIO/
startup/radiation/hydrostatic/frame/output/video paths actually exercised.
Review source map first; choose the smallest decisive reproductions. No full
library audit, automatic carry-over patches, string-replaced generated models,
Topwater changes, broad refactor or upstream posting/PR without authorization.

Deliver per lead: selected/affected versions, exact paths, minimal frozen
reproduction, independent expected versus actual, physical consequence, cause
confidence, confirmed/unaffected/unconfirmed status, focused regression and
minimal generic fork correction ONLY if confirmed. Case-specific configuration
stays here. Missing evidence is not a confirmed defect. Record exact pins and
remaining limitations in the KB; concise report to coordinator and Capytaine
partner. Independent reviewer distinct from author before physical reliance.

Stop at bounded lead/exercised-path findings and reviewed minimal corrections
or explicit limitations. Reproduction commands/resources/expectations must be
declared before their execution lease. Result-invalidating confirmed issues
must be corrected or the affected result withheld. Speculative unrelated issues
must not delay the first working minimum-mesh10s simulation indefinitely.

## TOP-001: held-state ownership

Status: completed (setup only; implementation held). Owner topwater_expert.
Activation/documented equation/code/evidence map completed
in topwater_expert/held-source-map.md. No solver implementation, simulations,
test runs, optimization or method experiments. Stop pending future user authority.

## NUM-001: focused research support

Status: completed (setup only). Owner numerical_methods_researcher. Proposal template and targeted initial
wave/IRF-padding question completed in its KB. No broad literature/optimization
survey or numerical lease. Next only a focused question requested by reference
owners; proposed experiments go through coordinator selection and implementation
owner. Stop at supported answer and cheapest rejection plan, not a new solver.

## Downstream milestone (not this setup's completion criterion)

The first freely floating WEC-Sim ten-second incoming-wave motion and isometric
video, followed by motion variation with mesh. Keep shared forcing/initial state
and exact producer identity; no PTO/mooring/invented constraints. Coefficients
remain diagnostics. Setup, fork creation and role confirmations do not claim
that simulation/video or the four defect reproductions already exist.

## REF-001: minimum-mesh end-to-end reference

Status: completed and independently audited 2026-10-24; checkpoint8d3755b pushed. Coordinator owns acceptance/integration;
capytaine_expert owns producer and wecsim_expert owns consumer, dynamics/video.
numerical_methods_researcher supplies independent startup/interface review.

Complete CAP-001/WEC-001 work needed for the frozen canonical minimum-mesh
case. Exclude mesh refinement and Topwater implementation/tests/optimization.
Required completion evidence:

- Reproducible native Capytaine files, official BEMIO H5 conversion, exact
  producer/input/dependency identities and joint executed interface signoff.
- Canonical hull, equilibrium/world CG, mass/full inertia, physical constants,
  headings/phase/units and six DOFs mapped without defaulted physical metadata.
- Ten seconds of freely floating six-DOF WEC-Sim output, native time grids and
  rotation representation retained, initialized at actual rest. No mooring/PTO
  or invented constraints; no discarded startup or imposed temporal ramp.
- Frozen upstream incident packet shared by elevation/forcing/display, with
  documented arrival and convolution prehistory/padding. Verify locally calm
  initial water and force history; no moving wavefront over globally forced
  motion. Label visualization as prescribed incident surface.
- Isometric video showing incoming waves and boat motion over the interval,
  plus six motion histories and finite/output-duration/initial-state checks.
- Independent physical review of consequential corrections and startup mapping;
  confirmed result-invalidating defects resolved or affected result withheld.
  Preserve failed runs and limitations; no coefficient convergence threshold.
- Development-branch checkpoint committed/pushed, issue #30 updated and fresh
  handoff with result paths, commands, limitations and exact terminal evidence.

Stop at the verified minimum case; subsequent motion variation with mesh remains
a separate follow-up. A compatible export or successful setup alone is incomplete.

## REF-002: one buoy-derived canonical strong-motion case

Status: completed for the minimum-mesh end-to-end workflow and requested stronger
motion/video,2026-10-24. Actual source Hm01.657m,60s,6.62deg peakroll and
appreciable other5DOF; all native/data/render checks pass independently.
Physics/mesh/time/nonlinear applicability is provisional, not a precision truth
qualification. See current handoff and stronger-case evidence. Preserve REF-001 unchanged.
Target one WEC-Sim case with roll about5-10deg and appreciable other5DOF.
Use the existing collected-buoy sea-state generator, genuinely higher-energy
observations and complete seeded directional realizations; do not scale the
baseline sea or invent a fitted spectrum. Internal candidate screening is
permitted; deliver one frozen canonical state/seed/initialheading, not a suite.
Initial box orientation follows the dominant incoming direction by default;
store propagation/from and global/body-frame conventions explicitly. Heading
offsets may be examined for choosing the one case, not silently adopted.

Linearity is an a-priori engineering operating-envelope assumption established
from full-hull buoyancy/righting/restoring, wetted material surface, waterplane
and contact curves against WEC-Sim's fixed mean assumptions. No nonlinear
trajectory model or native weak-nonlinear adapter is required. Quantify the
curves and label assumptions; do not claim nonlinear dynamic validation.
Same physical boat/mass/CG/inertia/constants; no damping/boat fitting to hit roll.

Ownership: topwater_expert owns bounded existing-generator neutral export ONLY;
Topwater solver implementation/tests/simulation/optimization remain HELD.
capytaine_expert owns geometry operating-envelope and required native producer
coverage. wecsim_expert owns faithful generated-bank consumption/linear dynamics
and final video. numerical_methods_researcher independently reviews envelope
math/assumptions. Coordinator owns selection/acceptance/serial execution.

First jobs: three genuine high-energy observations with seed42, whole-band
energy/replay/frame checks, and exact geometry static curves. Before dynamics,
freeze new sea, heading coverage, initialization/duration, numerical controls
and prospective operating-envelope interpretation. Do not silently inherit
REF-001's artificial test packet or assume a generator realization is calm water.


### REF-002 latest user steering (2026-10-24)

The full validation duration is60s. Final startup is TRUE locally calm water with
waves initially upstream, body rest, zero velocity and zero radiation memory;
viscous damping remains off. Rest in an established sea is insufficient. No hidden
force ramp. The nonzero stationary-history
prototype and its unpruned-bank geometry/preview work are preserved under their
original policies but superseded for the final requested case. Do not encode that
old prototype or resume its stopped geometry audit.

Prune unsupported incident bands only after checking the actual full radiation-
diffraction physics and numerical coverage. The generator already enforces the
project wavelength>=2beam envelope (~.53346Hz for9ft); current bands end.495Hz.
No FK-only approximation or universal beam-based diffraction cutoff is assumed.
Retained measured energy/MEM/phases must not be renormalized/fitted. Radiation-
kernel integration coverage and response transients need separate qualification.

Three native energy witnesses completed:17 coupled2.2 negative mean outgoing work,
17heave1 positive control,68samecoupled2.2 positive. Strong mesh sensitivity;
see producer witness receipt. Both official RM3 standard400s and provided
regularCIC200s regression also completed independently PASS; supplied CIC golden
errors<=4.157e-13m versus1e-5m tolerance. These do not qualify our boat operator.

Current diagnostic: new singlepacket preflight rejected by initial heave-force
gate; no body run. Selected unchanged accepted REF001 threepacket extension
to60s, with originalH5/CIC10. Incident extension independently EXECUTEDPASS:
exact spectralidentity/oldwaveprefix and first501six-force samples; full3001
forceconvolution independently reproduced. No retiming/scaling/new calm field.
Originaloperator body extension also completed32.182s and independently PASS:
bounded small motions, first10s q/v reproduced to roundoff. It is a diagnostic,
not the final strong-motion reference. Stored-data CIC work stays cumulatively
nonnegative in both this run (+7.390J final) and the earlier rejected Jan90
large-response run (+237.894kJ final). No actual net radiation-energy creation or
causal link to that response has been established. See the actual-radiation-work
review. Cap read-only next-step recommendation is complete:
a compact producer audit can reject, but only full6motion histories establish
mesh variation. No further producer execution allocated. No Cap/WEC radiation implementation
bug, valid cutoff or qualified replacement has been established. Topwater role
scope remains existing-generator data/export only; all solver work held.

Superseded user priority: progressively modify the passing official regularCIC
example. WEC owns bounded plan/first-step SOURCE readiness: sole physical change
PTO damping1.2e6->0, official H5/heave constraints/ramp100/CIC20/dt.1/200s retained,
native raw save/work accounting/baseline comparison. Researcher independently
reviews exact source/freeze before coordinator grants a single execution.
This PTO-zero proposal was cancelled and remains UNRUN. Later user regroup
delivered the mild and stronger incoming-buoy cases and videos. No suite,
new mesh/Topwater implementation or solver-energy bug campaign is allocated.
