# Checkpoint handoff

## Latest checkpoint: paired body-only uniform family rejected

Read docs/paired-body-study-001.md and docs/paired-body-protocol.md first. Branch
baseline/fluid-body; numerical source223050d. Goal ACTIVE/INCOMPLETE. No qualified
Capytaine reference, complete nonlinear RHS, TD motion run or usable speed exists.

- Controlled reflected half-seed, no matrix averaging/zeroing. Exact physical
  canonical hull/CG/full dryM/rho1025/g9.81/infinite depth preserved. Nested body
  counts54/216/864/3456; all7facets/moments/conforming edges verified. Actual
  Capytaine radii checked at every level, both resource preflights passed.
- Complete paired-body-study-002:4Topwater levels/24Capytaine calls, elapsed9.930s
  offline. Exact session17509 terminalexit1 by design, completeTrue/passedFalse.
  All18even entries fail reference prev/latest and candidate mesh/comparison;
  trendfails[0,2],[0,4]. Maximum ratios2091/1150(ref) and63.69/309.45(candidate).
  Failure masks/ratios and raw36entry matrices/U retained. Original strict gates
  unchanged; same-mesh correlation is not convergence. Previous fine reference
  also unqualified and used only diagnostically.
- Odd-parity normalized couplings reduced toCpt3.34e-7 /Top1.08e-14. This addresses
  seed triangulation symmetry, not full accuracy. All6χaux at every level retained
  for both solvers, normalized Cpt sources/bc retained. All8 complete A matrices
  reconstruct from savedχ; originalS/K checks pass beforeLU. No finite-wave or
  complete reference promotion. The general pressure-block producer still needs
  to save all7 known/aux fields on its next revision.
- Poor shape persists:48.15%panels aspect>1000 /0.2116%area. Exact matrix-increment
  accounting in paired-body-localization-001: latestCpt weighted increment68.7%
  below.1m depth,0.0172%in that sliver subset (Top67.2%/.1403%). This describes
  measured load increments only, not continuum-error location or harmless slivers.
- Preserved study001 API failure beforeCpt calls: wrongfaces_radii name, oneTop
  level complete, overallincomplete. Corrected faces_radiuses checked in freeze002;
  freeze001/002 numerical arrays bitwise equal. All native handles terminal.
- Durable evidence validation/evidence/paired-body-001.zip, summary nearby. SHA
  052df7c50865340e7ef4ebaddfc48adfef191b24057646bd395a889515dc98f2.
  Includes both freezes/studies, failedattempt, complete sources/physical inputs,
  auxiliary fields/failure masks, exact accounting code, JUnit and parent payload.
  Run restore_shared_inputs.py after extraction; five input paths SHA-verified.
  Fresh extraction/owned graphs checked. All261tests passed, exact session19563
  terminalexit0. Independent audit rebuilt matrices/masks/accountings; no blocker
  for preserving this failed investigation. No live numerical jobs.

**Exact next bounded diagnosis:** separate frozen-source pressure quadrature from
remaining discretization before remeshing. Use existing864-panel normalized Cpt
indirect sources and same physical surface, evaluateS at strict interior descendant
centroids (prospectively declare e.g.4/16/64subcells, all36 loads, limits/stabilization).
First reproduce native centroid loads, then preserve exact telescoping
A_f-A_c=(A_Q-A_c)+(A_f-A_Q). The remainder also contains finest-output quadrature
error; never label it uniquely pure source error. Topwaterχaux is NOT an indirect
source density. No newBEM solves or runtime/formulation/optimization campaign.
Full finiteA/B/excitation/phase, initial compatibility, nonlinear/time/domain and
complete runtime gates remain unresolved. Earlier pressure/field context follows.

## Prior checkpoint: unqualified P0 pressure component and error attribution

Read docs/mixed-pressure-001.md and both prospective protocols first. Goal
ACTIVE/INCOMPLETE: canonical boat, own nonlinear time-domain six-DOF response,
>=1s, <=1cm /0.1degree including uncertainty, <=200wall s/s including preparation.
No complete nonlinear RHS, motion run, usable speed or qualified reference exists.

- New runtime owners src/topwater/fluid_boundary.py and pressure.py use the same
  P0 Rankine/boundary path. Generic convex wet-facet graph/ear triangulation keeps
  contact vertices/concavity, exact physical facets and manifold closure; artificial
  hull closures carry no pressure. Owned pressure factors cannot pair another pose.
  Body-normal chi terms include rotation and Hessian, which remain mandatory inputs.
- Producer6213d96 mixed-pressure-001 completed, passed=false/exit1 deliberately.
  Calm exact-rest probe has762 panels(48body); preserved mild geometry closes2616
  (151body) but receives NO pressure call with invented zero derivatives. Calm
  accelerationmax3.12e-15 and combined pressure load independently reconstructs.
  All six raw A columns saved, but NOT six chi_aux traces: retain all seven chi
  fields in the next producer before claiming artifact-level rawA reconstruction.
- Timing: warmup.279s once; calmgeometry.100s/assembly-factor.073s/pressure.0062s.
  Projection36.17wall s/s for pressure component at provisional midpointdt.01 only.
  This excludes remaining derivative/RHS work and is NOT a speed pass.
- Added-mass diagnostics fail all four36entry comparisons. Attribution902c2be
  pressure-discrepancy-001 completed six Capytaine infinity-frequency calls on the
  exact48body panels; seven facets/edges/moments proved, FFF order. Topwater
  halfspace vs same-body Capytaine differs5.1–26.0% diagonally; closed vs halfspace
  differs5.3–19.2%; same-body Capytaine vs finest differs-56.3% to+30.3%. Finest
  reference itself remains unqualified. Multiple unresolved resolution/formulation
  and surface/domain effects, not an established absolute-error attribution.
- Complete evidence validation/evidence/mixed-pressure-001.zip and summary. SHA
  ec77f90b85b6a101efa81b65cf714100ab95b9c4f7ca44ce4b9f20ce66914236.
  Shared parent field payload appears once; run bundled restore_shared_inputs.py
  after extraction to restore both original paths with verified SHA. Fresh extraction
  and both complete producer/input graphs verified. Independent review rebuilt all
  144 comparisons and combined pressure load. All255 local tests passed, JUnit
  saved; exact test session32046 terminalexit0. Both native jobs terminalexit1,
  completeTrue/passedFalse, no live numerical processes. Branch baseline/fluid-body.

**Next bounded validation:** assess independent body resolution before changing
surface/domain or optimizing. Design/preflight controlled nested refinement of
the same physical wet facets (including all contact/feature edges), then run a
prospectively declared paired Topwater-halfspace/Capytaine study. Actual counts,
memory/time limits, all36 entries/uncertainties and at least three/four useful
levels must be declared before solves. Same-mesh agreement is not convergence.
Retain original reference gates; do not relax them. Preserve all auxiliary chi
fields in the next pressure producer. Body-interior resolution is currently tied
to sparse wet-facet triangulation; current closure is not a qualified spatial mesh.
Later address surface/domain controls and derivative recovery, then physical
motion comparison before major optimization. Previous field/source context follows.

## Prior checkpoint: exterior fields and surface-map candidate

Numerical producer **05cc94e**, branch **baseline/fluid-body**. The goal remains
ACTIVE/INCOMPLETE. No own coupled TD run, complete restart or usable speed exists.
Read docs/periodic-field-map-001.md and docs/periodic-grid-protocol.md first.

- periodic-fields-001 completed all 12 evaluations in 124.014 offline seconds.
  Exact session 86292 is terminal exit 1 intentionally: complete=true/passed=false.
  All 1,518 reference points are strictly exterior at z0; archived sources reused,
  no new radiation/diffraction solves. Maximum Green block 100 by 7,700.
- Full/half/quarter-amplitude maps pass geometry only: 145 contacts, seven inserted
  facet switches, 2,030 triangles, positive projected areas, no whole-triangle
  hull-interior overlap, every contact edge on a hull facet. Conservative opening
  clearance >=0.584 m. Mapped scalar potential is a pullback, not a compatible state.
- Both prospective source-mesh and waterline-limit trends FAIL. Combined field
  envelopes <=11.076 micrometres eta /0.000051836 m^2/s phi do not waive the failed
  trend requirement. Total angular envelope reaches 9.667 mm /0.034599 m^2/s.
  Full reference/restart, body compatibility, domain and nonlinear-limit flags false.
- field-sampling-001 is an archived array-only diagnosis: outer incident eta
  envelope 9.667406 mm versus total9.667362 mm; scattered maximum across all rows
  0.065735 mm /0.000301874 m^2/s. The angular defect is mostly interpolation of
  known incident waves. Preserve those analytically when correcting representation;
  do not blindly refine the far grid or retroactively pass the failed diagnostics.
- Complete evidence validation/evidence/periodic-fields-001.zip, summary nearby;
  SHA cc35d61aa1febeb6b9389d9f756304fa15a411fdbd14c01b6ae932eeabc45d65.
  Includes actual producers, controls, physical inputs, nested source/raw/freeze
  evidence and diagnostic code/identities. All 246 local tests passed; exact
  session88691 terminal exit0. No numerical jobs live. No runtime Python changed.
- Independent review rebuilt field signs/phases, all envelopes/trends, maps,
  contacts/whole-triangle overlap and source identities; no blocker for preserving
  this explicitly unqualified candidate. Original coefficient failures unchanged.

**Next bounded work:** correct representation of the analytically known incident
waves, then construct the minimum complete body-fitted fluid boundary and selected-
P0 pressure/acceleration compatibility calculation needed to measure initialization
sensitivity in motion units. Declare equations, normals, gravity counting, error
controls and a stopping point; review before numerical jobs. Body-boundary
projection, complete domain and nonlinear-limit evidence remain missing; mean
stern chine submergence is only ~2.8 mm. Radius/depth24 m are candidate controls,
not a qualified truncation; 32/48 m remain future domain checks. Use actual coupled
prep/RHS cost to reject impractical candidates before a long run; never infer
motion accuracy from a small residual. No major optimization/global refinement
campaign before physical validation. Earlier physical/source details follow.
One initial acceleration difference is not a one-second motion-error bound;
an unvalidated RHS sensitivity cannot be the sole qualified reference uncertainty.
Exact analytic incident fields also do not eliminate total-surface geometry error.
At0.1 Hz, k=0.040243/m and exp(-kH)=0.3807 at24 m /0.1449 at48 m. Do not use
a zero-flux24 m bottom as if it were infinite depth. Outer/bottom truncation and
their initial compatibility must be declared and checked before coupled work.

## Prior checkpoint: periodic source recovery

Source checkpoint **785c3d1**, branch **baseline/fluid-body**. Read
docs/periodic-source-recovery.md and docs/periodic-source-protocol.md first.
The broader goal remains ACTIVE/INCOMPLETE. No Topwater coupled TD run exists.

- Recovered all required Capytaine basis sources on the four unchanged meshes:
  84 solves and 384 reproduction checks passed, offline elapsed 329.246 s.
  artifacts/validation/periodic-sources-003.json/.npz/.source is complete.
  Exact session 32728 is terminal exit 0. This is not Topwater speed evidence.
- Saved source, hull potential/pressure, prescribed normal velocity and produced
  S/K traces for all seven bases per frequency. Raw S/K retained before cached K
  becomes LU. Independent pressure integration and geometry-derived rigid modes
  check signs/world-CG moments; all finite A/B entries and selected complex loads
  reproduce the pinned study. Largest tolerance ratio 0.0405471.
- Full reference qualification and nonlinear restart readiness remain false.
  No nonlinear-grid eta0/phi0 exists. Original coefficient-convergence failures
  are unchanged. Basis sources are Capytaine indirect densities, not interchangeable
  with Topwater Rankine densities, and never runtime force/trajectory drivers.
- Complete evidence: validation/evidence/periodic-sources-003.zip, SHA256
  4a317cbb91cc51348b0059b40b4a6437315e3656624ad2903743dec4721214e5; summary nearby.
  Includes complete nested raw/freeze evidence and unsuccessful attempts 001/002,
  each with zero completed physics calls. Session 50564 was interrupted for review
  corrections; session 69350 failed the symbolic-DOF array API. Both terminal exit 1.
- All 238 local tests passed (session 2470 terminal exit 0). The only post-run
  Python edit is an optional-Capytaine API-test skip guard; numerical producers
  are unchanged. Actual recovery still fails without the dependency. Source
  hashes were verified before that test-only change; archived producers preserved.
- No numerical jobs live. Preserve serial-job and immutable-producer discipline.
  Review and promotion scope is source reproduction only, never full qualification.
  Independent review rebuilt all basis/composed loads and phase/normal checks,
  verified identity/bundle/mesh evidence and found no remaining stage blocker.

**Next bounded goal:** freeze the minimum TD grid/domain and an admissible
reference-to-nonlinear geometry map, then sample the archived sources (avoid a
new 84-solve batch if possible) into total eta0/phi0 with measured field/mapping
uncertainty. Reject points inside the mean hull or above the reference fluid.
Check actual body-boundary compatibility, projection transients and amplitude/
wetting controls in motion units; a small boundary residual is insufficient.
Mean stern chine is only about 2.8 mm submerged. After a compatible packet exists,
screen actual coupled preparation/RHS cost before 1 s; reject a route projected
above 200 wall s/s including normal preparation. Major optimization stays held.
Declare controls/acceptance/stopping point first, review, checkpoint/push/update
issue #30, then write a fresh handoff. Never push main or merge without CI.

## Prior checkpoint context (retained physical/motion details)

- **Active outcome:** Topwater TIME-DOMAIN ONLY, canonical 24 ft boat in a frozen
  mild sea, small-motion nonlinear six-DOF run for >=1 s; errors plus uncertainty
  <=0.1 degree per rotation /1 cm per translation; <=200 wall s per simulated s
  including normal Topwater preparation. Read docs/time-domain-baseline-goal.md.
  Major optimization waits for a passing physical comparison. Early bounded cost
  screens reject impractical routes. Goal ACTIVE/INCOMPLETE; no own TD run exists.
- **Branch/current source:** baseline/fluid-body. This stage's source checkpoint
  a5ebf7e; inspect current git status/log and full manifest identities. No Topwater
  runtime Python was changed: this is external reference-only reconstruction.
  No PR, no push to main, no other worktree cleanup.
- **Current implementation:** unqualified P0 Rankine-panel core from 4e8f511,
  hydrostatics/instantaneous inertia only. No coupled TD wave path yet. The archived
  P1 runtime is a different implementation; recover generic math/tests selectively.
- **Frozen mild case:** validation/mild-periodic-case.json, case_id
  canonical-boat-mild-periodic-001. Controlled three-tone sea: amplitudes 20/10/5 mm,
  frequencies 0.10/0.25/0.50 Hz, propagation 0/45/90 degrees, phases 0, pi/2, pi,
  equivalent component Hm0=0.0648074 m. NOT an observed buoy realization. Full
  amplitudes/phases are frozen before any Topwater comparison. Amplitude 0.5/0.25
  checks are for linear-limit assessment, not replacement of the full case.
- **Initialization convention:** periodic prehistory, not rest-start. R=Exp(theta)R_eq;
  theta=Log(R R_eq^T) in mean world axes, omega_body=R^T J_left(theta) theta_dot.
  Body packet has exact pose and frame-correct velocity only; analytic acceleration,
  RAOs and future traces stay external diagnostics. Total incident/diffraction/
  radiation eta0/phi0 and a compatible nonlinear geometry map are still missing.
  Projection can launch adjustment transients, all charged against motion limits.
- **Prospective acceptance:** docs/mild-periodic-protocol.md defines the new
  forcing-specific output budget, not a relabeling of failed coefficient gates.
  Sum max(abs(P4-P3),abs(P3-P2)) over complex wave components for every physical DOF,
  plus reference sampling bound sum(abs(omega*P4))*dt/2, must fit 10% of user limits
  with declared non-increasing component trend. Constant-initial-pose signal minus
  twice that envelope must exceed user tolerance in a translation AND a rotation.
  This does not bound Topwater interpolation, initial fields, compatibility or
  nonlinear representation; those remain required before complete readiness.
- **Executed result:** artifacts/validation/mild-periodic-001.json/.npz/.source.
  Complete, overall passed=false and exit 1 by design because four readiness flags
  remain false. Motion-spatial/sampling gate PASSED, maximum reference-budget ratio
  0.4409778; informative-case check PASSED. Frozen-pose lower error exceeds target
  in surge/heave/roll/pitch; sway/yaw are weaker and require separate mode checks.
  No pressure/BEM/Topwater TD solve ran in this stage.
- **Durable evidence:** validation/evidence/mild-periodic-001.zip and companion
  summary/SHA256; includes all histories/derivatives, complete raw matrices and
  increments, source archive, body/wave-only JSON and nested complete raw-reference
  bundle. The body-wave packet is SimulationInputs v2, NOT a complete fluid restart.
  Missing periodic_fluid_fields/field_uncertainty/nonlinear_compatibility/linear_limit
  flags are all false. Full reference and actual baseline remain unqualified.
- **Checks/review:** 60 local validation tests passed: phasor sign/derivatives,
  non-canceling uncertainty, sampling allocation, finite-difference body angular
  velocity, conservative negative control and changed-case rejection. Independent
  fixed_topology_review rebuilt the artifact/provenance and all stated checks;
  no blocker for adopting this limited stage. It explicitly did NOT qualify a
  complete reference or simulation. Numerical producers stayed immutable while run.
- **Underlying data:** original study 003 (four meshes/all 132 calls) remains FAILED
  under 0.1%+1e-6 coefficient policy, including all six Ainf diagonals. Raw input is
  prospectively SHA-pinned in the case. Geometry is 308-panel graded seed, affine
  N=2/3/4/5, counts 1,232/2,772/4,928/7,700, same frozen physical surface/CG/full M/K,
  rho 1025 / g 9.81/infinite depth. Canonical source/Green settings/old gates unchanged.
- **Historical math:** b568edd:tests/capytaine_time_reference.py and analytic tests
  at 135c430. Cummins signs/endpoints sound, but old depth 4/1000 kg/m3 canoe pin is
  inapplicable. It computes radiation for supplied V/a, not coupled excitation/motion.
  No invented broadband IRF is needed for this exact three-tone periodic case.
- **Exact next bounded work:** freeze minimum TD fluid grid/domain and an admissible
  reference-to-nonlinear geometry map BEFORE claiming physical initial-state readiness.
  Reject reference points inside mean hull/outside fluid; no silent potential
  continuation. Then export/compose total periodic eta0/phi0 with source identities,
  field/mapping uncertainty and body-boundary compatibility checks in motion units.
  Capytaine indirect source fields are available via solve(keep_details=True);
  selected three frequencies need all six radiation modes plus the selected heading
  diffraction, with full phase/units verified against retained raw A/B/excitation.
  No future Capytaine pressure/body trajectory/acceleration may drive Topwater.
  After a complete compatible packet, run a bounded preparation/actual coupled-RHS
  cost screen before attempting one second; primitive timing alone proves nothing.
- **Stopping checkpoint:** compatible initial-state evidence and minimal TD feasibility
  result, or a concrete rejected candidate with unchanged motion limits. Declare
  controls before jobs; archive immutable producers/inputs, run serially, verify
  exact handles terminal, resolve review, commit/push and update issue #30/handoff.
  No major optimization or alternate-formulation campaign before actual validation.
- **Live jobs:** none. Array candidate command terminal exit 1 (expected readiness
  block), tests terminal exit 0. Earlier study/observer also terminal. No TD speed
  has been measured. CI remains PR-to-main only; no merge without current-head CI.
