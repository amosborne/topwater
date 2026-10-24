# Checkpoint handoff

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
