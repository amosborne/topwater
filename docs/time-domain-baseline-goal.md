# Validated, practical time-domain baseline

Historical broader Topwater goal, now HELD by the user's instruction to work
only on Capytaine plus WEC-Sim. Current reference-only goal and motion-based
policy are in capytaine-wecsim-workflow.md; conflicting prerequisite coefficient
gates below no longer control that workflow. No Topwater code work is authorized.

The user clarified this goal on 2026-10-23. Reference qualification is a
prerequisite, not the final deliverable. Topwater's runtime is time-domain only.

## Scope and deliverables

Deliver Topwater's own small-motion nonlinear fluid/body simulation of the
canonical 24 ft boat in a frozen mild sea state, with all six rigid-body motions.
Run for at least one second and show matching Capytaine-based time histories,
six motion/error plots, per-axis maximum errors and uncertainty, and a measured
complete runtime. Provide reproducible commands, raw traces, source/input
identities, independent review and blocking physical regressions.

Complete the missing coupled wave-response capability in the selected P0
Rankine-panel time-domain path. Reuse sound existing math and tests selectively;
the archived P1 diagnostic runtime is a different implementation and supplies no
qualification evidence for this boat. Supported boats continue to change through
input data, with one supported runtime execution path.

Capytaine supplies the matched linear-limit external reference; it cannot run
the identical nonlinear formulation. Nonlinear-versus-linear differences count
against the same motion-error budget. Capytaine-derived loads/RAOs/Cummins
coefficients are reference machinery, not a substitute runtime fluid solver.

## Acceptance

- Use the same physical hull, pose/CG conventions, full dry inertia, density,
  gravity, depth, incident frequencies/headings/phases/amplitudes and declared
  body/fluid initial state in both models. Freeze the mild case before comparing.
  A scaled catalog spectrum must be identified as a scaled validation case.
- Over the complete declared comparison interval of at least 1.0 s, each surge,
  sway and heave error must be at most 0.01 m, and each roll, pitch and yaw error
  at most 0.1 degree. Reference and comparison uncertainty must be charged against
  these limits; a difference alone cannot pass while uncertainty is unresolved.
- Establish reference and Topwater spatial, domain, temporal, startup/initial-state
  and amplitude controls appropriate to this case. Declare numerical controls,
  uncertainty estimates and acceptance before judging new results. Do not fit
  phase/amplitude to the simulation, hide startup error or use an interval with
  negligible forcing as proof of wave response. Preserve six-mode/coupling checks.
- The mild case must support a small-motion comparison. Check approach to the
  linear limit by amplitude refinement; do not call a single mild-case result
  general nonlinear or full-minute accuracy qualification.
- The complete Topwater run, including its normal numerical preparation and
  initialization, must cost at most 200 wall seconds per simulated second on the
  development machine. Report offline Capytaine reference generation separately.
  Accuracy and practical runtime are both required for goal completion.

## Validation before optimization

Measure setup and a few steps before launching an expensive full run. Reject
obviously impractical routes early rather than spending thousands of seconds
per simulated second. Bounded implementation fixes needed to execute the
validation are allowed. Major architecture/performance campaigns remain held
until an actual coupled Topwater-versus-Capytaine comparison passes. Then retain
that case as a blocking regression during any required speed work. Speed never
permits a relaxed motion-error gate or an unvalidated physics substitution.

## Next bounded work and stopping checkpoints

1. Preserve study 003 as failed under its original coefficient policy. Select and
   freeze a mild, nontrivial validation case and matching initialization. Determine
   which reference histories that case requires, and declare an output-error
   uncertainty protocol before generating or accepting them. Do not automatically
   repeat blanket coefficient refinement without a motion-error-budget rationale.
2. Audit the existing reference-only time-domain machinery at commits 135c430
   and b568edd. Its generic kernel/convolution/back-transform math has analytic
   tests, but its depth-4 m canoe binding and unqualified frequency band cannot
   serve this deep-water boat. It evaluates radiation forces for supplied motion
   histories; it does not already supply a coupled motion/excitation solver.
   Reuse the math with explicit current physical
   identity and band/spacing/tail/memory/time-grid checks where transients require
   them. Steady harmonic reconstruction alone does not validate arbitrary startup.
   A periodic reconstruction requires matching body, fluid and memory initial
   conditions. Local chine submergence is only about 2.8 mm at equilibrium, so
   a small wave-height/hull-length ratio alone does not establish linearity;
   record wetting changes and representation uncertainty during amplitude checks.
3. Implement and verify the minimum missing Topwater time-domain physics, with
   early measured cost screens. Compare actual coupled histories, resolve physical
   discrepancies and uncertainties, and obtain the first passing practical case.
4. Promote only after independent review, reproducible evidence, blocking
   regressions, commits/pushes on baseline/fluid-body and tracker updates. At each
   bounded checkpoint retain failures, confirm every exact job handle terminal,
   and write a fresh-context handoff. Failed investigations leave this goal open.

No old failed study is retroactively qualified. Any prospective reference-policy
revision must be explicit and reviewed; study 003's original gates/results remain
unchanged. The former 60-second/20x targets remain later requirements, not claims
established by this first one-second/200x milestone.

## Current controlled case

The first prospective case is frozen in validation/mild-periodic-case.json:
a controlled three-component directional sea, not a buoy realization. Its
periodic body-motion candidate and spatial/sampling assessment are documented
in mild-periodic-protocol.md. This explicitly periodic benchmark avoids an
unqualified broadband conversion, but it still requires matching fluid initial
fields and nonlinear compatibility/representation checks. Body traces alone
do not qualify the reference or Topwater.
