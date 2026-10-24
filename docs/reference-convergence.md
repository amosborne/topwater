# Capytaine reference qualification goal

Historical coefficient-first goal, REPLACED by the user's Capytaine/BEMIO/
WEC-Sim10s video and motion-convergence goal. Read capytaine-wecsim-workflow.md.
Preserve every old failure under its original rules; do not continue this old
campaign or apply its intermediate coefficient thresholds as new motion gates.

Objective: establish a qualified, reproducible Capytaine reference for the
canonical 24 ft boat, resolving convergence failures until the declared
acceptance gates pass. The completed foundation reset is the starting point.
This reference work is a prerequisite of the broader user-directed
[one-second time-domain baseline](time-domain-baseline-goal.md), not its endpoint.

## Scope and exclusions

Use the frozen calm-water physical surface, pose, CG, full dry inertia,
restoring matrix and conventions in validation/canonical-boat.json. Execute the
existing four-level study first. Diagnose failures and correct the reference
harness or numerical controls with evidence. Further refinement/resource
revisions require an explicit rationale and review before launch. Physical
surface changes, relaxed accuracy gates, Topwater solver changes/qualification,
new formulations, nonlinear validation, optimization and timing are excluded.

## Deliverables and acceptance

Retain immutable producer/input archives and raw results for every study:
all 36 A-infinity entries, full finite A/B, diffraction, Froude-Krylov, total
excitation and complex coupled six-DOF RAOs for the full declared case product.
Record commit, source hashes, dependency identities, mesh increments and
entrywise empirical uncertainty. Missing/stale/incomplete evidence fails closed.

The unchanged contract requires both finest increments within 0.001 relative
plus 1e-6 normalized absolute for every quantity, with the latest increment no
larger than the previous plus the absolute floor. Uncertainty is the envelope
of those increments; it is not a certified continuum bound. No unsupported
Richardson extrapolation. Independent adversarial review must resolve physical
scope, signs, units, uncertainty, provenance and gate coverage before promotion.

Deliver a reproducible accepted reference and blocking reference checks,
resolved review, updated issue #30 and fresh-context handoff. This qualifies
only the frozen-surface linear reference, not Topwater or nonlinear motion.

## Checkpoints and stopping condition

1. Audit readiness and producer/input ownership, then launch the existing study.
2. Preserve its complete or failed outcome and diagnose individual failures.
3. Declare and review each bounded correction, then run the next serial study.
4. Promote only after all unchanged accuracy gates and independent review pass.

Numerical jobs run serially; confirm the exact handle terminal before edits or
another job. At checkpoints commit/push baseline/fluid-body, update issue #30
and hand off exact state. Failed studies do not complete the goal. Resource or
numerical limits produce an incomplete checkpoint with the next required
decision; never silently omit a level or loosen acceptance.

## Readiness record

Initial clean worktree: baseline/fluid-body at
4d9428a (full commit recorded in the first study manifest). Existing freeze is
artifacts/validation-reset/frozen.json; the reference command verifies current
physical/validator identities before solving. Declared meshes: 68, 272, 1088,
4352 panels; estimated finest workspace 3,846,176,768 bytes, 8 GiB estimate cap,
1800 seconds checked between physics calls. This describes reset readiness;
the subsequent [study 001](reference-study-001.md) is complete and unqualified.

## Historical correction checkpoint after study 001

The original mesh family failed all response families. A reviewed contract
revision declares graded Gmsh 4.15.2 meshes on the same frozen physical surface.
All original physical arrays and input JSON are unchanged; the overall identity
changes because it also binds the contract and derived reference meshes.
Geometry preflight gives 2464/3538/5152/7604 panels. The independent mesh levels
are ordered refinement, not nested subdivision. No extrapolation is used.
Resource revision, geometry proofs and review resolutions are recorded in the
study report. The qualification goal remains active and incomplete.

## Current geometry checkpoint

Study 002 completed unqualified, exposing both spatial/trend failures and a
separately diagnosed vertex-order sensitivity. The reviewed farthest-first
control remains selected. The [fixed-topology preflight](fixed-topology-preflight.md)
then qualified the geometry of a single graded seed and its affine transfinite
refinements, retaining every physical facet and original input exactly.

The canonical contract now selects 308 seed panels with N=2/3/4/5, giving
1232/2772/4928/7700 panels. The full physical reference remains unqualified.
Stop this geometry checkpoint after the current canonical freeze and evidence
are reviewed, committed/pushed and handed off. The next bounded action is the
complete 132-call reference study on those frozen meshes, after current-state
identity/resource checks. No accuracy or resource gate is relaxed.

## Study 003 outcome and updated user target

Study 003 completed all 132 calls and failed its original coefficient/trend
criteria. It remains unqualified; see reference-study-003.md and the preserved
raw evidence. The user subsequently clarified the desired outcome: an actual
small-motion nonlinear Topwater run in a mild sea, at least one second, within
0.1 degree/1 cm of matched Capytaine-based time histories, and practical at
200x real time. Validation precedes major speed improvements; early cost
screens and bounded feasibility fixes are permitted.

The next reference decision must address the selected case's time-history
uncertainty and matched initialization before more blanket matrix refinement.
No failed study is retroactively qualified, and no current gate is silently
relaxed. Declare and review any prospective output-budget protocol first.
