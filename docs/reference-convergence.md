# Capytaine reference qualification goal

Objective: establish a qualified, reproducible Capytaine reference for the
canonical 24 ft boat, resolving convergence failures until the declared
acceptance gates pass. The completed foundation reset is the starting point.

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

## Correction checkpoint after study 001

The original mesh family failed all response families. A reviewed contract
revision declares graded Gmsh 4.15.2 meshes on the same frozen physical surface.
All original physical arrays and input JSON are unchanged; the overall identity
changes because it also binds the contract and derived reference meshes.
Geometry preflight gives 2464/3538/5152/7604 panels. The independent mesh levels
are ordered refinement, not nested subdivision. No extrapolation is used.
Resource revision, geometry proofs and review resolutions are recorded in the
study report. The qualification goal remains active and incomplete.
