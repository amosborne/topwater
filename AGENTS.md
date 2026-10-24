# Validation before experimentation

The active application fixture is the 24 ft boat in `tests/reference_boat.py`.
The former canoe case is retired. Preserve old experiments in Git history;
do not use their passing checks or timings to qualify this boat.

The single implementation under qualification is the P0 Rankine-panel core
extracted from `4e8f511`, through `topwater.inertia.instantaneous_added_mass` and
its `bem`/`boundary` owners. Status: **unqualified**; hydrostatics and instantaneous
added mass exist, coupled time-domain wave response does not. The archived P1
diagnostic-tank runtime is a different implementation. Record current commit
and source hashes in each manifest; the word "baseline" alone is insufficient.

New formulations, optimization experiments, nonlinear 6DOF qualification and
timing work are on hold until a validated fluid–body baseline exists. Work on
reference convergence, regression gates and corrections exposed by those gates
is the current priority. The active tracker is GitHub issue #30.
The user's reset instructions are preserved in `docs/validation-reset.md`.

Every goal must declare scope, exclusions, deliverables, acceptance criteria
and a stopping checkpoint before work begins. The current goal is qualified,
reproducible Capytaine reference convergence (docs/reference-convergence.md);
the foundation reset is complete. At each
checkpoint commit and push, update the tracker, resolve adversarial findings
and write a fresh-context handoff naming the next bounded validation goal.
A completed investigation may reject a candidate; it does not qualify it.

## Required order

1. Freeze the canonical calm-water physical problem: hull, equilibrium pose,
   world CG, full dry inertia, density, gravity, depth, axes and wave convention.
2. Establish a Capytaine mesh-convergence study on that same physical surface.
   Report all 36 added-mass entries and uncertainty for each response of
   interest. Two-level agreement and same-mesh correlation are not convergence.
   Do not apply Richardson extrapolation without evidence for its assumptions.
3. Validate the existing Topwater path against the qualified reference: all six
   added-mass responses, their couplings, and selected small-amplitude wave
   responses including complex phase. Separate reference uncertainty from
   Topwater spatial, domain and temporal error.
4. Make these comparisons blocking regressions after formulation changes.
   Missing dependencies, cases or reference data; skipped required comparisons;
   stale producer/input identity; and unresolved uncertainty cannot pass.
5. Only after those gates pass, return to nonlinear 6DOF validation and timing.

Record numerical tolerances and acceptance criteria explicitly before judging
a result. Residuals, symmetry, positive energy, fitted loads and local ODE
tolerances do not establish physical or full-minute motion accuracy. Never
silently relax a gate to accept a change.

## Development discipline

- Keep one simple runtime execution path. Supported boats change through data,
  not solver edits or reference-mesh assumptions.
- Use an adversarial review before promoting numerical changes. Review physical
  scope, signs, units, uncertainty, redundant code and unnecessary machinery.
- Run numerical jobs serially. Archive the actual producer graph and physical
  inputs before launch, keep them immutable while running, and confirm the
  exact process handle is terminal before changing them or starting another job.
- Preserve checkpoints and push development branches. Never push to `main`.
  Do not create PRs for bookkeeping.
- Keep one active implementation branch, `baseline/fluid-body`. Classify other
  branches as supporting or historical. Preserve commits before retiring refs;
  do not delete active worktrees or discard uncommitted work.
- Preserve generic mathematical and physics coverage when retiring an obsolete
  application. Do not rename the canoe or silently delete useful checks.
- Issue #30 owns promotion requirements and a concise optimization table.
  Each avenue has one status: open, in-progress or closed, with rationale and
  evidence. Held work is open, not authorization to execute it.
