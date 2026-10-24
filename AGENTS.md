# Validation before experimentation

## Current user direction: Capytaine plus WEC-Sim only

The user replaced the previous goal. Establish an end-to-end canonical
Capytaine -> BEMIO -> WEC-Sim reference: a10s six-DOF simulation, isometric video
and characterization of motion variation with solver mesh resolution. Start
the boat at rest in locally calm water; waves must roll in from upstream, not
appear under the hull. No global load/amplitude step or imposed ramp is the
accepted startup. Use the same traveling input for forces and visualization.
No Topwater-specific solver code, implementation, comparison or timing work.

Begin with the minimum hull geometry as input; refine solver panels on that
same surface as needed. Obtain the first working time history before pursuing
refinement. Gate the SIX time-domain motions, including startup, and characterize
their numerical variation. Intermediate added mass/damping/coupling entries
are diagnostics, not independent qualification gates for this workflow. Keep
all original failed studies and their original policies/results unchanged.
No more blanket coefficient/quadrature campaigns. Missing/corrupt evidence or
invalid solves cannot establish a result, but old0.1% coefficient thresholds
must not block this new motion-based workflow.

This direction supersedes conflicting historical scope/order below. Read
docs/capytaine-wecsim-workflow.md and docs/checkpoint-handoff.md first. The
previous Topwater goal and coefficient-first qualification path are held.

The active application fixture is the 24 ft boat in `tests/reference_boat.py`.
The former canoe case is retired. Preserve old experiments in Git history;
do not use their passing checks or timings to qualify this boat.

The single implementation under qualification is the P0 Rankine-panel core
extracted from `4e8f511`, through `topwater.inertia.instantaneous_added_mass` and
its `bem`/`boundary` owners. Status: **unqualified**; hydrostatics and instantaneous
added mass exist, coupled time-domain wave response does not.
An unqualified conforming fluid boundary and mixed pressure/body component now
exists in `topwater.fluid_boundary`/`topwater.pressure`; boundary derivative
recovery and a complete nonlinear RHS remain missing. Its coarse-mesh physical
comparisons fail; see `docs/mixed-pressure-001.md`.
The paired body-only uniform family also fails the unchanged gates; current
reference diagnoses are recorded in `docs/paired-body-study-001.md` and
`docs/frozen-source-quadrature-001.md`. Frozen-source pressure integration is
material; frozen-source Gauss12/16/20 now passes conditional output integration
(docs/frozen-source-gauss-continuation-001.md). This does not qualify source
discretization; all earlier native-mesh and finite-wave failures remain failed.
The common-output54/216/864 source diagnostic confirms all three integration
controls pass while source increments fail; see docs/common-output-sources-001.md.
Do not spend the long one-level paired-family continuation seeking qualification;
its retained prior increment is already far outside the unchanged policy. The
consumer resource rehearsal and fixed-seed003 asset/physical/quality preflight
are complete. Fresh six-mode Ainf export/native reproduction and conditional
output integration now pass on its1232-panel mesh; source-mesh/full reference
remain unqualified. All four source bases now reproduce native matrices;
18 remaining source solves completed in66.45s under declared limits. Next is
a fixed2772 output continuation: cold12/16/20 fields completed but three
couplings fail both unchanged integration gates. Resolve output error before
larger source jobs; no automatic continuation or qualification claim.
See docs/fixed-field-2772-001.md and docs/checkpoint-handoff.md.
Existing finite source data do not contain
Ainf densities or all9 frequency-heading diffraction cases.
Source/reference qualification remains false.
The archived P1
diagnostic-tank runtime is a different implementation. Record current commit
and source hashes in each manifest; the word "baseline" alone is insufficient.

The current user-directed goal is a validated, practical, time-domain-only
Topwater baseline on the canonical boat in a frozen mild sea: at least one
second, every rotation within 0.1 degree and translation within 1 cm of a matched
Capytaine-derived reference, including uncertainty, and at most 200x real time.
See docs/time-domain-baseline-goal.md. Necessary completion of the selected
time-domain wave/fluid-body path and bounded cost screens are in scope. Major
optimization/alternative-formulation campaigns remain held until an actual
coupled comparison passes; do not spend indefinitely on a 5000x route. General
nonlinear/full-minute qualification follows this case. The active tracker is
GitHub issue #30.
The user's reset instructions are preserved in `docs/validation-reset.md`.

Every goal must declare scope, exclusions, deliverables, acceptance criteria
and a stopping checkpoint before work begins. The foundation reset is complete;
reference qualification is now a prerequisite of the broader time-domain goal.
The user's latest outcome and validation-before-major-optimization order take
precedence over older scope/timing holds. At each
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
5. Establish the small-motion nonlinear, time-domain, one-second/200x case first.
   Measure cost early and permit bounded feasibility fixes, but major speed work
   and general nonlinear/full-minute qualification follow a passing physical
   comparison. Capytaine is external reference data, never the runtime solver.

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
- Automatic CI runs only for pull requests targeting `main`. Continue local
  testing on development branches. Never merge or squash-merge a PR unless all
  required CI checks passed for its current head; missing, pending, skipped or
  failed required checks do not permit a merge.
- Keep one active implementation branch, `baseline/fluid-body`. Classify other
  branches as supporting or historical. Preserve commits before retiring refs;
  do not delete active worktrees or discard uncommitted work.
- Preserve generic mathematical and physics coverage when retiring an obsolete
  application. Do not rename the canoe or silently delete useful checks.
- Issue #30 owns promotion requirements and a concise optimization table.
  Each avenue has one status: open, in-progress or closed, with rationale and
  evidence. Held work is open, not authorization to execute it.
