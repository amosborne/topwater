# Validation foundation and stopping checkpoint

The implementation under qualification is the **P0 Rankine-panel fluid-body
core**, extracted from `4e8f511`, on branch `baseline/fluid-body`. Its numerical
entry point is `topwater.inertia.instantaneous_added_mass`, using the `bem`,
`boundary` and `rankine` owners. Current source hashes and commit belong in
each generated manifest. The extraction SHA is provenance, not current proof.

Status: **unqualified**. Hydrostatics and instantaneous added mass are
implemented. Coupled time-domain wave response is not. The archived P1
diagnostic-tank runtime is a different implementation and is not substituted.

## Active qualification goal

The foundation reset below is complete. The active goal is
[qualified, reproducible Capytaine reference convergence](reference-convergence.md),
including bounded diagnosis/correction checkpoints until the unchanged gates
pass. A failed individual study leaves this qualification goal incomplete.

## Completed foundation reset

Restore a controlled validation-first foundation. Deliver the canonical-only
fixture, one implementation identity, the reviewed contract and executable
validation stages, concise issue #30, preserved branch/worktree inventory and
a fresh-context handoff. Preserve generic mathematical and physics checks.

Completion means these reset deliverables are reviewed, checked, committed and
pushed. It does **not** mean a reference or fluid-body simulation is qualified.
No Capytaine convergence study, new formulation, optimization, nonlinear
qualification or timing experiment is part of this goal. Stop at the handoff.

## Frozen contract

[canonical-boat.json](../validation/canonical-boat.json) declares the physical
problem, selected cases, numerical/resource limits, normalization and gates.
The canonical reference is at calm-water equilibrium with zero twist. It is
distinct from the observed SoCal sea realization used by the future nonlinear
case. Full dry inertia and restoring couplings are retained.

Reference mesh policy: 0.1% relative plus 1e-6 normalized absolute, both finest
increments and a non-increasing trend within the absolute floor. Report the
componentwise envelope of those increments as empirical mesh uncertainty.
No Richardson extrapolation or certified continuum bound is claimed.

Candidate policy: 1% relative plus 1e-5 normalized absolute, applied to every
physical entry with reference uncertainty charged against the same budget.
These are predeclared **linear engineering verification criteria**, not a
derivation of the eventual 0.1 degree / 1 cm minute requirement.

## Executable sequence

Use fresh output paths: existing reports, arrays and producer archives must
never be overwritten. Numerical jobs run serially with immutable sources and
inputs. The resource timeout is checked between physics calls, not a promise
to interrupt a native solve mid-call.

```powershell
uv sync --locked --group correlation
uv run --locked python -m tests.validation freeze --output artifacts/validation/frozen.json
```

`freeze` records actual inputs, wetted panels, equilibrium pose, world CG, full
dry mass matrix, restoring matrix, current producer identities and the complete
study preflight. Its successful exit only validates the freeze operation.
The inspected reset snapshot is retained in
[`validation/frozen-physical-inputs.json`](../validation/frozen-physical-inputs.json).
The actual calm-water wet surface has 17 panels; the declared reference levels
have 68, 272, 1,088 and 4,352 panels. The finest conservative workspace estimate
is approximately 3.85 GB, within the declared 8 GiB estimate cap. This is a
preflight, not measured resource usage or a convergence result.

That reset mesh family was executed and rejected in [study 001](reference-study-001.md).
The active contract now uses Gmsh 4.15.2 graded planar-facet meshes, preserving
the exact physical surface with a verified geometric reflection. Its four
independently remeshed levels are an ordered refinement family, not nested
midpoint subdivision. Geometry preflight produced 2464/3538/5152/7604 panels.
The estimate cap is explicitly revised to 12 GiB; the estimator and accuracy
gates are unchanged. See the study report for review and resource evidence.

**Active goal: Capytaine reference qualification.** Regenerate a fresh freeze
under the current contract before reference generation. Execute:

```powershell
uv run --locked --group correlation python -m tests.validation reference --frozen artifacts/validation/frozen.json --output artifacts/validation/reference.json
```

This command must run every declared mesh/case, preserve raw matrices and
complex excitation/RAOs, and decide qualified/unqualified under the contract.
A failed study is a completed investigation with an unqualified reference;
identify the next bounded convergence problem rather than weakening the gate.

Only after an accepted reference, the subsequent linear-validation goal uses:

```powershell
uv run --locked python -m tests.validation compare-inertia --frozen artifacts/validation/frozen.json --reference artifacts/validation/reference.json --output artifacts/validation/inertia.json
uv run --locked python -m tests.validation compare-waves --output artifacts/validation/waves.json
uv run --locked python -m tests.validation qualify --frozen artifacts/validation/frozen.json --reference artifacts/validation/reference.json --inertia artifacts/validation/inertia.json --waves artifacts/validation/waves.json --output artifacts/validation/qualification.json
```

**The last two commands must currently fail.** Wave comparison reports an
unsupported capability; qualification cannot override that with a supplied
passing artifact. Missing/skipped/stale/unqualified required results fail too.
The time-domain domain/step/startup/window/amplitude controls are explicitly
unresolved. A later reviewed contract revision must resolve them in the selected
implementation; these commands are not a claim that a wave solver exists.

CI keeps generic tests/package checks separate from **Physical baseline
qualification**. The latter freezes the current inputs and runs the real
qualification command without generating references. Its expected current
result is red/unqualified. This is a blocking promotion policy and status
check; GitHub branch-protection enforcement is not asserted or changed here.

After actual linear-wave gates pass, use separate goals for nonlinear 6DOF
minute accuracy and then performance improvement. Capytaine remains external
validation data, never the runtime solver or a duplicate pressure load.

## Provenance and retained coverage

| Material | Origin | Active purpose |
|---|---|---|
| Physical core, canonical boat/sea, independent panel references | `4e8f511` | Existing implementation and verification, without redesign |
| Generic sphere mesh and closed-boundary tests | `8c7767e` | Mathematical BEM/kernel checks, not another application boat |
| Explicit Capytaine units/constants correction | `de81051` | Validation helper, integrated selectively |
| Generic hydrostatics controls | Restored from archived checks on the actual boat | Six load derivatives, altitude, closure, submerged volume and failures |
| Full optimization graph | Tag `archive/optimization-2026-10-23` | Historical evidence, held |

The canoe fixture, image, runners and dedicated application tests are absent
from this branch. No historical worktree or uncommitted work is discarded.
See [branch purpose and preservation](branches.md).
