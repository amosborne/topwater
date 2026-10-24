# Reference study 003: declared, not yet run

## Scope and exclusions

Execute the complete canonical Capytaine convergence protocol using the reviewed
fixed-seed transfinite family. This is the next bounded study after the geometry
checkpoint in fixed-topology-preflight.md. The overarching qualified-reference
goal remains active and incomplete. Do not launch hydrodynamics until the geometry
checkpoint's fresh freeze, evidence and handoff are committed and pushed.

Preserve the canonical 24 ft boat, calm-water equilibrium, world CG, full dry
inertia and restoring matrix, density 1025 kg/m3, gravity 9.81 m/s2, infinite
depth, zero speed, axes and exp(-i*omega*t) wave convention. No Topwater solver
changes/comparisons, new formulations, nonlinear qualification, optimization or
timing experiments. No change to accuracy gates or resource estimator.

## Inputs, cases and controls

Canonical contract: validation/canonical-boat.json. Fresh selected freeze:
artifacts/validation/frozen-fixed-001.json (also deliver its complete durable
bundle and proof summary at the geometry checkpoint). Verify current hashes
and physical/input identity before launch; regenerate on any producer change.

The fixed graded seed has 308 panels. Logical levels 1/2/3/4 map to affine
subdivisions 2/3/4/5 and 1232/2772/4928/7700 actual panels. All use the same
seed geometry/grading map and farthest-first vertex convention. Levels are not
all nested. The geometry-only evidence does not establish hydrodynamic accuracy.

Execute all 132 physics calls: all 36 A-infinity entries; full six-mode radiation
A/B at 0.10/0.25/0.50 Hz; diffraction, Froude-Krylov, total excitation and coupled
complex six-DOF RAOs for every combination with 0/45/90 degree propagation.
Use Capytaine 3.0.0, unchanged Green settings and one native thread. Archive
actual producers and inputs before the job, then keep them immutable until the
exact process handle is terminal. No other numerical job runs concurrently.

The conservative workspace estimate is 11,703,102,464 bytes on the finest mesh,
below the unchanged 12 GiB estimate cap. Recheck current available OS memory
before launch; require at least the estimate plus 2 GiB available headroom.
Retain the timestamped observation. The 1800-second elapsed ceiling remains
checked between calls and after levels; a native call is not forcibly timed out.

## Acceptance and deliverables

For every normalized real/complex quantity, both finest increments must fit
0.001*abs(Q4)+1e-6 and the latest increment cannot exceed the previous plus1e-6.
Retain empirical U=max(abs(Q4-Q3),abs(Q3-Q2)) entrywise. No magnitude-only complex
comparison, omitted case, Richardson extrapolation or certified bound claim.

Deliver every raw matrix and response, exact frozen mesh/input/producers,
commit/source/dependency identities, increments and uncertainty, explicit
qualified/unqualified decision, independent adversarial review and reproduction
instructions. Make evidence self-contained and durable with bundle SHA256.
Missing/stale/incomplete/failed evidence cannot pass. This decision concerns
only the frozen-surface linear external reference; Topwater remains unqualified.

## Exact next command and stopping checkpoint

After inspecting the completed geometry checkpoint and rechecking resources:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation reference --frozen artifacts/validation/frozen-fixed-001.json --output artifacts/validation/reference-study-003.json
```

Use fresh output; do not overwrite an existing report, NPZ or source archive.
Stop when the complete declared study terminates or hits a declared validity or
resource failure. Retain partial arrays on failure and confirm the exact handle
is terminal before edits or another job. Diagnose every failed family/entry;
do not repeat a failed study unchanged, omit a finest level or relax a gate.
Resolve review findings, commit/push baseline/fluid-body, update issue #30 and
write a fresh-context handoff. A failed investigation leaves the overarching
qualification goal active; only complete passing evidence and final review
permit reference promotion.
