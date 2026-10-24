# Reference study 003: completed, unqualified

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
0.001*abs(Q4)+1e-6 and the latest increment cannot exceed the previous plus 1e-6.
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

## Outcome

All four meshes and all 132 physics calls completed. Exact producer commit:
2de87a8 (full commit and source hashes in the report). The original entrywise
coefficient policy failed in every response family; this study is unqualified.
No producer or physical input changed during the run. Session 73213 is terminal,
exit 1; the process observer 45767 is terminal, exit 0.

[Complete raw/frozen/producer/resource bundle](../validation/evidence/reference-study-003.zip)
and [summary/SHA256](../validation/evidence/reference-study-003-summary.json).
The bundle includes the complete bound freeze. The memory observer held exact
process 25868 through exit; observed peak working set was 3,019,403,264 bytes and
peak committed memory 4,453,543,936 bytes. These are resource observations, not
Topwater runtime measurements. The elapsed reference-generation time 448.064 s
does not measure an own-model simulation; no coupled Topwater runtime exists yet.
The complete freeze is the nested bound-freeze.zip; its SHA256 is explicitly
bound in the study summary and was independently verified.

| Family | Increment failures | Trend failures | Entries |
|---|---:|---:|---:|
| A-infinity |18|0|36|
| finite A |44|2|108|
| B/omega |31|1|108|
| diffraction |13|0|54|
| Froude-Krylov |2|0|54|
| total excitation |9|0|54|
| complex RAO |21|1|54|

All six A-infinity diagonal latest increments also fail (about 1.93 to 5.81
times budget). The failure is not solely a relative-error effect in small
couplings. Independent review recomputed coverage, array/mesh/source identity,
excitation addition, full-mass/restoring impedance signs, every gate and the
uncertainty arrays, and found no packaging or provenance blocker. The explicit
decision remains unqualified under the original policy.

An exploratory scaling of the existing nine steady single-wave RAOs at the
declared 0.01 m incident amplitude gives last-two-increment envelopes of at most
1.2883e-5 m translation and 0.001703 degree rotation. This is not a continuum
bound, transient/startup comparison, new mild-sea result or qualification.
It motivates assessing the user's actual motion-error budget before another
blanket matrix-refinement study. The original failed gates remain unchanged.

## User-directed next outcome

The active goal now includes Topwater's own small-motion nonlinear six-DOF
TIME-DOMAIN simulation in a frozen mild sea. It must run for at least 1 s, match
Capytaine-derived time histories within 0.1 degree / 1 cm in every axis including
uncertainty, and be practical at no more than 200x real time including normal
Topwater preparation. Validation comes before major speed improvements; reject
obviously impractical routes using early cost screens. See
time-domain-baseline-goal.md for scope, acceptance, reuse and stopping checkpoints.

Independent review approved this prospective ordering, with requirements to
freeze the mild case and full initialization first, establish all time-reference
uncertainty controls, preserve the original failed study, and keep early cost
screens separate from major optimization. The historical adapter is radiation
math only; periodic RAO synthesis and rest-start transients require different
initial-state treatment. These findings are incorporated in the current goal.
