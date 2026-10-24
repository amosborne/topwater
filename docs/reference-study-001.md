# Reference study 001: unqualified

All four declared midpoint meshes and all 132 physics calls completed. Every
convergence family failed the unchanged contract; no reference is qualified.
Elapsed study time was 155.752 seconds (resource accounting only, not a runtime
performance claim). Exact commit and source hashes are in the archived report.

The [summary](../validation/evidence/reference-study-001-summary.json) records
gate ratios and counts. The [complete evidence bundle](../validation/evidence/reference-study-001.zip)
contains reference JSON/NPZ, bound frozen JSON/NPZ and both complete producer
archives including input JSON. Extract it to a separate directory; do not run
historical producers in the current worktree. The old gate accepts neither
these failed results nor changed producers as a qualified reference.

| Family | Entries outside increment budget | Total entries |
|---|---:|---:|
| A-infinity | 36 | 36 |
| finite A | 106 | 108 |
| B/omega | 96 | 108 |
| diffraction | 52 | 54 |
| Froude-Krylov | 36 | 54 |
| total excitation | 50 | 54 |
| complex RAO | 46 | 54 |

These counts concern increment budgets; trend is a separate required check.
The finest heave A-infinity changed by about 1.7%, so merely imposing mesh
symmetry cannot resolve the diagonal convergence failure. Asymmetric original
diagonals also generate large odd couplings that decay with refinement.

## Independent review and correction checkpoint

Independent reviewer `reference_adversarial_review` found no physical sign/unit
defect in the six-mode mapping, CG, full dry inertia, restoring sign, impedance,
normalizations or complex increment policy. Findings to resolve before promotion:

- Native Green tabulation occurred outside the single-thread limit. The generator
  now loads the native library before discovery, limits initialization and solves
  together, and records/checks effective pools.
- Initial evidence was not saved before expensive initialization. It is now saved
  before native loading and again before tabulation.
- Local ignored artifacts were insufficient durable delivery. This complete
  failed-study bundle is committed; accepted evidence must likewise be bundled.

The control/gate suite passes 26 tests with a workspace-local temporary directory.
The initial default-temp invocation had 15 sandbox directory-creation errors;
those are not counted as passes.

## Next bounded correction

Reconstruct the connected coplanar facets of the same frozen wet polyhedron and
use pinned Gmsh 4.15.2 in the validation-only dependency group. Preserve all thin
chine strips, physical edges, waterline and outward orientation. Verify actual
reflection symmetry before half-surface meshing and reflection; still compute
and gate every response entry. Use explicit balanced/edge-graded mesh controls,
with geometry evidence, actual panel counts and preflight before physics launch.

Geometry acceptance is roundoff-level surface preservation: per-facet area,
oriented area and first/second area moments agree within 1e-10 relative plus
1e-12 SI absolute; every panel is contained in its original planar facet with
coordinate tolerance 1e-11 m. Require matching feature edges and waterline,
no gaps/overlaps or added lid. Prove physical arrays and input JSON unchanged
when regenerating the contract-bound freeze. No hydrodynamic accuracy tolerance
is relaxed. Independent review endorsed this direction subject to these proofs
and explicit contract revision before another hydrodynamic study.

## Graded-mesh readiness decision

[Geometry preflight evidence](../validation/evidence/graded-mesh-preflight-001.json)
records 2464/3538/5152/7604 panels and exact equality of all eight original
physical arrays and input JSON bytes. This first freeze correctly failed the
8 GiB resource preflight; no hydrodynamic calls were made. Its producers and
arrays are retained locally as artifacts/validation/frozen-graded-001.*.

The contract explicitly revises the estimate cap to 12 GiB (12,884,901,888 bytes)
for the next bounded study. The same conservative estimator requires
11,417,873,408 bytes for 7604 panels. Windows GlobalMemoryStatusEx observed
33,832,407,040 total and 15,873,208,320 available physical bytes on
2026-10-23 during geometry readiness (the observation was not timestamped). Recheck available memory before
launch. This is a workspace estimate, not an enforced OS cap. Elapsed limit
remains 1800 seconds; all accuracy gates remain unchanged.

Reviewer `reference_adversarial_review` accepted this resource revision and the
surface proof. Additional review findings resolved: the contract's geometry
tolerances are checked against the implementation constants; missing graded
mesh arrays cannot fall back to midpoint meshes; actual panel counts must
increase; and diagnostics explicitly call the new family independently
remeshed, not nested. The gate/geometry suite passes 32 tests. Accepted evidence
must still pass the complete hydrodynamic convergence study and final review.
