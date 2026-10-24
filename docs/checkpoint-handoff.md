# Checkpoint handoff

- Active goal: establish a qualified, reproducible Capytaine reference for the
  canonical 24 ft boat, resolving convergence failures until unchanged gates
  pass. ACTIVE AND INCOMPLETE. Topwater also remains unqualified. The foundation
  reset is complete; see docs/reference-convergence.md for scope and exclusions.
- Branch: baseline/fluid-body. Inspect current git status/log before continuing.
  Every freeze/study manifest owns its exact commit and source hashes. No
  Topwater runtime numerical source was changed in this investigation.
- Completed evidence: study001 (original midpoint family) and study002 (reviewed
  graded planar family) each completed all four levels and all 132 calls, and
  failed every convergence family. Durable bundles and summaries with SHA256:
  validation/evidence/reference-study-{001,002}.{zip,-summary.json}. Each bundle
  includes raw matrices, frozen inputs, actual input JSON and producer archives.
- Study002 meshes: 2464/3538/5152/7604 panels, all seven original physical facets
  retained; original eight physical arrays/input JSON exactly match reset.
  Independently remeshed levels are ordered refinement, not nested subdivision.
  Latest A-infinity increments reach 6.27 times budget; FK has three trend
  failures despite all its increments fitting budget. Read docs/reference-study-002.md.
- Completed diagnostic: all 48 calls in vertex-order-diagnostic-001. Equivalent
  cyclic ordering changes responses by up to 280.27/76.52/16.59 reference budgets
  for Ainf/A/B. Installed Capytaine uses only first-vertex distance as face radius.
  Durable diagnostic bundle/summary under validation/evidence; bound freeze in
  study002 bundle. Diagnostic is unqualified and cannot substitute for convergence.
- Reviewed correction: reference_mesh.farthest_first preserves physical triangles
  and winding while putting the farthest vertex first. Contract binds this choice;
  actual ordered triangles are frozen. Independent reference_adversarial_review
  approved adoption. This fixes representation sensitivity, not spatial convergence
  or native kernel approximation error. All 33 gate/geometry tests passed.
- Current contract: same physical problem and cases, unchanged 0.1% + 1e-6
  normalized budget, both finest increments plus trend, empirical envelope;
  no Richardson/continuum-bound claim. Explicit 12 GiB estimate cap; unchanged
  12-dense-array estimator, 1800 seconds between physics calls. Do not use
  partial memory observations from study002 as full final peak measurements.
- Fresh freeze: old frozen-graded-002 is stale after the reviewed ordering/control
  revision. Generate a fresh freeze before any subsequent current-source work;
  preserve all previous archives. Recheck original physical arrays/input JSON.
- Live jobs: none. Sessions 11777, 11118 and 52130 all terminal with nonzero
  exits for failed convergence/representation gates. No numerical job is pending.
- CI policy: automatic CI only on PRs targeting main. Continue local testing.
  Never merge/squash-merge without all required current-head CI passing. No PR
  was created; development checkpoints pushed only to baseline/fluid-body.
- Exact next bounded action: design and geometry-preflight a FIXED-TOPOLOGY GRADED
  reference mesh family, with NO hydrodynamic solves in that checkpoint.
  Decompose the same convex planar facets into fixed triangular/quad patches;
  use pinned Gmsh transfinite surfaces with shared compatible edge counts.
  Keep patch topology, diagonal pattern and continuous grading map fixed; sample
  it more finely over four nontrivial levels (target current 2500-7600-panel
  range, but report actual counts). Derive per-segment ratios from a fixed grading
  exponent; do not keep ratios fixed while changing counts. Disable independent
  unstructured remeshing/optimization, reflect the verified half, apply farthest
  first and require all existing geometry checks. Record patch definitions, curve
  distributions, quality statistics, physical equality and resource preflight.
- Stop at that reviewed geometry manifest; commit/push, update issue #30 and
  refresh this handoff. Then declare the next complete convergence study. The
  overarching qualification goal stays active until all gates and final review
  pass. No Topwater comparisons, wave implementation, new formulations, nonlinear
  work, optimization experiments or performance claims are authorized here.
