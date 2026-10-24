# Checkpoint handoff

- **Active goal/status:** establish a qualified, reproducible Capytaine reference
  for the canonical 24 ft boat, resolving convergence failures until unchanged
  acceptance gates pass. ACTIVE AND INCOMPLETE. Both hydrodynamic reference and
  Topwater remain unqualified. The fixed-topology GEOMETRY checkpoint is complete.
- **Branch/identity:** baseline/fluid-body. Numerical/source checkpoint for the
  selected family and freeze: 96925446db194ed4f46e45ae74a3a62cf2d27525. Inspect live
  git status/log; subsequent commits may contain documentation/evidence only.
  Every owned manifest has actual commit and source hashes. No Topwater runtime
  numerical source was changed. No PR was created and main was not pushed.
- **Prior failures retained:** studies 001/002 each ran all four meshes/all 132
  calls and failed the declared convergence gates. Vertex-order diagnostic 001
  completed 48 calls and found material sensitivity; the reviewed farthest-first
  convention is selected. Complete bundles and summaries are under
  validation/evidence. Read docs/reference-study-002.md for these diagnoses.
- **Completed geometry work:** a centroid-quad fan was rejected without a job
  because it would preserve poor shapes in the real thin chine strip. Candidate 001
  generated 382 seed panels, then correctly failed resource preflight for its
  declared finest 9,550 panels. Candidate 002 uses broad-facet seed sizes 1.5 m / 0.5 m,
  growth 0.3 and a facet-specific size target 20 times minimum planar width.
  It has 308 seed panels and fixed affine Gmsh transfinite subdivisions 2/3/4/5,
  producing 1,232 / 2,772 / 4,928 / 7,700 actual panels. Grading/patch topology is fixed;
  levels are not all nested. No hydrodynamic solve ran in this checkpoint.
- **Evidence:** docs/fixed-topology-preflight.md records decisions/review.
  validation/evidence/transfinite-preflight-{001,002}.zip retain actual seeds,
  meshes, producers, protocols and complete bound freezes. Companion summaries
  give SHA256, geometry and per-facet shape statistics. All patch connectivity,
  affine lattice, positive Jacobian, shared-curve and full seven-facet surface
  checks passed for candidate 002. Genuine acute-apex slivers remain; small area
  is not an accuracy claim. The next hydrodynamic study must measure their effect.
- **Selected canonical freeze:** artifacts/validation/frozen-fixed-001.json
  (complete durable bundle and proof at validation/evidence/frozen-fixed-001.*).
  Current source/archives and physical identity verified. Seed and every mesh are
  byte-exact to candidate 002; all original physical arrays/dtypes/input JSON
  match reset. Capytaine's actual face radii equal maximum vertex distances on
  every selected panel. The proof and verifier are in the bundle; no pressure
  or radiation solve is involved. Regenerate fresh on any producer/input change.
- **Controls/gates:** same canonical physical problem/cases/normalization and
  full dry inertia/restoring. Capytaine 3.0.0, Gmsh 4.15.2, one thread, unchanged
  1800s between-call budget, unchanged 12 GiB estimate cap/conservative estimator.
  Finest estimate 11,703,102,464 bytes. Reference acceptance stays 0.1%+1e-6
  normalized for every real/complex entry, both finest increments plus trend,
  empirical envelope uncertainty; no Richardson/certified-bound claim.
- **Regression/review:** 51 local gate/geometry tests passed. Missing selected
  meshes/seed, incorrect seed_count*N? mapping and stale transfinite producers
  fail closed. Independent fixed_topology_review recomputed all patch/surface
  proofs, audited archives/identities/physical bytes, and approved geometry/control
  adoption at 9692544 with no blocking findings. This is not reference qualification.
- **Live jobs:** none. Preflight 001 returned terminal exit1 (resource failure),
  preflight 002 terminal exit0 (geometry pass), canonical freezes and adoption proof
  terminal exit0. No other numerical process from this checkpoint remains live.
- **Next bounded action:** execute COMPLETE reference study 003, declared in
  docs/reference-study-003.md. Recheck current repository/producer/input identity,
  ensure fresh output, and record current available OS memory >= finest estimate
  plus 2 GiB (13,850,586,112 bytes). Archive actual producers/input before launch and
  keep immutable; run serially. Exact command:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation reference --frozen artifacts/validation/frozen-fixed-001.json --output artifacts/validation/reference-study-003.json
```

- **Stopping condition:** complete all 132 calls/all 36 Ainf entries, full finite
  A/B and every frequency/heading complex excitation/RAO, or stop on a declared
  resource/validity failure. Preserve raw data/increments/uncertainty and confirm
  exact process handle terminal before edits or another job. Resolve independent
  review and publish qualified/unqualified decision; commit/push, update issue #30
  and hand off. A failed study leaves the overarching goal active. Do not repeat
  unchanged failures, drop cases/finest meshes or loosen gates. Topwater comparison,
  wave implementation, nonlinear qualification, new formulations and timing held.
- **CI:** automatic CI only on PRs targeting main. Continue local tests. Never
  merge/squash-merge unless every required current-head CI check passed.
