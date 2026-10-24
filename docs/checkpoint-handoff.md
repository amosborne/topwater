# Checkpoint handoff

- Current goal: qualified, reproducible Capytaine reference for the canonical
  24 ft boat, resolving convergence failures until unchanged gates pass.
  Active and incomplete; Topwater remains unqualified. See reference-convergence.md.
- Active branch: baseline/fluid-body. Inspect live git status/log for the exact
  current source checkpoint; every freeze/study manifest also records its commit
  and complete producer hashes. No runtime solver sources changed.
- Completed study: reference-study-001 ran all four midpoint meshes/all132 calls,
  complete but unqualified in every convergence family. Durable raw results,
  bound freeze/input and producers: validation/evidence/reference-study-001.zip;
  summary and SHA256 alongside it. See docs/reference-study-001.md.
- Correction: validation-only Gmsh4.15.2 graded planar facets on identical physical
  surface, verified y=0 reflection. All thin chine strips retained. Original
  eight physical arrays and input JSON match reset exactly. Mesh levels are
  independently remeshed ordered refinement, not nested subdivision.
- Geometry preflight: 2464/3538/5152/7604 panels; all surface checks passed. Initial
  frozen-graded-001 failed the original8GiB cap before any hydrodynamic solves.
  It remains locally preserved; summary in validation/evidence/graded-mesh-preflight-001.json.
- Declared revision:12GiB estimate cap, same conservative estimator, finest
  11,417,873,408bytes;1800s between calls unchanged. All accuracy gates unchanged:
  0.1%+1e-6 normalized, both finest increments and trend, componentwise envelope.
- Independent reviewer reference_adversarial_review approved readiness after
  single-thread initialization, initial-report persistence, geometry tolerance
  binding and missing-graded-mesh failclosed corrections.32 gate/geometry tests
  passed. No convergence acceptance or physical qualification claim.
- Live jobs: none at this source checkpoint. Confirm exact sessions terminal
  before changing producers or launching another numerical job.
- Study002 is now terminal, all132 calls complete, every family unqualified.
  Durable complete bundle: validation/evidence/reference-study-002.zip; summary
  alongside. Source6fae539, controls unchanged during run; later commits only
  change automatic CI/merge policy. See docs/reference-study-002.md.
- CI now triggers only PRs targeting main; local tests remain required. Never
  merge/squash-merge a PR without all required current-head CI passing.
- Next action: run bounded vertex-order diagnostic, reviewed and declared in
  docs/reference-study-002.md. Command: uv run --locked --no-sync --cache-dir
  artifacts/uv-cache --group correlation python -m tests.validation.diagnose_order
  --frozen artifacts/validation/frozen-graded-002.json
  --output artifacts/validation/vertex-order-diagnostic-001.json
  It changes only cyclic triangle ordering;48 radiation calls, no qualification.
  Existing freeze remains current because this new standalone diagnostic does
  not change its required physical/validation producers or contract. Never
  overwrite checkpoints. Keep all producers immutable while running.
- Stop/checkpoint after the full declared study or resource/validity failure;
  preserve every raw matrix, diagnose failed entries and declare/review the next
  bounded correction. The qualification goal remains active until all gates and
  final review pass. Topwater comparisons/waves, nonlinear work and timing excluded.
