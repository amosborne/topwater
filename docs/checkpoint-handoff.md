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
- Next action: create fresh freeze under committed current contract at
  artifacts/validation/frozen-graded-002.json. Recheck available memory and verify
  all original physical arrays/input JSON match reset. Then serially run reference
  with that freeze to artifacts/validation/reference-study-002.json. Never overwrite
  an output/archive or reuse stale freeze/producer identities. Use workspace
  uv cache: uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation.
- Stop/checkpoint after the full declared study or resource/validity failure;
  preserve every raw matrix, diagnose failed entries and declare/review the next
  bounded correction. The qualification goal remains active until all gates and
  final review pass. Topwater comparisons/waves, nonlinear work and timing excluded.
