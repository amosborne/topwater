# Checkpoint handoff

Updated 2026-10-25. The latest user message controls scope. Read
[AGENTS.md](../AGENTS.md), [coordination](experts/coordination.md) and the
[shared physical contract](experts/shared-contract.md); consult relevant expert
states rather than repeating historical setup.

## Current authorization

Implement the approved canonical example/workflow plan. Native coupled dynamics
remain excluded and must fail clearly. The first working backend is WEC-Sim.
One flat `examples/simple-boat/` holds `boat.npz`, `sea.npz`, `result.npz`,
`video.mp4` and the newly requested `sea.log`, `result.log`, `video.log`.
A future generate-boat command will use boat.log; it is not implemented now.

Public commands accept an optional directory, defaulting to the current directory:

```powershell
uv run topwater generate-sea [DIRECTORY] --record RECORD_ID [--seed N]
uv run topwater solve [DIRECTORY] DURATION_SECONDS [--solver=native|wecsim]
uv run topwater render [DIRECTORY]
```

Keep the commands independent. The user explicitly declined streaming/pipelining
solve into render; validation orchestration may run them sequentially through
completed result.npz. Current WEC integration is fixed at0.02s; surface output
uses uniform timestamps targeting0.08s (751 frames for the selected60s case).
Future native integration may be adaptive; output cadence is a separate choice.

The current case preserves the accepted24ft/9ft/2500kg open V-hull and stronger
measured46086/2023-12-23T00:40Z realization, seed42. Public generate-sea has now
reproduced EVERY stored sea array bit-exactly and created canonical sea.log.
No complete result.npz/video.mp4 has yet been accepted or published.

## Current implementation and evidence

Latest pushed Topwater checkpoint:5ef862a (recovery). New logging/triangulation
changes are reviewed and tested but awaiting the next source checkpoint.
Maintained WEC-Sim fork9952634f6f2e2746a1626484e27c3003d1e1d5a9 is pushed and
pinned; its two newest fixes close sidecar handles and print paths literally.

- Current full workflow/logging/geometry suite:228 passed, one optional libigl
  skip,14.98s; artifacts/workflow/logging-geometry-tests-attempt1.
- Real BEMIO lifecycle probe passes: original KH leak reproduced; two imports
  retain32 numeric arrays exactly each; malformed sidecars leak no handles;
  15 Windows fixture files rename/delete. bemio-handles-attempt3.
- Display-only MATLAB probe passes three Windows/percent/UNC paths with zero
  warnings; no wave/body setup or new engine. wave-display-check-attempt3.
- Public canonical sea generation/log replay passes with unchanged numeric
  input. canonical-sea-log-attempt1.

All these jobs are terminal, hashes/processes verified and leases released.
No numerical job is currently allocated. Execution ledger, not this snapshot,
owns current liveness. Do not re-poll terminated sessions79250/28263.

## Remaining numerical work

Actual canonical attempt2 terminated exit2 after1397.677s. Native integration,
startup checks and complete radiation/diffraction synthesis finished. Export
failed at frame277/time22.16s with one roundoff-degenerate triangle, area
3.04945e-20m2. This is a geometry construction defect, not demonstrated flooding.
Attempt1 had masked its primary error with a leaked KH.dat cleanup failure;
that independent importer/cleanup bug is fixed and verified.

The completed field and failing frame are retained privately at:

`C:/Users/alex/AppData/Local/Temp/topwater-wecsim-mt5icmiq/surface-recovery/`

completed-field.npz SHA256:
`66039a8f78cfb8f9dcdec5d565e2197e9a06b01361521b19e0d5d46279613549`.
The original numerical data and motion must remain immutable.

An export-only internal-diagonal flip repairs the near-collinear polygon while
preserving all vertices, contact edges, orientation, area and existing validation
tolerances. Shared triangulation/physical gates are unchanged. Independent
Capytaine review approves it. Saved frames276/277/278 pass in8.423s; frame277
inputs/vertices/boundary chains are bit-identical, with only two faces changed.
Evidence: surface-frame-verification-attempt1/frames/report.json.

Next, under a fresh root-owned serial lease:

1. Run the frozen helper artifacts/workflow/check_cached_surface_geometry.py
   over all751 cached frames, with a fresh output subdirectory. This uses geometry
   only, no body/BEM/wave solve and no wet-area diagnostic. Retain any failing
   frame and original error; do not weaken validation or change motion/forcing.
2. Complete the public solve/render workflow with the new logs; verify full
   geometry, total-minus-incident metrics, warning triage and all decoded video
   frames/stills. Include logs in the canonical case; no solver intermediates.
3. Track result/video with Git LFS, commit/push baseline/fluid-body, update own
   issue30 and refresh this handoff after acceptance. Goal remains active.

## Physical and file requirements

Public NPZ files contain numeric arrays only, no pickle, generic metadata or
solver marker. Boat is an uncapped convex physical hull with one inferred
opening. Internal caps carry no mass/pressure. Result contains a complete physical
hull snapshot, CG, absolute poses and mandatory packed variable-topology total
water geometry with independent surface times. Render reads only result.npz.
See [the contract](npz-format.md), [mesh review](mesh-tooling.md) and each
[command page](cli/solve.md).

WEC total surface is incident+diffraction+radiation from its mean-hull model and
already computed native modal velocities. Reconstruction/replay never reintegrates
or changes motion. Selected provisional policy is cutoff24rad/s, spacing at most
.025rad/s, full-history causal projection and zero direct terms. Exact exterior
trace continuation handles vacated mean-interior targets; no interior Green
function evaluation or nonlinear contact accuracy is claimed.

Independent [method review](experts/numerical_methods_researcher/surface-reconstruction-review-2026-10-25.md)
accepts that selected numerical approximation. Its material signed-lag-versus-
causal sensitivity is42.7mm maximum/12.5mm RMS at tested points, not a true error
bound. Tiny frequency-grid sensitivity does not establish physical accuracy.
No mesh/time convergence or full-complex near-field equivalence is claimed.

Fatal errors prevent complete supported output and preserve previous data;
latest-attempt logs are replaced even on failure. WEC fixed-mean restoring and
represented wet-set departures warn at the agreed10% engineering threshold.
The initial represented total-elevation gate is separate from full-hull incident
eta/velocity bounds; no total disturbed-fluid velocity bound is established.
See [logging](cli/logging.md) and [actual warning triage](cli/warning-triage.md).
Expected warnings are not automatically harmless or accuracy-qualified.

## Execution and integration

Use one numerical owner, frozen actual source/input/dependency identities and the
central artifacts/numerical.lease plus docs/experts/execution-state.json. Claim,
separately verify, run, verify exact terminal identity, record idle under lease,
separately read back and release the lease as the final action. UTF8-sig readers
handle the ledger BOM. Source writers never run numerical work independently.

Existing shared MATLAB CAPYTAINE_REFERENCE uses engine29796/launcher40036,
started2026-10-25T17:20:15Z. These are checkpoint identities, not liveness proof.
The command window9245660 was rechecked MINIMIZED during attempt2. Keep it
minimized/headless; never open model/figure windows or quit the shared engine.
Owned KH descriptor recovery completed; importer probes left no open handles.

Active project branch is baseline/fluid-body. Preserve other worktrees. The
maintained fork branch is topwater/reference-v7.1.0. Never push main, merge or
contact upstream without separate authorization. Own issue30 updates and these
development checkpoints are authorized by the user/project instructions.

## Historical evidence

Earlier renderer/reference checkpoints, failed studies and their original
policies remain preserved in their files and the immutable
[5ef862a handoff](https://github.com/amosborne/topwater/blob/5ef862a/docs/checkpoint-handoff.md).
[Historical WEC workflows](../reference/wecsim/README.md) explicitly identify the
retired scene-cache/two-path render CLI. Historical queued work never authorizes
new work or overrides the current scope above.
