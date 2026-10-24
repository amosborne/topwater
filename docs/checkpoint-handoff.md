# Checkpoint handoff

- **Current goal and qualification status:** bounded validation-foundation reset.
  Reset deliverables are complete; the Capytaine reference, selected fluid-body
  implementation, nonlinear accuracy and performance remain **unqualified**.
- **Active branch and commit:** `baseline/fluid-body`; verified numerical/source
  checkpoint `3bd62bc5407fef1211032f74b634bd93bd3fee5e`. Subsequent checkpoint
  commits contain only durable documentation, evidence and inventory. Inspect
  live `git status` and `git log` before continuing.
- **Implementation under qualification:** existing P0 Rankine-panel core,
  through `topwater.inertia.instantaneous_added_mass` and its `bem`/`boundary`/
  `rankine` owners. Extraction origin `4e8f511`. Hydrostatics and instantaneous
  inertia exist; a coupled time-domain wave path does not. The archived P1
  diagnostic runtime is a different implementation and is not substituted.
- **What passed:** 199 owned unit, mathematical, fixture and gate-policy tests.
  CI passed Ubuntu Python 3.12/3.13/3.14, Windows Python 3.12 and packaging at
  the source checkpoint. Physical-input freeze succeeded. Canoe application
  imports/runners/fixture/image are absent from active code and CI. Generic
  hydrostatics and mathematical BEM coverage are retained.
- **What failed or remains unproven:** missing-reference inertia comparison,
  unsupported-wave comparison and overall qualification returned their required
  nonzero exits. Physical qualification CI is intentionally red. No Capytaine
  convergence study ran; no physical accuracy or speed claim follows.
- **Authoritative evidence and commands:** [issue #30](https://github.com/amosborne/topwater/issues/30),
  [reset checks](reset-checks.json), [validation sequence](validation.md),
  [independent review](reset-review.md), and
  [CI source-checkpoint run](https://github.com/amosborne/topwater/actions/runs/40055853746).
  Full owned local artifacts are `artifacts/reset-units-3bd62bc.*` and
  `artifacts/validation-reset/{frozen,waves,inertia,qualification}.*`.
- **Numerical controls and acceptance:** [canonical contract](../validation/canonical-boat.json).
  Actual frozen input snapshot: [frozen physical inputs](../validation/frozen-physical-inputs.json).
  Physical identity: `6b9b330b5d1b666875b675a1136623c8f234465c5c54ea14fbe65e48be3553c3`.
  The calm wet mesh has 17 panels; study levels have 68/272/1088/4352 panels.
  Finest estimated workspace is 3,846,176,768 bytes, below the 8 GiB estimate cap.
  Four declared levels, all 36 A-infinity entries, 0.10/0.25/0.50 Hz radiation
  and all combinations of 0/45/90 degree waves, full complex excitation/RAOs.
  Reference 0.1% + 1e-6 normalized; candidate 1% + 1e-5 normalized, charged for
  reference uncertainty. These are linear policy choices, not minute accuracy.
- **Uncommitted work and other worktrees:** root is required clean at final
  push. [Branch inventory](branch-inventory.json) records historical dirty and
  untracked paths with hashes; those worktrees were not reset or removed. The
  old optimization graph is pushed at `archive/optimization-2026-10-23`.
  Supporting retirement work is pushed at `archive/retire-canoe-preparation`.
  All harness source changes were selectively integrated; do not merge entire
  historical branches. One untracked benchmark file was preserved under
  `artifacts/retirement-preserved/benchmark_rankine.py` instead of deleting it.
- **Live jobs:** none remain from this reset. Every launched local command has
  returned terminal status; no Capytaine convergence or optimization job is live.
- **Adversarial review:** fresh reviewer `baseline_adversarial_audit` found no
  remaining blocking issue after explicit identity, generic coverage, physical
  entry scaling, uncertainty/resource controls, stale-validator checks and
  stable physical-identity fixes. Approval is for the reset only.
- **Next bounded goal:** execute the frozen **Capytaine reference convergence**
  protocol and publish an explicit qualified/unqualified decision with measured
  uncertainty. Do not start Topwater wave qualification, new formulations,
  nonlinear qualification or timing in that goal.

## Exact next action and stopping condition

Inspect current repository, contract and owned freeze sources. If unchanged,
run serially, with fresh output and immutable producers:

```powershell
uv run --locked --group correlation python -m tests.validation reference --frozen artifacts/validation-reset/frozen.json --output artifacts/validation/reference-study-001.json
```

If the local freeze is unavailable, regenerate it using the documented `freeze`
command and inspect its identity before reference generation. Do not reuse stale
evidence or overwrite checkpoints.

Stop once the declared study completes or hits a declared resource/validity
failure. Retain all raw matrices, case coverage, mesh increments and uncertainty;
record a qualified/unqualified reference decision. An unqualified result is a
completed bounded investigation, not baseline qualification: identify the next
bounded reference problem and hand off in a fresh conversation. Do not silently
drop a finest level, loosen criteria or start another formulation.
