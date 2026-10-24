# Branch purpose and preservation

There is one active implementation branch: **`baseline/fluid-body`**. It is an
unqualified validation foundation, not an accepted nonlinear simulator.

| Role | Reference | Use |
|---|---|---|
| Active implementation | `baseline/fluid-body` | Canonical boat, selected core, contract and blocking validation machinery |
| Released checkpoint | `main` | Existing package history; not a qualified simulator; no direct pushes |
| Preserved optimization graph | `archive/optimization-2026-10-23` at `4e8f511` | Complete prior experimental work, held |
| Supporting reset checkpoint | `validation/retire-canoe` at `8c7767e` | Reviewed generic-test decoupling; selectively extracted, not a competing runtime |
| Supporting harness work | Checkpoint named in the handoff | Selectively integrated into the active branch |
| Other branches | Listed in `branch-inventory.json` | Historical; no implied authorization to resume them |

[The inventory](branch-inventory.json) records each ref/SHA, worktree path and
dirty/untracked paths with hashes at its stated snapshot. Inspect live Git
state before use; the handoff records the current implementation checkpoint.

No worktree was removed or dirty work reset. Existing historical worktrees can
contain uncommitted work and must remain intact. One untracked benchmark file
was preserved, with its hash verified, under
`artifacts/retirement-preserved/benchmark_rankine.py`; its committed origin is
also retained in `8c7767e`.

Do not spend the reference-convergence goal on wholesale branch cleanup.
Before later retiring any ref, preserve its exact commits remotely and inspect
the owning worktree. Never delete active worktrees or discard uncommitted work.
