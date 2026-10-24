# Checkpoint handoff

Updated 2026-10-24. Read this before choosing work; latest user instructions take
precedence. Stable startup/routing rules are in [AGENTS.md](../AGENTS.md).

## Authorization and phase

Startup/delegation cleanup is complete and independently reviewed. The
development-branch checkpoint records this cleanup; no implementation task is active.
Reference simulation work awaits a future user request; queued assignments do
not authorize resumption. Topwater solver implementation, tests, numerical runs
and optimization remain HELD during the external-reference phase.

The pending milestone is minimum-mesh Capytaine -> BEMIO -> WEC-Sim: ten seconds
of six-DOF motion from actual rest, incoming upstream waves, isometric video,
then motion variation with mesh. Physical requirements are in the
[shared contract](experts/shared-contract.md); no motion/video result exists.

## Preserved work and execution

- Project branch: baseline/fluid-body. Preserve untracked
  reference/capytaine/export_case.py: unfinished/unrun, SHA256
  a1feb6e8581695c1a7562bf1e129bbd609f96b54ab7dcac25da9b5ecdaa92942.
  Do not discard, stage as accepted export code or launch without owner review.
- [Execution state](experts/execution-state.json) records idle/no lease. No
  numerical jobs are authorized by this startup cleanup. Old MATLAB engine
  observations are historical; recheck only when relevant to authorized work.
- Earlier quadrature session36762 was user-interrupted, terminal exit1; preserved
  partial data is not a task to resume or a numerical failure to relabel.

## Owners and evidence

Four project roles were natively activated at setup checkpoint 9011cfe.
[Activation/reuse](experts/activation.md) distinguishes historical receipts from
current client capabilities. No repeated activation campaign is needed.
[Assignments](experts/assignments.md) records completed setup and queued follow-up;
read only the relevant expert's README/current-state when assigning work.

WEC-Sim fork/pin is recorded in [wecsim-fork.json](experts/wecsim-fork.json).
No engine correction or numerical validation exists. The four marine-gyro leads
remain unconfirmed on the selected release; detailed reproduction plans belong
to the WEC-Sim expert. Producer/consumer executed signoff remains pending.

Verification: four TOMLs parsed with required fields and no overrides; changed
JSON/local Markdown links and whitespace checked; exporter hash unchanged. An
independent fresh-context read-only review passed after correcting lease release
ordering. A temporary non-numerical CreateNew probe rejected a duplicate claim.
No new documents, solver changes, numerical jobs or activation campaign.

Next permitted action after this cleanup checkpoint: respond to the user's next
request. If reference work is requested, bound the relevant queued assignment
and follow the execution claim/recovery procedure before any numerical launch.
