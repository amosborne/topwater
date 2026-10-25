# Project working instructions

## Fresh context

1. Read docs/checkpoint-handoff.md first. Reconcile it with the latest user
   request; the latest request controls scope. A milestone or queued assignment
   is not permission to resume it.
2. Check `git status --short --branch` and docs/experts/execution-state.json.
   Preserve existing changes/worktrees. Recorded process status is a checkpoint,
   not proof of current liveness; reconcile it before numerical work.
3. Read docs/experts/coordination.md. For technical work, also read
   docs/experts/shared-contract.md and only the relevant expert's README.md and
   current-state.md. Follow historical evidence links only when needed.
4. Delegate a bounded task to the appropriate loaded role. Do not repeat setup
   or invoke every expert just to prove activation. If a required role is absent,
   follow docs/experts/activation.md rather than claiming it is loaded.

## Ownership

The coordinator owns requirements, scope, priorities, acceptance and integration.
Use the project roles in .codex/agents, inheriting parent model/settings without
overrides; do not routinely repeat their technical investigations.

- capytaine_expert: geometry, solver meshes, native hydrodynamics and exports.
- wecsim_expert: BEMIO, maintained WEC-Sim fork, dynamics, startup and producer
  visualization adapters.
- visualization_expert: shared PyVista/VTK visualization, standard scene inputs,
  cameras, lighting, animation and measured rendering performance.
- topwater_expert: Topwater implementation, feasibility, tests and performance
  when that phase is authorized.
- numerical_methods_researcher: focused primary-source method proposals,
  assumptions, tradeoffs and cheap rejection plans; no independent solver work.

Each assignment states scope, exclusions, deliverable, acceptance and stopping
point. Experts maintain concise current facts and versioned evidence, returning
outcome/evidence/physical consequence/unresolved issue/next action. Independent
reviewers must not author the consequential physical correction they review.

## Constraints

Current phase and authorization belong in the handoff; physical requirements
belong in the shared contract. Roles confer ownership, not permission or evidence.
Supported boats change through data. Preserve generic mathematical coverage and
old failed studies under their original policies. External reference loads or
trajectories never substitute for independent Topwater physics.

Numerical jobs run serially under one coordinator-assigned execution owner;
follow the claim/recovery procedure in docs/experts/coordination.md. Freeze actual
producer/input/dependency identity until the exact process is terminal. No
numerical job merely to demonstrate activation.

Keep the active project implementation branch baseline/fluid-body; preserve other
worktrees. At substantive authorized checkpoints, commit/push this development
branch, update issue #30 and refresh the handoff. Never push main or create
bookkeeping PRs. Upstream issues/PRs, outside messages and merges require separate
user authorization. CI runs only for PRs targeting main; run appropriate local
checks. Never merge without required checks passing for the PR's current head.
