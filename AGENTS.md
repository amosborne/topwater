# Project coordination and expert routing

## Current boundary

Use the durable project experts below; current task/status belongs in the
handoff. The project milestone is an end-to-end external reference:
Capytaine -> BEMIO -> WEC-Sim, ten seconds of
six-DOF motion from rest in locally calm water, waves visibly rolling in from
upstream, an isometric video, and motion variation with solver mesh resolution.
No Topwater solver-specific code work is authorized in this reference phase.
Read docs/experts/coordination.md and docs/experts/shared-contract.md. Temporary
status and exact next action belong in docs/checkpoint-handoff.md and the expert
current-state files, not this stable routing document.

## Delegate technical ownership

Use the project expert definitions in .codex/agents, inheriting the selected
parent model/settings without overrides. Main thread coordinates scope,
physical requirements, priorities, bounded assignments and acceptance; it should
not routinely repeat expert implementation investigations.

- capytaine_expert owns geometry/solver meshes, native hydrodynamics, exports,
  versions/conventions and producer-side mesh studies.
- wecsim_expert owns BEMIO, the maintained WEC-Sim fork, time-domain simulation,
  startup/excitation/radiation memory, model construction, postprocessing and video.
- topwater_expert owns later solver implementation, feasibility, tests and complete
  performance measurement. Establish knowledge now; implementation remains HELD.
- numerical_methods_researcher owns focused primary-literature/method proposals,
  assumptions/tradeoffs and cheap rejection experiments; no independent solver
  implementation or broad optimization campaign.

Each expert reads docs/experts/<role>/README.md and current-state.md at startup
and updates durable evidence after substantive work. Keep current state concise,
separate from detailed historical notes. Findings cite source version, exact
reproduction and limits; named expertise is not evidence or permanent learning.
Return compact outcome/evidence/physical consequence/unresolved issue/next action.
Independent reviewers are separate from authors of consequential fixes.

## Physical and numerical scope

The canonical24ft boat is shared DATA, neutral input at
reference/capytaine/canonical-boat.json. Keep minimum full hull geometry separate
from derived solver meshes, preserving physical surface/equilibrium/world CG/full
dry inertia (including products), density1025, gravity9.81, infinite depth and
zero forward speed. The old canoe and archived P1 implementation do not qualify
this boat. Preserve generic mathematical coverage and all old failed studies.

The first incoming-wave result has NO imposed temporal amplitude ramp or global
load/amplitude step. Waves start upstream; actual excitation and displayed
incident-wave field must share arrival timing/directions. Retain startup; reaching
steady state within10s is not required. Do not animate an arriving front over
unrelated steady/globally forced body motion, or label prescribed incident water
as a computed scattered fluid field. First obtain minimum-mesh end-to-end motion,
then characterize changes in all SIX time-domain histories with refinement.
Intermediate A/B/Ainf/coupling differences are diagnostics, not independent
coefficient thresholds. Never relabel old failures under their original policies.

Review the exercised workflow proportionally. Independently reproduce marine-gyro
leads on the selected WEC-Sim version, then resolve confirmed result-invalidating
issues or withhold the affected result. Unrelated speculation and whole-library
or coefficient audits must not replace the first working reference milestone.
BEMIO belongs to the WEC-Sim expert; both reference experts approve its interface.
Generic engine fixes belong in the small maintained fork; case-specific geometry,
forcing and rendering stay in this project. Do not blindly copy old workarounds
or patch generated models with strings when a proper engine/library fix is feasible.

## Execution and checkpoint discipline

Every task declares scope, exclusions, deliverables, acceptance and a stopping
checkpoint. One coordinator-assigned execution owner schedules numerical jobs
SERIALLY using docs/experts/execution-state.json. Archive actual producer/input/
dependency identity before launch. Keep them immutable until the exact process
handle is terminal; a timeout or stale lock is not proof. Read-only work may run
concurrently. No numerical run merely to demonstrate agent activation.

Preserve active worktrees, current uncommitted work and checkpoints. Keep one
active implementation branch baseline/fluid-body; never push main or create
bookkeeping PRs. Commit/push, update issue30 and write a fresh-context handoff at
checkpoints. Supported boats change through data, not solver/reference assumptions.
CI runs only for PRs targeting main. Continue appropriate local tests. Never merge
or squash-merge without required checks passing for the PR's current head;
missing/pending/skipped/failed required checks cannot permit a merge.

The broader Topwater physical-comparison/runtime goals are held. External
reference loads or trajectories must never substitute for independent Topwater
fluid physics when that implementation phase is later authorized. Major speed
campaigns remain held until the required physical comparison; local residuals,
energy checks, fitted loads and ODE tolerances do not establish motion accuracy.
