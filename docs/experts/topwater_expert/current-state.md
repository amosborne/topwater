# Current state: Topwater expert (implementation HELD)

Follow the [handoff](../../checkpoint-handoff.md) for current authorization.
Setup is complete; no reference or Topwater implementation work is authorized
by a role activation or queued assignment.

2026-10-24 native activation confirmed as `topwater_expert`, actual agent ID
`/root/topwater_setup`, with native developer instructions and inherited settings.
Documentation-only setup completed; see [held evidence/code map](held-source-map.md)
for exact read paths, governing equations, source checkpoints, evidence limits and
missing capabilities. No source inspection, tests, numerical jobs or MATLAB runs.
That setup changed only this expert's knowledge files. No role setup blocker
was identified; consult execution-state.json for scheduling, not this history.

Role setup and durable knowledge only. User expressly forbids Topwater solver
code work during external-reference phase. No simulation, optimization, solver
edits or formulation experiments are authorized by creating this role.

Historical selected implementation: P0 Rankine-panel core extracted from4e8f511,
instantaneous_added_mass with bem/boundary owners. Read-only code map:
src/topwater/inertia.py, bem.py, boundary.py, rankine.py, fluid_boundary.py,
pressure.py. Current record says hydrostatics and instantaneous inertia exist;
complete coupled wave-response RHS and derivative recovery are missing. That
record must be rechecked before future implementation, not adopted as timeless.
Archived P1/canoe runs do not qualify this24ft boat.

Historical evidence: docs/validation-reset.md, mixed-pressure-001.md,
time-domain-baseline-goal.md (HELD), reference-convergence.md (SUPERSEDED), and
relevant validation/evidence packets. Preserve useful generic math tests and all
old failures. External hydrodynamics/trajectories are never runtime substitutes
for independently calculated Topwater loads and motion.

Completed bounded assignment: read the role, shared contract and existing
documentation evidence; recorded the held-state/code/evidence map and confirmed
ownership and limits. Setup is at its stopping checkpoint. Await assignment;
no Topwater code changes or numerical jobs are authorized.
Later, only when authorized: feasibility review of researcher proposals,
bounded implementation/integration, uncertainty and complete measured runtime.
