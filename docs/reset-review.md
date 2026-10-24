# Independent reset review

Reviewer: fresh agent `baseline_adversarial_audit`, independent of the reset
implementation. Review scope: organization, qualification claims, physical
contract, source ownership, executable gates and generic coverage. No Capytaine
study or numerical formulation experiment was run by the reviewer.

| Material finding | Resolution |
|---|---|
| Historical checks implied a baseline that did not exist | README, AGENTS, issue #30 and contract identify the unqualified P0 core and missing time-domain waves |
| Canoe results and conflicting canonical initial states | Canoe retired from active source/tests/CI; one calm-equilibrium shared case; old moving-twist fixture remains only in history |
| Experimental branch carried 477 files of accumulated work | Lean branch from main with selective provenance; complete graph preserved by pushed archive tag; no worktrees deleted |
| Generic hydrostatics coverage could be lost | Restored six world derivatives, altitude, submerged volume/CB, vertex crossing, exact wet closure and overload checks on the actual boat |
| Implicit reference constants/conventions | Explicit density/gravity, angular units, modes, phase, propagation and elevation amplitude; unsupported declarations rejected |
| Full whitening could conceal physical entry errors | Diagonal mass scaling for entrywise A/B gates; full inertia retained in RAOs |
| Refinement/uncertainty and limits were insufficiently precise | Four declared meshes; both finest increments and trend required; empirical envelope disclosed; no unjustified extrapolation, NaN limits or silent finest omission |
| Stale validation logic could reuse old results | Physical and validation producer identities checked; changed CLI/gate code invalidates evidence |
| Whole-report hash prevented equivalent fresh freezes | Stable identity binds contract/input and named array dtype/shape/bytes; timestamps/commit remain provenance; exact source checks retained |
| Generic green CI could be mistaken for qualification | Separate real physical-qualification job is intentionally red while required evidence/capability is absent; no branch-protection claim |

The reviewer found no remaining blocking source/protocol issue for the bounded
reset after these fixes. Approval is limited to this reset foundation, subject
to its final checks, checkpoint and handoff. It does not qualify the reference,
fluid-body implementation, nonlinear simulation or performance.

[Reset check evidence](reset-checks.json): 199 tests passed. Physical freeze
succeeded; missing-reference inertia, unsupported-wave and overall qualification
commands returned their required nonzero exits. CI passed all four unit-platform
jobs and packaging; the physical qualification job failed at the expected real
promotion command. No convergence study was executed.
