# Topwater: validation reset and bounded goals

## Direction

Establishing a simple, correct fluid–body baseline that matches Capytaine is the highest priority. The active application fixture is the canonical 24 ft boat. Retire the canoe from active application fixtures and validation; preserve its history in Git.

There is currently no demonstrated Capytaine-validated, full nonlinear 6DOF baseline for the canonical boat and observed Southern California sea state. Do not describe narrow checks, historical results or experimental implementations as that qualification.

Preserve the optimization work, but hold new formulations, optimization experiments, nonlinear qualification and timing work until the fluid–body validation gates pass. Corrections exposed by those gates are the immediate development work.

Keep reorganization limited to what supports physical validation. Avoid another prolonged investigation of project organization.

## Additions to AGENTS.md

- Name the single implementation under qualification, its commit and its qualification status. The word "baseline" alone is insufficient.
- Each goal must define scope, exclusions, deliverables, acceptance criteria and a stopping checkpoint before work starts.
- Keep one active implementation branch. Classify other branches as supporting work or historical; preserve commits before retiring references. Do not delete active worktrees or discard uncommitted work.
- Issue #30 owns baseline qualification and promotion requirements, plus a concise optimization table. Each avenue has one status: **open**, **in-progress**, or **closed**, with rationale and evidence.
- After numerical or formulation changes, the required physical regression gates must pass before promotion. Missing dependencies or cases, skipped required comparisons, stale source/input identity and unresolved reference uncertainty cannot count as a pass.
- Preserve generic mathematical and physics coverage when retiring obsolete application cases. Do not merely rename the canoe or silently delete useful coverage.
- At each checkpoint, commit and push the work, update the tracker, and produce a concise handoff sufficient for a fresh conversation.
- A completed investigation may conclude that a candidate failed. Record that result explicitly; it does not qualify the baseline.

These additions supplement the existing validation-before-experimentation instructions, including serial numerical jobs, immutable producer/input ownership while running, adversarial review, one simple runtime path and no pushes to main.

## Immediate bounded goal

> Restore a controlled, validation-first development foundation for Topwater.
>
> Establish one clearly identified fluid–body implementation under qualification for the canonical 24 ft boat. Retire the canoe from active application fixtures and validation. Preserve optimization work in Git, organize branches by purpose, and remove ambiguity about which implementation is being developed.
>
> Rewrite issue #30 around the evidence required to qualify and change the baseline. Keep optimization avenues concise, with status, rationale and evidence links. Explicitly state which physical validations remain incomplete.
>
> Freeze a reproducible canonical validation contract covering physical inputs, reference conventions, selected cases, numerical controls, uncertainty assessment and acceptance criteria. Define the Capytaine convergence study and subsequent Topwater comparisons for all six added-mass responses and selected small-amplitude wave responses.
>
> Use an independent adversarial reviewer to audit the organization, qualification claims and validation contract. Resolve material findings.
>
> **Completion checkpoint:** the repository, tracker and branch inventory identify one implementation and an executable validation sequence; retirement is verified; the review is resolved; all changes are committed and pushed. Produce a fresh-context handoff naming the exact next validation goal.
>
> New formulations, optimization experiments, nonlinear qualification and performance claims remain on hold.

## Validation bar

1. Freeze the physical problem: hull, equilibrium pose, world CG, full dry inertia, density, gravity, depth, axes, mode order, wave direction, phase origin and elevation-amplitude convention.
2. Establish Capytaine mesh convergence on the same physical wetted surface. Report all 36 added-mass entries, selected finite-frequency radiation and excitation results, coupled wave responses and their measured mesh uncertainty. Three levels permit a trend assessment; they do not automatically prove convergence. Use Richardson estimates only when justified.
3. Validate the actual Topwater implementation against the qualified reference. Include all six responses and couplings, selected small-amplitude wave cases and complex phase. Match physical problems and distinguish reference uncertainty from Topwater spatial, domain, temporal and startup error.
4. Make those comparisons blocking regression gates after formulation changes. Record criteria before judging results; do not relax them to obtain a pass.
5. Only then return to nonlinear 6DOF qualification and performance work.

Capytaine frequency-domain results are external validation data. Interpret them consistently against Topwater's time-domain behavior; they are not permission to substitute a frequency-domain runtime solver or add fixed hydrodynamic loads on top of the same pressure solution.

The eventual motion requirement remains at most 0.1 degree in each roll, pitch and yaw and 1 cm in each surge, sway and heave throughout a full 60-second simulation, against a meaningfully converged reference. Local tolerances and short runs do not prove this.

The eventual performance target remains at most 20 wall seconds per simulated second on the development machine, with real time preferred. Timing follows qualification and uses matched physical, numerical and resource contracts.

## Subsequent goals and fresh-context checkpoints

Use separate goals and conversations rather than a single indefinite development thread. Replace "exhaustively research and experiment" with bounded questions, declared tests and explicit decisions.

| Goal | Deliverable and checkpoint |
|---|---|
| Capytaine reference convergence | Execute the declared canonical-boat study; publish measured uncertainty and an explicit qualified/unqualified decision. If unqualified, identify the remaining convergence problem and the next bounded study. |
| Linear fluid–body validation | Compare Topwater with the qualified reference, resolve discrepancies, and establish blocking physical regressions. State any remaining failed gates explicitly. |
| Nonlinear 6DOF qualification | Demonstrate the required motion accuracy over the prescribed duration, including relevant domain, spatial and temporal controls. |
| Performance improvement | Test one justified improvement at a time against the qualified baseline. Retain it only when correctness and complete cost support promotion. |

An investigation can finish with a negative result. The broader qualification requirement remains open until its acceptance criteria actually pass. Do not call an unqualified reference or failed implementation complete merely to close a conversation.

A context checkpoint can also occur while a qualification goal remains incomplete: freeze and record the exact state, then continue that same goal in a fresh conversation. Clearing conversational context does not require pretending the goal was achieved.

## Handoff required at every checkpoint

Keep durable information in repository files and issue #30. A fresh conversation should need those sources and the following short handoff, rather than the prior conversation's narrative.

```markdown
# Checkpoint handoff

- Current goal and qualification status:
- Active branch and commit:
- Implementation under qualification:
- What passed:
- What failed or remains unproven:
- Authoritative evidence and reproduction commands:
- Numerical controls and acceptance criteria:
- Uncommitted work or other worktrees that must be preserved:
- Live jobs and exact handles, or confirmation that none remain:
- Adversarial review findings and resolution:
- Next bounded action and stopping condition:
```

The receiving agent must inspect current repository and external state before continuing. If source or inputs have changed, stale evidence cannot qualify the new state.