# Expert coordination

The main coordinator owns requirements, the physical case, priorities, bounded
assignments, acceptance decisions and integration. Experts own their solver's
technical investigation, implementation, reproduction and durable evidence.
Do not repeatedly inspect work already delegated; request a concise return and
challenge unsupported claims. A named expert still needs evidence and independent
review of consequential corrections.

Route Capytaine production to capytaine_expert; WEC-Sim AND BEMIO to
wecsim_expert; Topwater implementation to topwater_expert (currently held);
literature/method proposals to numerical_methods_researcher. Reference experts
jointly approve the producer/consumer contract, with BEMIO owned by its consumer.

Each task declares scope, exclusions, deliverables, acceptance and a stopping
point before execution. Return at most a short outcome/evidence/physical
consequence/unresolved issue/next-action summary. Put code maps, derivations,
reproduction logs and regression details in the owner's knowledge base.

One execution owner, assigned by the main coordinator, schedules numerical jobs
serially. Read docs/experts/execution-state.json BEFORE a numerical launch.
Only that owner may launch or alter an active job. Archive the actual code,
inputs and dependency identity before launch; keep them immutable until the
EXACT process/session handle is terminal. A timeout or stale lock file is not
terminal proof. Read-only investigations may proceed concurrently. An idle
MATLAB shared engine is not permission to start competing jobs.

No global coefficient-gate, performance or whole-library audit campaign. The
immediate reference outcome remains the first minimum-mesh ten-second incoming
wave simulation/video, followed by motion variation with mesh. Reproduce the
known WEC-Sim leads and review the paths actually exercised; resolve confirmed
result-invalidating issues or withhold the affected result. Speculative unrelated
issues do not block the entire workflow indefinitely.

Maintain current-state separately from historical notes. Cite source commit and
exact evidence; label hypotheses and decision needs explicitly. Update persistent
knowledge as part of completing each task. New chats read the short entry point,
not an entire chronological journal. No permanent model-learning claim.

Before physical corrections are relied on, assign an independent reviewer who
did not author the same fix. Generic WEC-Sim fixes go in its small maintained fork;
case geometry/forcing/render configuration stays here. Never silently carry a
marine-gyro workaround forward. Prepare upstreamable patches, but do not publish
issues/PRs, send messages or merge without separate user authorization.

Preserve current work. Commit/push development branch baseline/fluid-body, never
main; update issue30 and the fresh handoff at checkpoints. No bookkeeping PRs.
CI is PR-to-main only; local tests remain required. No PR merge without current
head required CI success. Model and reasoning settings inherit from the parent;
these role definitions deliberately contain no overrides.
