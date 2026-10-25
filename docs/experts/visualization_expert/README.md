# visualization_expert: start here

Read [handoff](../../checkpoint-handoff.md), [current state](current-state.md),
[coordination](../coordination.md) and the relevant [shared contract](../shared-contract.md).
Work only within a bounded coordinator assignment.

[Role definition](../../../.codex/agents/visualization_expert.toml).

Own the shared producer-independent PyVista/VTK visualization workflow and its
rendering performance. Producer experts own physics and source-specific adapters;
Topwater expert reviews package integration. Use standard rendering capabilities,
cached prepared inputs and a persistent scene, rather than a custom renderer.

Keep current facts concise; retain source versions, timings, images, limitations
and reproductions in linked evidence. A configured role is not proof that this
session can spawn it. See [activation](../activation.md).
