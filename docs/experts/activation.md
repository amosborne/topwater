# Activation and reuse

Check the current client's exposed agent types when delegation is needed. If it
supports a named project role, use that exact agent_type directly for a bounded
assignment, inheriting parent model/settings. Do not reactivate all four roles
or run numerical jobs just to prove setup. If the role is unavailable, report the
capability gap; do not silently start a separate CLI coordinator.

Roles live in .codex/agents/<role>.toml with name, description and
developer_instructions. Current authorization is in
[the handoff](../checkpoint-handoff.md), task scopes in
[assignments](assignments.md), and numerical ownership in
[execution-state.json](execution-state.json). A role or lease does not authorize
work. Follow coordination.md for scheduling and recovery. Reusable memory is the
version-controlled knowledge base; live thread IDs are session-scoped.

Historical setup on 2026-10-24 used Codex CLI 0.164.0 because that earlier host
chat lacked agent_type. Native discovery and all four invocations were verified;
see [client verification](client-verification.json),
[activation](native-activation.json) and [receipts](native-spawn-receipts.json).
Those are evidence for the versions then tested, not a current host limitation.
The role instructions were subsequently simplified; recorded historical role
hashes and activations do not certify the edited files' hashes or activation.
Do not rewrite historical receipts to imply otherwise.

Official guidance: [custom agents and delegation](https://learn.chatgpt.com/docs/agent-configuration/subagents).

## Visualization role added 2026-10-24

User requested a dedicated PyVista/visualization expert. Project role file
`.codex/agents/visualization_expert.toml` and its knowledge base now exist;
TOML required fields validate and no model/settings override is configured.
The current session's native spawn call returned
`unknown agent_type 'visualization_expert'` for an actual bounded read-only
contract review. Definition is present; current-session availability is not.
Do not claim activation or start a separate CLI coordinator to bypass this.
Recheck exposed agent types when a refreshed client/session loads the project.
Existing `/root/visualization` (topwater_expert) owns the initial library work
until an explicit stable-checkpoint transfer to the new role.
