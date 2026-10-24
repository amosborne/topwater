# Activation and reuse

Installed client: Codex CLI0.164.0. Official supported standalone project files
are .codex/agents/<role>.toml, with name, description and developer_instructions.
The four project definitions omit model/reasoning overrides, inheriting their
native parent settings. No permission/provider/MCP changes were introduced.

Verified in this setup:

1. TOML parsing and exact required fields for all four definitions.
2. codex app-server --strict-config --print-startup-config exits0.
3. codex debug prompt-input renders all four names in its ACTUAL spawn agent_type
   selector, not merely text in AGENTS.md.
4. A native CLI coordinator invoked each exact custom agent_type. All four
   returned role/developer-instruction/knowledge/scope confirmations and updated
   their own durable knowledge. See native-activation.json and client-verification.json.
5. Native activation session31290 ended exit0. The canonical task IDs in the
   receipt are scoped to that activation's agent tree; they are not global UUIDs
   or guaranteed permanent live services. Forensic session metadata additionally
   binds persistent UUIDs to each role in native-spawn-receipts.json. No numerical jobs were used as a demo.

The pre-existing host chat's collaboration wrapper does not expose agent_type.
Do not claim it hot-reloaded this schema. The installed native CLI DOES expose
the selector and was the verified mechanism. Fresh native clients discover these
project definitions; use the named native role for subsequent assignments.
An app/new chat whose tool schema supports custom roles can delegate directly.

Reproduce an activation from the project directory, without selecting another
model or running a numerical case:

```powershell
codex exec --approve-for-me "Spawn capytaine_expert using exact agent_type capytaine_expert. Read its start-here/current-state and shared contract, confirm ownership and limits concisely, then stop. No numerical lease, solver edits or model overrides."
```

Replace both role names for another expert. For several roles, use a native
coordinator and limit simultaneous children to respect the active client limit;
sequential invocation is sufficient. The original activation reported max3
children, and the fourth started after the first returned. Do not treat the
four permanent definitions as requiring four concurrent threads.

Reusable context is the version-controlled knowledge base. Follow-up task prompts
need only the role's start-here path, bounded outcome and execution lease state;
do not dump a giant chronological conversation or routine solver investigation
back into the main coordinator. Read-only task planning is not numerical permission.

Current owned assignments are in assignments.md. Numerical execution remains
idle/unleased in execution-state.json. No new client restart or reload was
needed for the tested fresh CLI invocation; this does not imply schema changes
inside the already-running host conversation.
