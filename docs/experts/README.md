# Project expert agents

Read [coordination](coordination.md), [shared physical/interface contract](shared-contract.md)
and [current setup checklist](setup-checklist.md). The main thread coordinates;
solver depth, reproductions and evidence belong to the named owners.

| Role | Start here | Current boundary |
|---|---|---|
| capytaine_expert | [Capytaine](capytaine_expert/README.md) | Native producer/interface for minimum-mesh reference |
| wecsim_expert | [WEC-Sim/BEMIO](wecsim_expert/README.md) | Maintained fork, reproduction leads, consumer/startup/motion/video |
| topwater_expert | [Topwater](topwater_expert/README.md) | Knowledge ready; solver implementation HELD |
| numerical_methods_researcher | [Research](numerical_methods_researcher/README.md) | Focused primary-source proposals only; no broad survey |

Agent definitions are under .codex/agents. They omit model/reasoning overrides.
Supported format and native selector discovery are recorded in
[client-verification.json](client-verification.json). Actual activation evidence is in [native-activation.json](native-activation.json)
and [raw spawn/session receipts](native-spawn-receipts.json); all four returned.
See [activation/reuse](activation.md), [owned assignments](assignments.md) and
[maintained fork pin](wecsim-fork.json). Do not equate
file creation, selector discovery or library loading with a completed simulation.

The user's original setup instruction is preserved in
[setup-instructions.md](setup-instructions.md). Durable knowledge lives in this
repository, not permanent model learning or an unrecorded conversation.
