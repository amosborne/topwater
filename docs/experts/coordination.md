# Expert coordination

Follow the startup order in [AGENTS.md](../../AGENTS.md). The latest user request
controls scope; [the handoff](../checkpoint-handoff.md) records current
authorization and phase. Queued ownership never authorizes execution.

## Keep one owner for each fact

- AGENTS.md: stable rules, startup order and expert routing.
- checkpoint-handoff.md: current authorization, phase, preserved work and next action.
- assignments.md: bounded task definitions, owners and task status.
- execution-state.json: execution owner, frozen manifest, process and terminal evidence.
- shared-contract.md: physical inputs and producer/consumer requirements.
- Each expert's current-state.md: technical findings, blockers and evidence links.

Update the owning document rather than append conflicting corrections elsewhere.
Historical setup instructions, checklists and activation receipts are background,
not recurring tasks or current capability claims. Recheck versions when evidence
is used; durable knowledge is files, not permanent model learning.

The coordinator assigns scope and acceptance; experts own investigation and
implementation. Give each writer explicit files/responsibility and preserve
others' changes. One writer edits a shared file at a time. Read-only work may run
in parallel. Return concise findings with evidence; do not flood the coordinator
with raw logs. A separate reviewer checks consequential physical corrections.
Generic WEC-Sim fixes belong in its maintained fork; case configuration stays here.
Capytaine and WEC-Sim experts jointly approve the interface; WEC-Sim owns BEMIO.

Correctness reviews follow the actual Topwater workflow and the supported inputs
that can reach it. A case's zeros or symmetry do not exempt an executed physical
boundary from review. Distinguish consumed results, deliberately discarded
intermediates and branches the adapter never enters; do not expand into unused
library features merely because they share a file.

Every fix needs a minimal executable example with the same inputs before and
after the fix, actual incorrect/corrected results, and an automated regression.
Put the problem and reproducer in the issue, using plain language. The PR should
reference that issue rather than duplicate it, then concisely state the change,
affected area/files, test command, coverage and remaining deficiencies. The issue's
before/after example may be linked from the PR. Preserve user-edited wording.
Internal tracking details must not obscure what a human reviewer needs.
For fork fixes, state whether the change is an upstream PR candidate and why.
Keep that technical assessment separate from permission to submit upstream.

## Numerical execution and recovery

A lease schedules already-authorized work. Before any numerical launch:

1. Coordinator selects one execution owner and bounded assignment. Check the
   handoff, execution-state.json and actual recorded process identity. Unknown
   process status blocks another launch; a timeout or idle MATLAB engine is not
   terminal proof.
2. Atomically claim the local lease from the repository root:

   ```powershell
   $leaseRoot = 'C:\Users\alex\topwater\artifacts'
   New-Item -ItemType Directory -Path $leaseRoot -Force | Out-Null
   $leasePath = Join-Path $leaseRoot 'numerical.lease'
   $leaseFile = [IO.File]::Open($leasePath, [IO.FileMode]::CreateNew, [IO.FileAccess]::Write, [IO.FileShare]::None)
   $leaseFile.Dispose()
   ```

   If creation fails, stop and reconcile the existing claim; never overwrite it.
   All project worktrees/coordinators use this same local claim and the ledger
   at C:/Users/alex/topwater/docs/experts/execution-state.json, not worktree copies.
   Other clones/machines must not run numerical jobs concurrently.
3. Under that claim, update execution-state.json to `claimed`: record owner
   session/task, assignment, timestamp and archived manifest path containing exact
   producer/input/dependency identities, command, resources and expectations.
   Freeze those files before launch. Read back the ledger and verify `claimed`,
   owner, assignment and manifest before launch. Any failed preparation or ledger
   command stops the launch; orchestration must inspect each nested command's
   exit status rather than treating a completed tool wrapper as success. Python
   readers use UTF-8-sig because PowerShell may write a UTF-8 BOM.
   Launch only this job, then record `running`
   with its tool session handle or PID plus start time, host and output paths.
4. On completion, verify the exact process and any numerical child processes are
   terminal; record exit status/evidence and result location. Record `idle` with
   last-job evidence retained while still holding the claim, then release it with
   `Remove-Item -LiteralPath $leasePath` as the final operation. Do not write the
   ledger after release. Inputs/producers stay immutable until terminal.

After a crash or new context, a lease or non-idle record requires reconciliation
before launch or producer edits. Identify the actual process using the recorded
handle/start time/command, not PID alone. If launch or completion was not recorded,
inspect the owner's session/logs and relevant host processes. Missing information
means `unknown`, not idle. Resume monitoring through the same owner, or transfer
ownership explicitly after reconciliation. Remove an abandoned claim only with
evidence that no associated numerical process remains; never infer this from age.

## Guidance basis

[OpenAI's AGENTS.md guidance](https://learn.chatgpt.com/docs/agent-configuration/agents-md)
describes instruction discovery. Its
[subagent guidance](https://learn.chatgpt.com/docs/agent-configuration/subagents)
supports bounded specialist work, concise returns and care with parallel writes.
Roles inherit parent permissions; instructions are not separate sandboxes. The
serial numerical lease and physical review requirements are project policies.
