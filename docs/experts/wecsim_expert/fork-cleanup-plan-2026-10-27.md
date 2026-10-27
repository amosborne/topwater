# WEC-Sim fork cleanup completion record

Completed 2026-10-27 after validated reference adoption. Independent current-state
audit passed; this record documents the completed cleanup. Scope:
`amosborne/WEC-Sim` and the audited local WEC trial checkouts. No upstream changes.

## Retained runtime and history

The adopted runtime/reference is **6b774228e1f6b31dee473e626204b0ef65e5ec38** over
[v7.1.0/a23d299](https://github.com/WEC-Sim/WEC-Sim/releases/tag/v7.1.0), with merged
PR6/7/9/14. Topwater adoption checkpoint7e1b854 is pushed. Native motion286f3a82,
private eight-file delivery4622b301 and actual default-path readiness3712ab33
provide the qualification chain; canonical boat/sea inputs are unchanged.

Only TWO remote branches remain:

| Branch | Exact retained revision |
|---|---|
| `main` (upstream mirror) | `0753b2e47f2457c078751dcfe5d251d1767b80ab` |
| `topwater/reference` (production integration) | `6b774228e1f6b31dee473e626204b0ef65e5ec38` |

The published annotated tag `archive/topwater-legacy-9952634` was verified to
resolve to **9952634f6f2e2746a1626484e27c3003d1e1d5a9**. It replaces the retired
legacy branch; the v7.1.0 release tag is preserved. Default runtime tracking is
limited to `topwater/reference`.

## Completed remote cleanup

- PR8/10 were natively archived after their authoritative superseded descriptions
  were finalized. Archived search confirms BOTH; closed state alone was not used
  as proof. Merged PR6/7/9/14 and closed fix issues3/4/5/12 remain.
- Obsolete issues1/2/11/13 were deleted and deletion verified. GitHub can return
  not-found or deleted responses (404/410); no claim that every later response
  remains404. Their essential history is preserved in the committed
  [retired-workflow record](../history/wecsim-retired-workflows-2026-10-27.md).
- Nine remote heads were atomically deleted with expected-SHA leases: the six
  retired `fix/*` branches, unused `dev` and `gh-pages` mirrors, and the old
  `topwater/reference-v7.1.0` branch. Legacy local branches were also removed
  from the default and common stores after archive verification.
- Actions was disabled ONLY for our WEC-Sim fork and read back `enabled=false`.
  Topwater Actions remains enabled. The fork has zero registered workflows;
  qualification used local tests. Inherited workflow source stays unchanged so
  `main` remains a mirror. Before deliberately re-enabling fork CI, handle the
  inherited documentation deployment, which can recreate `gh-pages`.

Evidence: [cleanup execution](../../../artifacts/workflow/wecsim-cleanup-execution1/),
[remaining heads](../../../artifacts/workflow/wecsim-cleanup-execution1/branches-after.json),
[legacy tag](../../../artifacts/workflow/wecsim-cleanup-execution1/legacy-tag.json),
[archived search](../../../artifacts/workflow/wecsim-cleanup-execution1/archived-search-confirmation.json).

## Preserved source and completed local cleanup

[Preservation](../../../artifacts/workflow/wecsim-cleanup-preservation1/) contains
ONE shared-store bundle, six exact dirty/untracked source snapshots and patches,
24 ignored Simulink file snapshots, Git configuration and the original
`shallow-boundaries.txt`. Recovery requires the bundle AND that shallow boundary.
An independent bare restore recovered 24 refs, all four unpublished trial commits,
and passed connectivity checks without network/alternates. A plain clone without
that boundary is not the verified recovery method.

After recovery verification, 11 obsolete worktrees and 13 common-store trial refs
were removed. Original failed numerical studies and their policies remain intact.
The common parent is clean/detached at releasea23d299; its linked6b validation
checkout remains, so that shared object store is retained.

WEC-related checkouts kept:

| Path under `artifacts/` | Purpose |
|---|---|
| `dependencies/topwater-wec-sim` | Clean production default, tracked `topwater/reference` at6b |
| `dependencies/wec-sim-v7.1.0` | Original stock releasea23 |
| `dependencies/wecsim-rebuild-v7.1.0` | Clean detacheda23 common-store parent |
| `dependencies/wecsim-reference-validation` | Clean linked validation6b |
| `worktrees/standard-wecsim-input` | Qualified953 Topwater source provenance |

[Local removal proof](../../../artifacts/workflow/wecsim-cleanup-execution1/local-worktrees-removed.json)
and [remaining refs](../../../artifacts/workflow/wecsim-cleanup-execution1/local-refs-after.json)
record exact paths and identities. No extra archive branches or per-worktree
bundles were created. Future maintenance follows the
[fork strategy](../../wecsim-fork-strategy.md); upstream publication still requires
explicit user permission. This cleanup does not resolve physical issues31-34.

[Independent completion audit](../../../artifacts/workflow/reference-adoption-completion-review-v1/audit.json),
SHA3968adb4, verifies actual API state, source/output hashes, process termination
and recovery data. No further cleanup operation is queued.
