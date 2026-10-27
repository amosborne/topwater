# WEC-Sim fork cleanup proposal

Status: authorized by the user's active validation/adoption-and-cleanup goal.
Execute cleanup after new-reference adoption; no branch or issue deletion performed yet.
Inventory checked 2026-10-27 after the user's squash merges of PR9 and PR14.
Scope is amosborne/WEC-Sim; no upstream changes.

## Intended end state

Two long-lived fork branches: main (upstream mirror) and topwater/reference
(our release-based integration). Preserve release tags and one immutable legacy
archive tag when retiring the old runtime branch. Future fixes use temporary
issue branches and focused PRs, removed after merge.

Current reference:6b774228e1f6b31dee473e626204b0ef65e5ec38, based on official
v7.1.0/a23d299 plus merged PR6, PR7, PR9 and PR14. main matches current upstream
main0753b2e. The existing Topwater pin remains9952634 on the old branch.

## PRs

| PR | Disposition |
|---|---|
| 6: Capytaine file handles | Keep merged history and final authoritative description |
| 7: Literal wave filenames | Keep merged history and final authoritative description |
| 9: Inertia ordering | Keep merged history and final authoritative description |
| 14: Spectrum phases/weights | Keep merged history and final authoritative description |
| 8: Imported excitation history | Use GitHub native Archive pull request after finalizing its superseded description |
| 10: Imported directional force dispatch | Use GitHub native Archive pull request after finalizing its superseded description |

Archive PR8/10 with GitHub's actual per-PR archival feature; closed alone does
NOT count. Native archival removes public visibility and makes discussion read-only,
while users with triage access or higher retain access and can unarchive. Finalize
their one authoritative superseded description first. No separate PR export or
archive-tag substitute is needed solely to achieve PR archival.

Verified precedent: Topwater PR13-29 are genuinely archived (`is:archived` search
returns17); merged11/12 are not. WEC-Sim currently has0 archived PRs. REST state
`closed` alone cannot distinguish archival; earlier contrary explanations were
incorrect. Official instructions:
https://docs.github.com/en/communities/moderating-comments-and-conversations/archive-pull-requests

The retired elevation-import approach was abandoned for Topwater; this does not
claim PR14 fixes its distinct unused defects. Remove links to deleted issues from
the superseded descriptions before archival. Record exact PR heads for branch cleanup.

## Issues

Keep closed issues3/4/5/12: they document merged fixes, executable before/after
examples and upstream candidacy. Keep each main body authoritative.

Delete issues1/2/11/13 after preserving one historical export:

- 1: imported excitation history, abandoned integration route.
- 2: imported directional force sum, abandoned integration route.
- 11: convolution performance proposal, abandoned approach.
- 13: directional weighting, duplicate folded into issue12.

Capture bodies/comments, referenced commits, relevant test evidence and associated
PR identities in one durable archive, then redirect current expert-note links to
that archive. Preserve original failed-study results under their original policies.
Do not create substitute archive issues. Issue deletion is permanent:
https://docs.github.com/en/issues/tracking-your-work-with-issues/administering-issues/deleting-an-issue

## Remote branches

Keep:

- main
- topwater/reference
- topwater/reference-v7.1.0 TEMPORARILY, until the pin migration gate below.

Delete after the archive/reference checks:

- fix/capytaine-sidecar-handles (merged PR6)
- fix/literal-wave-path-display (merged PR7)
- fix/rigid-hydro-inertia-order (merged PR9)
- fix/full-directional-spectrum (merged PR14)
- fix/imported-excitation-history (closed PR8)
- fix/native-imported-excitation (closed PR10)
- dev (unused fork mirror; exactly matches upstream c6b02e39)
- gh-pages (unused fork mirror; exactly matches upstream f1c51e87)

Fork Pages is disabled. Actions permissions remain enabled: inherited docs.yml can
recreate gh-pages on a later main/dev push. Retire the fork's documentation
publishing workflow through fork settings before deleting that branch; upstream
workflows and sources remain untouched. Do not confuse zero currently registered
workflow entries with Actions being disabled.

Squash merging means topic heads need not be ancestors of the new reference.
Use verified PR merged status and recorded exact head identities, not only
`git branch --merged`, when classifying branches as safe to retire.

## Old-runtime retirement gate

Do not delete topwater/reference-v7.1.0 yet. Topwater's recorded dependency and
checkout still select9952634. The new standard-input adapter is committed only in
an isolated Topwater worktree9539192; its38 Python and21native validator/handoff tests pass. Full new-route canonical
qualification remains unfinished; the first full attempt failed MATLAB startup
before physics because sandbox preferences were unwritable. Existing expanded-path
boat validation is not proof of the replacement adapter.

After the replacement is validated and its dependency-pin change is authorized:

1. Select and commit the new reference SHA and qualification evidence in Topwater.
2. Create one immutable archive/topwater-legacy-9952634 tag at exact9952634.
3. Verify the tag/archive can recover historical producer identity.
4. Delete the old runtime branch.

No rebase to a new upstream release or production pin update is part of the first
cleanup phase. Keep the v7.1.0 release tag unchanged.

## Local checkouts

Audit local dirty/untracked work, unpublished commits, active references and
archived source requirements before removing any local branch/worktree. Preserve
the active standard-input adapter worktree and original failed-study evidence.
Use one git bundle plus necessary dirty patches/source snapshots for obsolete
trial work; avoid creating many visible archive branches. Local cleanup is included in the authorized goal after this preservation audit,
not a side effect of remote pruning.

## Execution order after validated adoption

1. Export the obsolete issue/branch/PR evidence once; verify recovery.
2. Finalize superseded PR8/10 main descriptions, then natively ARCHIVE both and
   verify `is:archived` returns them. Closing alone is insufficient.
3. Delete obsolete issues1/2/11/13.
4. Retire unused fork docs deployment and delete the eight remote branches listed.
5. Verify four merged PRs/four closed fix issues, expected branches/tags and no
   accidental upstream/reference/pin changes.
6. Handle local checkouts and the old pinned branch only after their separate gates.
