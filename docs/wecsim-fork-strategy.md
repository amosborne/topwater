# WEC-Sim fork maintenance strategy

**Do not create upstream issues or PRs, or contact upstream maintainers, without
specific, explicit user permission.** Preparing, reviewing or merging changes in
our fork does not authorize upstream submission.

## Current integration

Our fork is `amosborne/WEC-Sim`. The user merged PR6/7/9/14 into the single
release-based integration branch,
[`topwater/reference`](https://github.com/amosborne/WEC-Sim/tree/topwater/reference),
now **6b774228e1f6b31dee473e626204b0ef65e5ec38** over official
[v7.1.0/a23d299](https://github.com/WEC-Sim/WEC-Sim/releases/tag/v7.1.0).
`main` is an upstream mirror, not a second production backend. Topwater uses
ONE exact validated commit, recorded in [wecsim-fork.json](experts/wecsim-fork.json).
A merge to the integration branch does not automatically move that pin.

| Merged fix | Generic purpose and upstream candidacy |
|---|---|
| [PR6](https://github.com/amosborne/WEC-Sim/pull/6) | Close CAPYTAINE sidecar handles on success and failure. Generic importer defect. |
| [PR7](https://github.com/amosborne/WEC-Sim/pull/7) | Print filenames literally rather than interpreting path characters as formatting. Generic display defect. |
| [PR9](https://github.com/amosborne/WEC-Sim/pull/9) | Correct inertia-product order at the active Rigid Hydro File Solid. Generic mapping defect; Drag/Flex are outside this correction. |
| [PR14](https://github.com/amosborne/WEC-Sim/pull/14) | Honor supplied FullDir phases and use consistent direction weights for elevation and forces. Generic spectral-input defects. |

These are upstream candidates, not permission or a guarantee that the current
commits suit upstream's current version. Before an authorized submission, check
its actual contribution target, current source and existing issues/PRs; reuse a
relevant existing issue instead of creating a duplicate. State measured coverage
and missing compatibility evidence honestly.

The adopted Topwater adapter uses WEC-Sim's EXISTING public `bodyClass(hydroStruct)` spectral
input. WEC-Sim computes excitation and motion. Boat/sea configuration and derived
coefficient preparation stay in Topwater; the fork needs no Topwater-specific
coefficient view/accessor or optional excitation-IRF writer extension. Prefer
existing native interfaces. A proposed capability extension needs a demonstrated
failure of those interfaces and a reviewed engineering rationale.

Adopted2026-10-27: the default pin is6b774228. The active Topwater adapter matches
reviewed953; full native force/motion, fresh total-field delivery and actual default-path
preflight have independent acceptance. The canonical eight files now match
the accepted private result, with boat/sea unchanged and verified rollback retained.
Cleanup remains the next authorized phase; the old995 branch stays until its
archive is verified. Current evidence and phase belong in the
[handoff](checkpoint-handoff.md) and [expert state](experts/wecsim_expert/current-state.md).

## Future fixes and promotion

1. Open a focused issue in OUR fork, with a plain problem/consequence, smallest
   executable same-input before/after example, independent expected result and
   acceptance criteria. Keep one authoritative description; rewrite it when
   scope changes rather than appending progress reports.
2. Create one focused topic branch from the intended integration base. Keep
   logical fixes separable; a coupled correction can use one issue/PR. Declare
   dependencies explicitly. Case adapters belong in Topwater.
3. Open a PR with repository **amosborne/WEC-Sim** and base **topwater/reference**
   explicitly selected; GitHub defaults may otherwise target upstream. Reference
   the issue and briefly describe changes/files, runnable tests, actual coverage
   and deficiencies. Use accepted [issue4/PR6](https://github.com/amosborne/WEC-Sim/pull/6)
   and [issue5/PR7](https://github.com/amosborne/WEC-Sim/pull/7) as writing models.
4. Review and test the PR's exact head, reproduce the defect before the fix, then
   squash into our reference when separately authorized. Retain tested identities
   and evidence; remove the temporary branch after merge. Do not rewrite tested
   published history merely to make the branch list tidy.
5. Promote an exact combined pin only after the affected selected-path checks
   pass. For a physical-path change, verify actual native force application,
   startup/rest/radiation, motion and required output compatibility under the
   [shared contract](experts/shared-contract.md). Check the adopted DEFAULT path,
   not only a private checkout override. Preserve rollback provenance and the
   previous validated pin until adoption is complete.

Tests must expose the actual defect, not only repeat the implementation. SLX
changes require intended-property/unrelated-content checks and an actual native
component witness. Preserve file-format compatibility claims separately from
runtime compatibility. Required local selected-path regressions do not imply
coverage of every upstream MATLAB/platform/variant combination; disclose untested
coverage rather than claiming it or requiring the whole matrix for a small proposal.

## Release maintenance and retirement

Upgrade only for a concrete reason and with separate authorization. Evaluate a
chosen OFFICIAL RELEASE in an isolated candidate, replay the minimal generic
patches, retain before/after regressions and requalify the exercised Topwater path.
Do not automatically follow main/dev, merge untested development ancestry into
the release integration, or silently move a runtime pin. An upstream fix accepted
on dev remains a local backport until the chosen release contains and verifies
it. Rebase disposable proposals if needed; do not rewrite the validated production
reference. Record semantic backport differences when an upstream proposal uses
a newer base, and retire a local patch only after its replacement passes checks.

The imported-elevation/convolution route is retired for Topwater. Its results,
limitations and abandoned fixes remain in the committed
[retired-workflow record](experts/history/wecsim-retired-workflows-2026-10-27.md);
PR14 does not claim to fix those unused paths. AFTER validated adoption, follow
[fork cleanup](experts/wecsim_expert/fork-cleanup-plan-2026-10-27.md): keep merged
PR history/issues3/4/5/12, natively archive PR8/10, delete obsolete issues1/2/11/13
and retired branches, and preserve release/legacy pins plus dirty/unpublished
study source before local cleanup. GitHub issues/PRs remain the per-fix history;
no second custom status ledger or archive-branch proliferation is needed.
