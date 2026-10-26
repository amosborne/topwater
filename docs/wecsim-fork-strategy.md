# WEC-Sim fork maintenance strategy

Reviewed 2026-10-25. **Keep the working pin unchanged. Prepare small upstream
patches separately, and keep Topwater case configuration outside the fork.**
The existing fork needs submission/compatibility rework, not an immediate
replacement of its selected-path corrections.

The authorized issues-first checkpoint is complete: issues were enabled only
on our fork, duplicates were checked, and the five work items below were
created and their bodies read back. The coordinator created the new clean
combined branch `topwater/reference` at official v7.1.0; runtime sources and the
existing Topwater pin remain frozen at `9952634`. No topic branch, PR, local
fork change, numerical run or upstream publication was made by this issues task.
Further implementation and execution need bounded coordinator assignments.
The first subsequent bounded rebuild is complete for issues #4 and #5:
[internal PR #6](https://github.com/amosborne/WEC-Sim/pull/6) and
[internal PR #7](https://github.com/amosborne/WEC-Sim/pull/7), each independently
based on v7.1.0 and targeting our `topwater/reference`. All nine focused tests
passed on Windows/MATLAB R2025b; stock defects were reproduced. See
[evidence and limits](experts/wecsim_expert/fork-rebuild-regressions-2026-10-25.md).
The user squash-merged PR #6 as `76adb29201291ffdec3e805aba01151d762328b9`;
PR #7 remains open. The held Topwater pin is unchanged.
Issue #1 subsequently became [internal PR #8](https://github.com/amosborne/WEC-Sim/pull/8),
also independently release-based and unmerged. Independent physical/source review
and eleven Windows/R2025b synthetic tests passed, with stock defects reproduced.
Its generic finite-record policy warns and zero-extends instead of imposing the
old Topwater-specific short-record fatal gate or negative-time clamp. See
[history rebuild evidence](experts/wecsim_expert/imported-history-rebuild-2026-10-25.md).
**Do not create any upstream issue or PR, or make any outside contact, without
the user's specific, explicit permission.** Preparing or reviewing local fork
patches does not authorize submission or contact.

## Verified starting point

| Item | Inspected identity/status |
|---|---|
| Working fork | `amosborne/WEC-Sim`, `topwater/reference-v7.1.0`, clean at `9952634f6f2e2746a1626484e27c3003d1e1d5a9` |
| New combined branch | [`amosborne/WEC-Sim:topwater/reference`](https://github.com/amosborne/WEC-Sim/tree/topwater/reference), `76adb29201291ffdec3e805aba01151d762328b9` over release `a23d299`; handle fix PR #6 merged by user, no promoted Topwater pin |
| Base / latest upstream release | `v7.1.0`, `a23d2993807410fe78cc63546c1ba9d81912024f`; release published 2026-06-30 |
| Upstream default branch | `main`, `0753b2e47f2457c078751dcfe5d251d1767b80ab`, 2026-09-30 |
| Upstream development branch | `dev`, `c6b02e397903dbc9651bd7c6a3d4cf75c85486dc`, 2026-10-12 |
| Local changes | Five commits containing six logical fixes; selected linear imported-wave/rigid-body evidence, not all-variant validation |

Branch and release identities were read from the live GitHub API, not inferred
from cached search pages. See the [release](https://github.com/WEC-Sim/WEC-Sim/releases/tag/v7.1.0),
[main commit](https://github.com/WEC-Sim/WEC-Sim/commit/0753b2e47f2457c078751dcfe5d251d1767b80ab)
and [dev commit](https://github.com/WEC-Sim/WEC-Sim/commit/c6b02e397903dbc9651bd7c6a3d4cf75c85486dc).

Current main is three commits ahead of the release. In particular,
[merged PR #1653](https://github.com/WEC-Sim/WEC-Sim/pull/1653) changed the same
body library and updated regression fixtures. An upstream refresh therefore
needs a semantic model merge and new result checks; it is not a harmless text
rebase. Do not replace its library with our release-based SLX.

At both inspected main and dev revisions, the relevant `bodyClass.m`,
`waveClass.m` and `readCAPYTAINE.m` blobs still match the original release:
`b1dee271ab7160a6af5cb5ee319ca07af73dfc63`,
`b606b0632eb2ccdb484ebd09fa7f8d4d37fbe282` and
`55b3e037d10340aedb50b95c216587c9d2c9d979`, respectively. The upstream body
SLX is now blob `8116d315eb5e6a5c47649d970f3422c6a7350d80`, versus release blob
`34655aed0de2eb10a4033d5f580293180379ce91`. Read-only in-memory XML inspection
confirmed the three unreordered inertia expressions and the shared
`userDefinedFe` input remain. This is source evidence, not execution on upstream
main/dev.

## Disposition of each logical fix

“Upstream candidate” means a generic defect with useful evidence; it does not
mean the existing commit is ready to submit unchanged.

| Fix / existing commit | Disposition and reason | Evidence and current upstream/duplicate status | Independent reproducer and regression target before submission |
|---|---|---|---|
| Imported-excitation history and lag alignment; part of `2f43f54` | **Upstream candidate.** Using only the simulation-cropped elevation loses excitation contributions from past/future samples. Explicit lag times also avoid assuming a symmetric kernel. Split this from directional aggregation. Review boundary/ramp behavior and API compatibility: the local insufficient-support error must not silently establish Topwater's startup policy for all upstream users. | [Stock crop witness](../artifacts/reference-imported-excitation-stock/result.json) and [fixed impulses/asymmetric kernel](../artifacts/reference-imported-excitation-fixed-attempt2/result.json). Current [body source](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/objects/bodyClass.m) still convolves the cropped record. No exact open duplicate found. [#1268](https://github.com/WEC-Sim/WEC-Sim/issues/1268) is related imported-elevation usage, not this correction. | Synthetic six-channel excitation kernels with impulses before, inside and after the simulation interval; asymmetric lag support, two time steps, zero/nonzero ramp, long and insufficient input records. Compare to an independently written time-indexed convolution, not a Topwater force table. Preserve ordinary regular-wave behavior. |
| Imported directional-component aggregation; other part of `2f43f54` | **Conditional upstream candidate.** Omitting later imported components is generic, but the supported group semantics and representation need agreement. Aggregation and dispatch must be reviewed together. | [Stock two-component witness](../artifacts/reference-multiple-imports-stock-attempt5/result.json): second contribution omitted, maximum discrepancy 32 in the synthetic fixture; fixed all-DOF fixture error zero. Main/dev preprocessing still chooses the first wave object. Open [#1663](https://github.com/WEC-Sim/WEC-Sim/issues/1663) concerns heading-bin weights in stationary full-directional spectra, a different path; do not claim it duplicates this bug. | Two distinct imported signals/headings with analytically known six-component totals. Check one-component equivalence, wave ordering, matching/mismatching time grids, and an actual native applied-force bus. Establish expected mixed-group behavior rather than copying our rejection policy by default. |
| File Solid products-of-inertia order; `fcd7c36` | **Upstream candidate, SLX rework required.** WEC's public order `[Ixy,Ixz,Iyz]` must map to File Solid's `[Iyz,Izx,Ixy]` at the boundary. Keep signs, public properties and residual added-mass mathematics unchanged. | [Stock witness](../artifacts/reference-inertia-stock-attempt2/result.json), [fixed tensor roundtrip](../artifacts/reference-inertia-fixed/result.json), [semantic diff](../artifacts/reference-inertia-fixed/library-semantic-diff.json). Current upstream XML still has all three direct expressions. MathWorks specifies the [File Solid order](https://www.mathworks.com/help/releases/R2025b/sm/ref/filesolid.html). Historical [#981](https://github.com/WEC-Sim/WEC-Sim/pull/981) introduced products; [#1425](https://github.com/WEC-Sim/WEC-Sim/pull/1425) changes `adjustMassMatrix` theory/docs, not these three boundary expressions. Neither is this patch. | Distinct nonzero `Ixy`, `Ixz`, `Iyz`, then full dry-plus-added tensors, including existing nonsymmetric residual handling. Verify compiled Solid plus residual reproduces the requested matrix. Check all three modified block paths; selected rigid-path success does not establish drag/Flex compatibility. |
| Imported total dispatched once, stale clone cleanup and mixed-group rejection; `f962103` | **Integration-specific implementation as written; generic solution uncertain.** It assumes the preceding aggregation stores one total, so each branch must not apply that total again. The underlying multi-import defect deserves an upstream solution, but this helper's aggregate-once representation, rigid-only scope and rejection policy are design choices. Do not submit it independently of aggregation. | [Fixed native graph tests](../artifacts/reference-imported-dispatch-fixed-attempt2/result.json): fresh/repeat/regular-to-import/stale clones, one positive total branch, mixed group rejected. [Boundary explanation](experts/wecsim_expert/accepted-input-2026-10-24.md) records the repeated identical-source defect. Current upstream rigid and Flex imported sources still read the same body array. No exact open duplicate found. Topwater's later precomputed-force binding is a separate case adapter, not proof of generic dispatch semantics. | Start from untouched upstream, two imported components and the actual parent sum. Decide between per-component branch storage and a single aggregate branch. Prove applied force equals the analytic total exactly once, then exercise repeated initialization, mode changes and intended mixed groups. Do not extend a rigid-only fix to Flex without tests. |
| CAPYTAINE sidecar handle lifetime; `f07a8df` | **Upstream candidate; small independent patch.** Close the successful KH existence probe and guard each owned text-read handle on success/failure. No coefficient parsing or normalization changes. | [Real regression](../artifacts/workflow/bemio-handles-attempt3/report.json): original one-handle leak reproduced, two imports preserve 32 numeric arrays each, malformed sidecars leave no handles, 15 Windows file rename/delete checks. Main/dev retain the leak. [#1660](https://github.com/WEC-Sim/WEC-Sim/pull/1660) is already merged into dev but changes other BEMIO readers; open [#1596](https://github.com/WEC-Sim/WEC-Sim/pull/1596) is pressure processing. Neither fixes `readCAPYTAINE`. | Valid, malformed and missing sidecars; unchanged numeric structs; open-file inventory before/after; unrelated preexisting handle preserved. Add multibody sidecar naming and non-Windows coverage. Use `openedFiles` on R2025b and a version-compatible inventory helper for older supported MATLAB. |
| Literal imported filepath display; `9952634` | **Upstream candidate; one independent display statement.** Fixed `%s` formatting prevents backslashes/percent signs being interpreted as format syntax. | [Real three-path check](../artifacts/workflow/wave-display-check-attempt3/report.json): Windows/percent/UNC paths verbatim, zero warnings, no wave setup or body run. Main/dev retain the original formatting. No exact open duplicate found; search hit [#1184](https://github.com/WEC-Sim/WEC-Sim/issues/1184) concerns GUI/variant behavior instead. | Capture actual `listInfo` output on the three paths, assert no warning and exact path text. Add supported string as well as char input; verify other wave-type displays remain unchanged. No numerical simulation needed. |

The search enumerated all nine currently open issue/PR entries and used targeted
all-state searches for inertia/products, imported excitation, `elevationImport`,
`readCAPYTAINE`, `KH.dat`, `fclose` and escaped-character terms. This found no
exact open duplicate of these fixes, not proof that no older discussion exists.
Recheck immediately before any authorized submission. API status takes priority
over the cached web page that still showed #1660 open.

## Keep the adapter boundary clear

Topwater's incoming spatial train, prepared six-load table, mean-frame geometry,
full-inertia conversion, startup gates, native bus audits, source-field
reconstruction and shared-engine lifecycle belong in Topwater. In particular,
`src/topwater/data/wecsim/run_workflow_case.m` redirects one copied instance's
From Workspace binding to `topwaterForceTable` using Simulink APIs. That binding
is explicitly integration-specific; it should not become another global WEC
library behavior. Generic WEC defects remain in the maintained fork until
upstream fixes are adopted. Do not remove existing fixes merely because the
current one-placeholder-wave adapter bypasses some imported-wave calculations.

## Branches, updates and release policy

The issues-first stage created these **open issues in our fork**; all bodies
were read back and matched their prepared Markdown files. Each begins with our
repository scope, the new release base, the explicit future internal PR base
repository/branch, and the upstream-permission rule.

| Authoritative fork work item | Scope |
|---|---|
| [#1](https://github.com/amosborne/WEC-Sim/issues/1) | Imported excitation history and lag alignment |
| [#2](https://github.com/amosborne/WEC-Sim/issues/2) | Coupled directional aggregation and once-only dispatch; generic semantics unresolved |
| [#3](https://github.com/amosborne/WEC-Sim/issues/3) | File Solid inertia-product order and compatible SLX packaging |
| [#4](https://github.com/amosborne/WEC-Sim/issues/4) | CAPYTAINE sidecar handle lifetime |
| [#5](https://github.com/amosborne/WEC-Sim/issues/5) | Literal imported-path display |

Only the fork's `has_issues` setting was enabled. No other repository setting,
topic branch or PR was created. The issues include old commit provenance,
independent reproducers, local-artifact limits and acceptance/missing evidence;
they do not certify that the old patches are ready to submit upstream.

The required working sequence is:

1. Track a WEC code modification in an issue in **`amosborne/WEC-Sim`**.
   Topwater adapter/reference-qualification issues #31-34 remain in
   **`amosborne/topwater`**; do not move or duplicate them into the fork merely
   because WEC-Sim is involved.
2. Implement one focused logical fix on its own branch, preserving separable
   commits and recording any dependencies. Aggregation/dispatch may be one
   coupled proposal.
3. Open an internal PR with **base repository `amosborne/WEC-Sim`** and the
   explicitly named combined integration branch **`topwater/reference`** as its
   base. Both repository and branch must be checked before publication: fork UI/CLI
   defaults can select `WEC-Sim/WEC-Sim` upstream instead. The head repository
   is our fork and the head branch is the focused fix branch.
4. Review and test the proposed change, integrate it into our combined branch,
   validate the resulting exact combined commit for its declared scope, then
   promote that exact Topwater dependency pin. Keep the previous validated pin
   and tag immutable.
5. Only after specific user permission, create the equivalent upstream issue
   or submit a clean PR. Recheck upstream source and contribution target at
   submission time; reuse an appropriate existing upstream issue rather than duplicate it.
   Upstream publication is a distinct action from an internal fork issue/PR.

The issue-creation step is complete; implementation/PR work for #1, #4 and #5 is
recorded above. No step changes the frozen live baseline. Keep every internal PR narrowly
scoped so its logical fix can later be replayed onto a clean upstream proposal
branch without the combined branch's unrelated changes or local adapter layer.

Keep `topwater/reference-v7.1.0` as the frozen legacy baseline during the
transition. Never force-rebase it or change a running checkout. After the new
`topwater/reference` combined branch is reviewed/tested, Topwater is repinned to
its exact validated revision, and the old `9952634` SHA is preserved under an
immutable archive tag, delete the old branch. That deletion is a later terminal
step, not authorized during this issues-first task. Keep the exact fork SHA,
base, dependency versions, selected-path evidence and adapter source identity
in [wecsim-fork.json](experts/wecsim-fork.json) and run manifests. Historical
reports continue to refer to their original SHAs.

The user requires **individual upstream proposal branches plus one combined
Topwater integration branch**. This is the mandatory branch model, not an
optional way to package the changes. The new combined endpoint exists; begin
code integration only after authorization and terminal proof, in a separate
candidate checkout:

- Inspect upstream release/source changes when a concrete maintenance need or
  authorized submission justifies it. Do not automatically track, fetch, merge
  or move a pin with upstream main/dev. Future production base upgrades target
  an official release, need a good reason, and are evaluated as candidates first.
- Each independent logical fix **must have its own internal PR branch** against
  our combined branch. Split `2f43f54`; directional aggregation/dispatch may be
  one coupled proposal. A later authorized upstream PR uses a clean proposal
  branch based on an exact upstream commit suitable for the then-current target.
  Do not submit the combined integration branch as an upstream PR.
- Our **combined branch `topwater/reference`** must assemble the patches needed
  by Topwater while upstream proposals are pending. Keep generic patch commits and any
  required integration-only commits in clearly separate layers. Case/model
  adapter work stays in Topwater wherever feasible; an unavoidable fork-only
  integration layer must be identified and excluded from upstream PR branches.
- Keep the combined branch's base explicit. A PR topic based on newer upstream
  `main` or `dev` must **not be merged wholesale into the release-based combined
  branch**: that would import untested upstream ancestry. Backport only the
  intended logical patch onto the release-based combined candidate and record
  `proposal branch/base/commit -> combined base/backport commit`, dependencies
  and any semantic differences. A future combined-base upgrade targets a chosen
  official release and requires a separately evaluated candidate and promotion
  gates.
- Rebuild disposable combined candidates when useful; do not rewrite the
  published `9952634` baseline to tidy history. If maintainers request an
  updated PR branch, retain the old tested SHA/evidence. Prefer reviewed patch
  replay/backports onto the release-based candidate over upstream merges into
  the pinned branch. Rebase only disposable/topic branches.
- Keep an official release base for production. Main/dev may supply a future
  upstream PR target and source comparison, not automatic runtime upgrades.
  A release label alone does not authorize upgrading; promote only an exact
  tested SHA.
- For conflicts, replay the intended semantic edit onto the new library using
  Simulink APIs. Never resolve a binary conflict by choosing our whole SLX or
  upstream's whole SLX without checking the lost changes. Drop a local patch
  only after upstream equivalence and its regression are verified.

Every run must name an exact immutable integration commit plus matching
Topwater adapter/dependency identities in its manifest. Give each promoted,
validated combined revision an immutable tag; never move an existing validation
tag or reinterpret its evidence after a branch advances. Mutable PR/combined
branch heads are work queues, not sufficient run provenance. Existing runs
remain tied to `9952634` until a separately validated candidate is promoted.

There is exactly **one production combined pin** at a time. A new release
candidate is an alternative pending validation, not a second production
backend. A fix accepted on upstream dev remains a local backport
until the chosen release contains it and its equivalence is verified. Backport
by semantic patch replay against the older release API/model, rather than
assuming an unchanged newer-dev commit can be cherry-picked safely.

Fork issues and internal PRs are the authoritative per-fix status and history.
This document supplies the initial technical inventory; `wecsim-fork.json`
records the runtime pin and evidence links. A combined-revision PR/release
record lists its fork issues/PRs, exact patch commits and semantic backport
mapping as needed. Keep upstream acceptance, release inclusion and Topwater
promotion explicit in that GitHub history. No separate status ledger, new
manifest schema or tracking automation is needed.

Choose the eventual upstream PR target with maintainers. Recent independent
fixes have landed on main, while #1660 and current #1668 target dev; the
[contribution guide](https://wec-sim.github.io/WEC-Sim/dev/developer/getting_started.html)
does not establish one universal target for every correction.

## Compatibility and acceptance gates

The original library was saved in R2020b; our proper R2025b save migrated its
serialization. Existing review verified 146 system connectivities, 39 embedded
scripts and 531 migrated port counts, with the three intended inertia edits
and documented migration changes. That selected-path semantic review does not
prove backward loading or compile all variants.

For upstream, recreate the small model edits on its current library and retain
the R2020b save/export convention through supported Simulink tooling, not ZIP
text surgery. Inspect the exported semantic diff and reopen it. Upstream's
[documented compatibility policy](https://wec-sim.github.io/WEC-Sim/dev/user/getting_started.html#matlab-requirements)
and [inspected CI matrix](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/.github/workflows/run-tests-main.yml)
cover R2024a/R2024b/R2025a/R2025b on Windows and Ubuntu; R2020b is the library
save convention, not a claim of fully tested R2020b runtime support.

Local candidate promotion requires the following checks, proportional to the
affected path and executed only when scheduled. Minimal resource/display PRs
use their focused regressions; they do not automatically require another body
simulation. Consequential dynamics/model upgrades require the selected-path
physical and public-workflow checks:

1. Minimal stock-fails/fixed-passes regressions above, source review and exact
   patch identity; no borrowed Topwater force table for native-import tests.
2. For SLX edits, reviewed systems/scripts/ports/variants diff and load/compile
   on our supported local R2025b selected path, plus any other changed variants
   for which support is claimed. State untested variants and MATLAB versions.
3. Topwater selected path: importer normalization, full mass/CG/frame mapping,
   every applied six-load sample, zero state/memory/startup, CIC and Newton-Euler
   audits, clean shared-engine restoration and owned-file lifetime. Compare all
   six motion histories to the preserved baseline. An upstream radiation-path
   change can legitimately alter results and needs explanation/independent
   review; do not overwrite fixtures to make a regression pass.
4. Complete public solve then render, valid NPZ, truthful logs, video decode and
   physical-warning review. Interface success does not certify mesh/time
   convergence or nonlinear accuracy. No full-compatibility or public-workflow
   success is asserted by this strategy document.

Upstream's broader compatibility matrix is a separate coverage goal. We do not
currently have local R2024a/R2024b/R2025a plus Windows/Ubuntu matrix evidence for
these patches. A small local PR proposal/review does not require reproducing
that entire matrix first. For an authorized upstream PR, include the focused
reproducer, local results, model serialization changes and explicit missing
coverage; use upstream/available CI and maintainer review to obtain the wider
checks. Verify oldest-supported loading and the advertised matrix before
claiming broad upstream compatibility, rather than labeling unrun checks passed.

Roll back by selecting the previous exact integration pin and matching adapter
revision for a new scheduled run. Preserve failed candidates, manifests, logs
and original outputs. Never swap dependencies midway through a live process.

## Incremental work after review

First finish and checkpoint the current end-to-end baseline. Then, if the user
authorizes rework, prepare the sidecar-lifetime and literal-display patches with
their independent tests; these have the smallest scope. Next split and prepare
the history correction and its boundary-policy tests. Prepare the inertia patch
with compatible serialization and the full tensor/matrix gates. Finally resolve
the multi-import representation and submit aggregation/dispatch together only
if native upstream reproducers establish the intended generic behavior.

Upstream submission order after separate user approval: handles and display as
independent PRs; history; inertia when compatibility evidence is ready;
directional aggregation/dispatch after design agreement. Upstream acceptance
reduces the local stack at the next tested dependency promotion, not by silently
changing this baseline.

Open questions: native mixed-import semantics; acceptable missing-history/ramp
compatibility policy; maintainer-selected target branch and model-save tooling;
access to the advertised MATLAB/OS test matrix. Stop here: five fork issues
created and verified; no topic implementation, PR, rework or upstream contact
performed. The legacy branch remains until the tested replacement, repin and
archive-tag conditions above are satisfied.
