# Checkpoint handoff

Updated 2026-10-25. Latest user steering controls scope. Read AGENTS.md and
experts/coordination.md; consult shared-contract.md and the relevant expert state.
Historical queued work is not authorization to restart it.

## Current priority and completed workflow

The user's order is working end to end, then correctness, then speed. The full
public workflow now WORKS. No further performance changes are authorized now.
Native coupled dynamics remain unimplemented. No solver/render streaming.

LATEST SCOPE: user requests rigorous correctness reviews of Capytaine, BEMIO and
WEC-Sim by their respective experts, LIMITED TO OUR CODEPATH. Supported inputs
that reach a boundary count even if this boat's zeros/symmetry hide the defect;
the active Rigid Hydro inertia correction is in scope regardless of this case's
numeric impact. Never-entered Drag/Flex/multi-import-only branches are excluded.
Three read-only audits are COMPLETE: capytaine_path_audit, bemio_path_audit and
wecsim_surface_examples (WEC dynamics/lifecycle). No numerical work ran for them.
Both source corrections are FROZEN and independently reviewed: lifecycle by
the separate BEMIO/WEC expert; finite-boundary guard by the Capytaine expert.
Root's current-path-checks-attempt2/session61218 is TERMINAL SUCCESS. Immutable
old/fixed snapshots reproduce two primary-error failures and two setting failures;
corrected backend suite21PASS and MATLAB suite13PASS (2setting +11finite guards).
Exact process trees/frozen hashes verified, idle read back and lease released.
No BEM/boat dynamics or shared-engine call.
Attempt1 stopped at pytest default-temp permissions; it is terminal/released,
preserved, and attempt2 uses owned basetemp with cache disabled. Do not edit
snapshot/drivers until terminal/released. Authors ran no tests. Capytaine
owner reconciled stale current-state documentation, without physics edits.
Distinguish discarded intermediate calculations
from functions not executed; prepared forces are the current applied authority.
Stop broad generic fork/export/all-variant fixture expansion; preserve unfinished
sources and failed evidence. No pin change/merge. Optional PRs stay separate.

Every fix PR requires a minimal executable PRE/POST example and actual results.
Examples are now published/read-back verified in PRs6/7/8 with tested heads
unchanged. User found issue4/PR6 unreadable; both have been rewritten in plain
language, leading with Windows cleanup failure and measured1openfile ->0 example.
Use that human-first standard for future issues/PRs; see coordination.md.
Latest issue-writing preference: give an upstream-candidate assessment with a
reason instead of a submission disclaimer. Issue4 and PR6 now say candidateYES:
generic importer cleanup with no Topwater-specific behavior. Both read back
verified. This assessment does not authorize upstream contact or publication.
User then refined the format: issues own problem/reproducer; PRs reference the
issue and state change, affected area, test command/coverage and deficiencies.
Issue5/PR7 have now been rewritten and read-back verified in that format; PR7
code head206bac3 is unchanged. Preserve user edits to issue4 and PR6.
User squash-merged PR6 on2026-10-26T00:27:39Z as76adb29201291ffdec3e805aba01151d762328b9
into topwater/reference. Root only advised after exact-head/test review; did not
merge. The user subsequently merged PR7 at2026-10-26T00:36:44Z as
25b58bda96ab0272cbfb4b0755539d5bb8ca227e. Current topwater/reference contains
handle+display fixes; only PR8 remains open among those three. Topwater's pin995
stays held. Issue1/PR8 have also been rewritten/read-back verified in the agreed
issue/PR format, preserving head805bb29 and stating that prepared loads replace
its zero result. No upstream action or root merge occurred.

Canonical attempt4 is TERMINAL SUCCESS: public60s WEC-Sim solve2393.743s,
public render226.932s,3001 motion samples,751 surface frames,751 decoded H264
frames at800x800/12.5fps/60.08s. Evidence:
artifacts/workflow/canonical-pipeline-attempt4/{report.json,terminal.json,output.log}.
Frozen source/input hashes and the exact recorded process tree were verified
terminal by independent CIM; idle was read back and the lease released.
Do NOT re-poll session61734 or restart this expensive run.

Root viewed decoded0/30/51.44/60s frames: full hull/water coverage and readable
HUD. Basic implementation/visual delivery passes; physics accuracy does not.
Actual total-minus-incident RMS:5.987mm full exterior/2.242% incident RMS;
7.993mm one-beam/2.960%; resolved peak65.087mm at44.56s. Full warning triage
is docs/cli/warning-triage.md. Its remaining correctness questions own issues.

User added automatic solve output result.png, then shared Y scales: translations
share one symmetric metre range, rotations a separate symmetric degree range.
New result_plot.py +solve_outputs.py publish NPZ/PNG after private preparation,
with rollback on handled failures. Code/CLI tests286PASS1optional skip11.23s;
latest shared-axis source checks44PASS5.07s. All jobs terminal/frozen/process
verified and leases released. The PNG was generated FROM the completed Result,
without rerunning physics: result-plot-backfill-attempt2/report.json,2160x1170,
all3001 samples, native-q agreement2.63e-14. Original seven data/video/log files
are BITEXACT. Root viewed final shared-axis PNG. Original result.log truthfully
records the earlier solve; it was not retroactively edited to fake PNG creation.

Public generator replay matches every sea array; unavailable-native check in a
PRIVATE copy preserved all completed canonical files/logs. Public-replay-attempt1
is terminal/released. Runtime is much slower than root's15-20minute estimate;
report measured times, not real-time claims. JIT experiment fully reverted.

No numerical job is live. Checks61218 and inertia export30729 are CLOSED/released.
Verify ledger/lease and exact processes before any additional numerical claim.
Existing shared MATLAB CAPYTAINE_REFERENCE engine29796/launcher40036 was not
quit or restored onscreen. Keep its window minimized/headless. TEMP/TMP must
be normal AppData/Local/Temp before engine discovery/connect.

## Public contract

One flat examples/simple-boat/ contains exactly eight final files:
boat.npz, sea.npz, result.npz, result.png, video.mp4, sea.log, result.log, video.log.
Future generate-boat/boat.log is not implemented now. Commands:

```powershell
uv run topwater generate-sea [DIRECTORY] --record RECORD_ID [--seed N]
uv run topwater solve [DIRECTORY] DURATION_SECONDS [--solver=native|wecsim]
uv run topwater render [DIRECTORY]
```

Directory defaults to cwd; filenames are fixed. Native default fails clearly;
independent native coupled dynamics remain excluded. No solve/render streaming
or pipelining feature is authorized. Current WEC integration dt=.02s; surface
output targets .08s (3001 motion samples/751 frames for60s, endpoints included).
Future integration may be adaptive independently of saved output cadence.

Preserve the accepted24ft/9ft/2500kg physical boat and measured stronger sea
ndbc-46086-20231223T0040Z, seed42. Actual generate-sea replay reproduced EVERY
stored sea array bit-exactly and generated canonical sea.log.

Public NPZ files have numeric arrays only, exact keys, allow_pickle=False, no
metadata/version/solver marker or incident-elevation cache. Boat is an uncapped
convex physical triangular hull with exactly one inferred opening. Internal caps
carry no mass/pressure and never define a waterline. Result includes its physical
hull snapshot/CG, absolute poses and mandatory packed total-water geometry with
independent surface times and frame-local connectivity/offsets.

LATEST APPROVED GEOMETRY: preserve the producer water mesh; DO NOT cut/remesh
triangles on export. Hull-interior portions are invalid fluid and implicitly
excluded. A separate shared surface_geometry routine computes the physical
hull/total-surface intersection. WEC output remains incident+diffraction+radiation
from its already-computed motion; physics/forcing/motion are unchanged. The
format also permits the native solver's own evolving mesh without a linear
field decomposition. Current top-down render has its own documented limits.

Delta metrics integrate exact represented-P1 moments over exterior support using
temporary integration polygons, not an exported contact mesh. One-beam selection
uses ORIGINAL source-triangle centroid distance to the independent waterline;
this is an explicitly discrete band approximation. Exact3D wet-set diagnostics,
startup gates and physical opening/cap checks remain. Native nonlinear output
must not inherit fixed-mean linear applicability criteria indiscriminately.
See docs/npz-format.md, cli/*.md and experts/shared-contract.md.

## Latest verified checkpoint

Pushed Topwater **bda901d**: eight-file example, PNG integration, logs and strategy.
LFS upload342MB verified; Git LFS fsck passes and the development branch matches
origin. Prior source checkpoint537cf08; runtime source SHA:
_capytaine_surface.py2858b847efea11a97c50bd9d12c744bf3069723e32e5fa401309c4e2ee0f352a.

- Complete workflow checks:270PASS/1optional libigl skip,14.11s;
  artifacts/workflow/uncut-workflow-tests-attempt2.
- Shared contact+geometry:34PASS; logical log/CLI/render integration59PASS,
  including actual two-frame encode/decode. Rich Windows flush fragmentation is
  fixed, console bytes/handler configuration preserved. Video summary timings
  are rounded; source and encoded durations remain distinct.
- Three-frame replay0/277/600: original world-grid vertices/connectivity BITEXACT,
  retained arrays unchanged, physics/remesh calls forbidden. Full diagnostics
  took4.2946s; opening+waterline alone.1436s. These are sparse smoke metrics,
  NOT canonical60s disturbance statistics or proof of real-time throughput.
- Paired diagnostic shortcuts:3.902s before/2.509s after (1.55x), packed output
  BITEXACT, wet-value differenceZERO, metrics within1e-12;
  artifacts/workflow/uncut-diagnostics-compare-attempt1. Independent source review
  approved gates/ownership/area rules. All those jobs terminal/frozen/process
  verified, leases released before current public run.

No author is currently editing runtime. Root owns integration/execution. Experts
may research/write non-runtime docs only. Native physics, forcing changes and
new numerical campaigns remain excluded while the public run completes.

## Correctness issues and parallel research

User explicitly requested independent issues and WEC-Sim-only first-line scope:
#31 causal radiation-surface phase; #32 reference mesh/irregular frequencies;
#33 reconstructed field/prepared-body loads; #34 fixed-mean linear applicability.
All four created/updated and first-line readback verified. #30 now links to them
rather than owning their investigation checklists. Issue writes were explicitly
authorized; no upstream messages/PRs/merges were authorized.

Numerical researcher completed READ-ONLY #31 investigation:
experts/numerical_methods_researcher/causal-phase-issue31-investigation-2026-10-25.md.
No sign/normalization fault was proved. With the same band/quadrature, cosine
inversion folds the raw negative-lag part into positive lags. Sharp bandlimits
can produce negative ringing even for causal data; the77mm future contribution
is not proof of noncausal BEM physics. Proposed later discriminator: ONE Heave
near/far saved-transfer check plus known-causal analytic control, separating band
from lag-window effects. No calculation, import, runtime edit or job was run.
This proposal is NOT a new numerical allocation during the current lease.

Keep provisional qualifications: complex near-field recovery not passed;
42.7mm max/12.5mm RMS signed-lag-versus-causal sensitivity is NOT a true-error
bound. Mesh/irregular-frequency risk, field/load interpolation discrepancy and
restoring/wet-set warnings remain unresolved. Successful files/tests do not
certify physical accuracy. Final actual logs still need triage after delivery.

## WEC-Sim fork strategy request

User was not aware how fork changes had progressed. Root disclosed five commits
(six corrections) beyond official v7.1.0: forcing history/directional sum/once-only
dispatch, File Solid inertia mapping, BEMIO handles, literal path display.
The forcing/inertia categories are consequential model changes with SELECTED
path evidence, not universal WEC-Sim validation. Resource/display fixes do not
change numerical parsing or dynamics. Upstream remains untouched.

Current clean pinned fork: amosborne/WEC-Sim, topwater/reference-v7.1.0,
9952634f6f2e2746a1626484e27c3003d1e1d5a9; basea23d299(v7.1.0).
Keep the runtime pin unchanged under the user HOLD; future fork work uses an
isolated checkout. Causal reconstruction is in
Topwater's adapter, not these fork changes.

User now explicitly asks WEC-Sim expert to define a maintenance strategy and
assess possible rework. /root/wecsim_surface_examples completed a bounded
read-only assignment; owns docs/wecsim-fork-strategy.md (and concise KB state).
Classify every logical patch as generic upstream candidate/Topwater integration/
uncertain; check CURRENT upstream branch/release/issues/PRs; propose independent
repros, minimal submission patches and tracking/rebase policy. Pay attention to
R2020b->R2025b SLX migration noise/compatibility and precomputed-force assumptions.
NO fetch/checkout/rebase/cherry-pick/fork edits, numerical work or upstream writes.
Do not infer permission to submit PRs or rework the frozen fork. Plan/review first.

## Current bounded fork rebuild assignment

The first two fixes are published as separate internal PRs in OUR fork:
https://github.com/amosborne/WEC-Sim/pull/6 (handles, issue4,
0d364b1f6154c5b8bee99d95f4fbf7e62b024159) and
https://github.com/amosborne/WEC-Sim/pull/7 (literal paths, issue5,
206bac360222b86ebe33a903225f0dae118c7093). Both target topwater/reference and
were independently based on exact clean releasea23d299. User merged PR6 only;
no further root merge is authorized.
Isolated checkout: artifacts/dependencies/wecsim-rebuild-v7.1.0. Old checkout
and Topwater runtime pin9952634 remain unchanged; combined reference is now25b58bda.

Root serial stock/fixed checks passed on Windows/MATLABR2025b: all9 fixed pass;
stock7 expected fail/2 controls pass. All32 numeric fields EACH sphere/RM3 have
exactly equal values (not byte identity). Independent source/harness review
completed. Attempt1 wrapper and attempt2 helper-discovery failures are preserved;
attempt3 is terminalSUCCESS, exact processes/hashes verified, idle read back and
lease released. Only README/messages changed after successful tested blobs.
Evidence/limits: experts/wecsim_expert/fork-rebuild-regressions-2026-10-25.md.
No full WEC suite, other MATLAB/OS or dynamics simulation was run for these PRs.
Issue1 history/lag correction is now internal PR8:
https://github.com/amosborne/WEC-Sim/pull/8, branch fix/imported-excitation-history,
805bb29b1aa5479988b61d5187f52809656cc82c independently atop a23d299. Source and
independent scalar oracle were reviewed before root's serial execution. Windows/
R2025b history-regressions-attempt1 passed all11 candidate tests; stock10fail/
1regularCIC control pass, no incomplete tests. Selected stock witnesses are real
preprocessing numerical crop/lag failures. All exact processes/hashes were verified
terminal, idle read back and lease released; session78983 is CLOSED. Only README/
message changed after tested source15a38c68. New test/helper/wrapper static clean.

Its generic policy retains complete FILE history, zero-extends beyond it with a
named warning, rejects invalid internal data, preserves the original lag grid and
applies the legacy absolute-time SOURCE ramp once. Existing already-ramped direct
calls stay valid. No new negative-time clamp or mandatory calm-start/padding gate.
Evidence/limits: experts/wecsim_expert/imported-history-rebuild-2026-10-25.md.
No dynamics, multi-import, full-suite, MATLAB/OS matrix or performance qualification.
Inertia remains outstanding on the active rigid path; multi-import work is out
of current scope. No upstream write; user merged PR6, combined branch76adb292,
Topwater/oldcheckout995 and public artifacts unchanged. User later merged PR7;
current combined reference25b58bda includes both resource/display corrections.
Pin promotion and legacy branch retirement remain on hold.

The broader issue3 source assignment is STOPPED by latest scope steering above.
Preserve uncommitted tests/fileSolidInertiaOrderTest.m and tests/helpers/ in the
isolated fix/file-solid-inertia-order branch. No generated SLX was promoted.
Root's inertia-export-attempt1 supported API/save/export/reopen succeeded as
R2020b, but strict audit FAILED with53 unexplained entries. Candidate SHA:
82ca792005f4969ce578a48a0517b7633ff509d901e35db7928435115e1059f6, private
under that job's private-library/exported/. Exact processes and frozen identities
were verified terminal; idle read back and lease released. No native tensor
probes ran. Do not relax the audit or resume export/native development now.

New NARROW source assignment after d043f8d: WEC expert owns
artifacts/workflow/rigid-inertia-probe-v1/ for the ACTIVE Rigid Hydro Solid only.
Same inputs before/after, six independent wrench columns with direct initial
acceleration, identical dry/production-adjusted bodies in two modes; six columns
x two modes means twelve tiny native runs when later approved. Infer the native
6x6 from measurements; label added residual composition separately as algebra.
A single roll-torque option gives the smallest issue3 example. No model load,
MATLAB call, test/job or source-library edit is allocated to author. Preserve
the interrupted broad fixture. Root must obtain independent source review before
any serial probe. No Drag/Flex, new SLX export, branch promotion or pin change.

Completed read-only issue2 proposal is retained as OPTIONAL generic work:
experts/wecsim_expert/multiimport-rebuild-proposal-2026-10-25.md and own-fork
issue2 comment6405172942. No implementation or numerical job is allocated.
Its needed future gate is simulated downstream native sum/bus evidence; old
receipt was graph-only. Current priority is the actual-codepath/minimum-patch
read-only audit, not completing all five generic fork issues.
## Current fork process and next actions

User requested unwind/rebuild through own-fork issues, focused branches and
internal PRs. Issues-first is COMPLETE: amosborne/WEC-Sim issues1-5 cover history,
aggregation/dispatch, inertia/SLX, handles and literal display. Expert enabled
ONLY that fork's issues and verified every body. Subsequent isolated fixes and
internal PRs are recorded above; no upstream writes occurred. Root CREATED
clean topwater/reference at official latest releasev7.1.0/a23d299. The user has
since merged handle fix PR6 as76adb292; it is not yet Topwater's runtime pin.
Future internal PR base repo MUST be amosborne/WEC-Sim, branch topwater/reference.

The user explicitly says DO NOT UPDATE TOPWATER YET. Keep9952634/current local
checkout. Ultimately delete old topwater/reference-v7.1.0, but only after the
new reference is reviewed/tested, Topwater repinned and old exact SHA retained
under an immutable archive tag. Do not leave the old branch as a second active
baseline. Release upgrades require a good reason and separately tested candidate;
no automatic main/dev movement. No upstream issue or PR without specific explicit
user permission. Current strategy: docs/wecsim-fork-strategy.md, independently
reviewed; generic candidates are NOT automatically submission-ready.

1. Working workflow/artifacts/source/docs are committed and pushed as bda901d.
   Eight files, LFS objects and final shared-axis PNG verified. Authorized issue30
   now records completion, actual times and independent correctness/fork links.
2. Keep issue31's next stored-data/analytic-control experiment queued, not
   automatically allocated; user prioritizes fork organization and working delivery.
3. Continue fork rebuild only through bounded expert assignments and reviewed
   own-fork issues/PRs. User's immediate steps were issues-first and clean reference
   branch; those are done. Never silently replay coupled conditional changes,
   switch the held Topwater pin, merge without authorization/checks, or contact
   upstream. The required future branch retirement remains outstanding.
4. The complete initial workflow has implementation evidence; the broader active
   goal includes the newly requested fork reorganization. Do not declare that
   future reorganization complete just because the example now runs.

Use artifacts/workflow/claim_check.ps1, run_check.py and finish_check.ps1 for
serial jobs. Claim, separate readback, run, verify terminal/frozen/process state,
record idle under lease, separate readback, then release lease as FINAL action.
PowerShell scriptblock invocation avoids local script execution-policy issues.
Native DLL/FFmpeg/MATLAB jobs and authoritative CIM need unsandboxed execution.

Historical failures and original policies remain in their evidence and Git.
Attempt3 was intentionally stopped by user after320.329s; incomplete result.log
has no END. Attempt2 retained the complete immutable field cache at
C:/Users/alex/AppData/Local/Temp/topwater-wecsim-mt5icmiq/surface-recovery/completed-field.npz,
SHA66039a8f78cfb8f9dcdec5d565e2197e9a06b01361521b19e0d5d46279613549.
Do not alter/delete it. Earlier clipped export passed751 frames in28.5min but is
superseded by the uncut contract; do not rerun it. Prior handoff chronology is
preserved at commit537cf08 and older immutable checkpoints; no study/policy was
deleted by condensing this handoff.
