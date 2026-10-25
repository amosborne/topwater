# Checkpoint handoff

Updated 2026-10-25. Latest user steering controls scope. Read AGENTS.md and
experts/coordination.md; consult shared-contract.md and the relevant expert state.
Historical queued work is not authorization to restart it.

## Current priority and completed workflow

The user's order is working end to end, then correctness, then speed. The full
public workflow now WORKS. No further performance changes are authorized now.
Native coupled dynamics remain unimplemented. No solver/render streaming.

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

No numerical job is currently live; verify ledger/lease before another claim.
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

## Current fork process and next actions

User requested unwind/rebuild through own-fork issues, focused branches and
internal PRs. Issues-first is COMPLETE: amosborne/WEC-Sim issues1-5 cover history,
aggregation/dispatch, inertia/SLX, handles and literal display. Expert enabled
ONLY that fork's issues and verified every body. No topic branches, PRs, source
edits or upstream writes. Root CREATED clean topwater/reference at official
latest releasev7.1.0/a23d299. This new branch contains NO local fixes yet.
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
