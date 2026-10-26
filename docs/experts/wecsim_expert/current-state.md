# Current state: WEC-Sim / BEMIO expert

Updated2026-10-26 UTC. Currentauthority: [handoff](../../checkpoint-handoff.md),
[physicalcontract](../shared-contract.md), [serialprocedure](../coordination.md).
No author execution allocation or automatic follow-up. Root owns integration.

Latest integration checkpoint: [PR #10](https://github.com/amosborne/WEC-Sim/pull/10)
is a draft at `eb722540`, explicitly depending on PR #8. Exact-head four
group/native-bus and eleven history controls pass; retained preprocessing error
is 18 before and 0 after. The strict one-header archive and fresh API reopen
passed. Integration with PR #8/#9 and whole-boat force/startup/motion remain open.

The fixed T320 coarse/refined coefficient controls passed all unchanged finite-
probe gates at three genuine headings plus seam. The subsequent full-band
25-heading Capytaine candidate completed 3606 genuine solves and preserved every
old numerical quantity exactly. Non-author review accepted coefficient handoff
only. Mesh-resolution (3180) and irregular-frequency (2040) warnings remain
unqualified. Official actual H5 readback source is being prepared; no production
settings/pin/public physics promotion occurred. Current evidence and limits:
[qualification checkpoint](kernel-and-inertia-checkpoint-2026-10-26.md).

Active-rigid inertia [PR9](https://github.com/amosborne/WEC-Sim/pull/9) is ready
for review at pushed head `67e7c93f7039881a4af9b6f44a4f6232c029429c`, branch
`fix/rigid-hydro-inertia-order`, based on combined25b58b. The corrected R2020b
library and portable reproduction helper are included. Strict extraction audits
passed control0/candidate1 with no other changes; generic helper reproduced the
same candidate bytes; fresh absolute-path API reopen passed without a save.
Exact-head Windows/R2025b checks passed8layout and3native tests (12simulations),
plus the public two-simulation minimal example. Cap independently accepted the
measured evidence; processes/frozen identities were verified terminal and leases
released. New full-matrix claim is the tested1e-7 tolerance, not historical
roundoff values: that run's detailed matrix report was temporary. Retained
minimal-example raw data and all acceptance paths are in
[PR9 qualification](inertia-pr9-qualification-2026-10-26.md).
The tested README predates validation; the PR body carries the newer status.
Do not amend this tested head merely to update prose. No merge, dependency repin
or upstream write occurred; held995 and public physics artifacts are unchanged.

Root's native-excitation component checks passed: backend plus offline comparison
36/36, and MATLAB startup diagnostic/gate2/2. After the reviewed test-only cwd
teardown correction (`01f60c0c`), startup attempt2 again passed2/2 with a clean
MATLAB log and no warnings. Exact frozen sources and process termination were
verified; leases released. The first attempt's Windows deletion warnings remain
preserved; production math did not change. These component checks do not qualify
actual-case forces or new fork dispatch. Evidence:
`artifacts/workflow/native-startup-tests-attempt2/`.

User authorized WEC-Sim-owned excitation and motion. The unexecuted integration
source now emits physical directional origin-elevation records, preserving both
signed branches and padded absolute time; no Python force is supplied or rebound.
MATLAB performs native initial-force certification and applied-bus identity.
The agreed H5 handoff retains the exact original full-band radiation kernel/Ainf,
with genuine dense raw excitation and its explicit low-frequency bridge.
Offline force comparison remains separate from production. Isolated candidate
`fix/native-imported-excitation` is based on combined `25b58b`, then explicit
history dependency `805bb29` replayed as `0335022`; the separate group/rigid
dispatch patch and supported API preparation sources await review/execution.
Held runtime995, public artifacts and published topic branches are unchanged.
See `artifacts/workflow/native-excitation-comparison-v1/README.md` for gates.

Earlier standalone active rigid Hydro Simscape inertia evidence PASSES the measured single-
torque witness and twelve six-basis before/after runs. Corrected native dry
matrix error is5.33e-15; native adjusted matrix plus unchanged algebraic residual
matches M+A to3.55e-15. Eight layout controls pass. Independent review accepted
these narrow claims; exact processes/frozen files were verified terminal and
leases released. This does not test dynamic residual/CIC forces or boat accuracy.
See [measured proof](inertia-adjustment-proof-2026-10-25.md). Failed earlier
expectation/link trials remain preserved; no further inertia job is authorized.
The broader interrupted Drag/Flex/export sources stay preserved and stopped.
That historical witness did not promote a library or public artifact; the later
focused PR9 binary qualification is recorded above.

Issues-first work is EXECUTED: enabled issues only in `amosborne/WEC-Sim`,
checked duplicates, and created/read-back verified fork issues
[#1](https://github.com/amosborne/WEC-Sim/issues/1),
[#2](https://github.com/amosborne/WEC-Sim/issues/2),
[#3](https://github.com/amosborne/WEC-Sim/issues/3),
[#4](https://github.com/amosborne/WEC-Sim/issues/4), and
[#5](https://github.com/amosborne/WEC-Sim/issues/5).
The user squash-merged handle fix PR #6 into `topwater/reference` as
`76adb29201291ffdec3e805aba01151d762328b9`, then display PR #7 as
`25b58bda96ab0272cbfb4b0755539d5bb8ca227e`, based on release `a23d299`.
The old branch/runtime pin `9952634` remains unchanged and on USER HOLD. The old branch is to be
deleted only after its reviewed/tested replacement, exact repin and immutable
archive tag. [Strategy and scope](../../wecsim-fork-strategy.md). No new fork
code/topic branch/PR or upstream write is authorized by this KB update.

Bounded resource/display rebuild is locally tested in isolated
`artifacts/dependencies/wecsim-rebuild-v7.1.0`, origin explicitly our fork.
Issue #4 branch `fix/capytaine-sidecar-handles` is local commit
`0d364b1f6154c5b8bee99d95f4fbf7e62b024159`; issue #5 branch
`fix/literal-wave-path-display` is `206bac360222b86ebe33a903225f0dae118c7093`.
Each topic was independently based on exact `a23d299`; old clean runtime
`9952634` is unchanged. Windows/R2025b attempt 3 passed
all nine fixed tests; stock failed seven and passed two controls. All 32 numeric
fields for EACH sphere/RM3 fixture have exactly equal stock/fixed values by
`assert_array_equal`, not byte identity or a dtype check. Fresh processes used
repository fixtures/direct display, without Topwater cases or body simulation.
Exact MATLAB processes were terminal/CIM/frozen-identity verified; lease released.
[Evidence](../../../artifacts/workflow/fork-regressions-attempt3/report.json).
Only candidate README validation and commit messages were then amended; all four
production/test Git blobs match the executed `54bca31`/`e8491f4` candidates.
Broader OS/MATLAB compatibility remains untested. Root published both branches
and verified internal [PR #6](https://github.com/amosborne/WEC-Sim/pull/6) and
[PR #7](https://github.com/amosborne/WEC-Sim/pull/7) with these exact heads,
base `topwater/reference`, and no cross-repository target. The user subsequently
merged PRs #6/#7; current reference is `25b58bda`. Topwater's pin remains on user hold. CI does not target
this reference branch; local checks above are the validation evidence.

## Generic CLI integration checkpoint

Current-path lifecycle corrections passed independent source review and root's
serial Windows/R2025b checks (current-path-checks-attempt2):
combined native/cleanup errors retain the primary transcript plus
`EngineCleanupError`/recovery semantics; the guard snapshots/restores and checks
root `ErrorIfLoadNewModel`. Focused Python log/recovery and MATLAB setting
success/failure witnesses reproduce the old failures. Corrected Python backend
suite21PASS; MATLAB13PASS (2setting +11finite hydrodynamics). Exact processes and
frozen inputs were independently verified terminal; lease released. The finite-
hydrodynamics helper is included in backend resource copying. No physics,
fork pin, model or public artifacts changed.
[Scope and executable pre/post examples](current-path-lifecycle-corrections-2026-10-25.md).

Issue #1 history/lag-only candidate is locally verified and frozen at
`fix/imported-excitation-history`, commit
`805bb29b1aa5479988b61d5187f52809656cc82c`, directly based on `a23d299`.
It retains the full first imported record, original lag sample grid and literal
legacy source-time ramp; uses zero extension outside FILE bounds with one named
support warning; preserves already-ramped legacy calls. No calm-history clamp,
mandatory padding, aggregation or dispatch changes. Independent review approved
the executed source; new test/helper/wrapper static checks were clean.
Windows/R2025b [attempt 1](../../../artifacts/workflow/history-regressions-attempt1/report.json)
passed all eleven fixed tests; stock failed ten and passed the regular-wave
control, with no incomplete tests. Real-preprocessing crop/lag numerical
witnesses failed cleanly on stock. Exact processes/CIM/frozen identities were
verified and the lease released. Only README validation and commit message were
amended; production/test/helper blobs match executed `15a38c68` exactly.
Scope is excitation-only: no dynamics, multi-import, full-suite, MATLAB/OS
matrix or performance qualification. Other topics, combined base, held runtime
pin and public artifacts are unchanged. Root published the topic as internal
[PR #8](https://github.com/amosborne/WEC-Sim/pull/8), targeting our
`topwater/reference`. No merge or upstream write occurred.

The cosmetic `waveClass/listInfo` path fix passed a real three-path display
probe: verbatim Windows/percent/UNC paths, zero warnings, no wave setup/body run
or new engine. Helper style is clean; 113 upstream findings are unchanged.
Fork `9952634f6f2e2746a1626484e27c3003d1e1d5a9` is pushed and pinned. No wave
physics changes. See [verified evidence](wave-path-display-2026-10-25.md).
Canonical public attempt 4 solve AND independent render succeeded: 3001 motion
samples and all 751 video frames decoded. Solve wall time 2393.743 s; render
226.932 s. Frozen identities and exact process terminal status were verified by
CIM; lease released. See
[run evidence](../../../artifacts/workflow/canonical-pipeline-attempt4/).
The working end-to-end result remains a provisional linear reference, not a
physics/mesh/time-accuracy qualification.

Current outputs are [result.npz](../../../examples/simple-boat/result.npz),
[video.mp4](../../../examples/simple-boat/video.mp4) and
[result.png](../../../examples/simple-boat/result.png), with command logs beside
them. PNG was produced separately from completed Result data without resimulation.
Full suite passed 286 tests/one optional skip in 11.23 s; latest shared-Y plot
style checks passed 44 tests in 5.07 s, with actual PNG/native-q difference
2.63e-14. Canonical sea generation remains array-bit-exact. Artifact checkpoint
bda901d is committed/pushed; the342MB LFS upload and object integrity are verified.
This note does not authorize another run.

The first full canonical CLI attempt ran and failed during surface export after
1989.285 s; no result/video was published. A leaked BEMIO KH probe handle caused
Windows scratch removal to mask the original export exception. Root confirmed
and closed only the owned MATLAB file ID. Maintained-importer scoped closure and
package primary-error preservation passed the allocated regression: old leak
reproduced, two new imports preserved 32 numeric arrays each, malformed KH and
Hydrostatics left no handles, and 15 Windows fixture files renamed/deleted
(reported 1.672 s). Importer correction commit
`f07a8dfe1d81d9430a22c299a83d91252afb3a33` is included in the current fork pin.
Recovery/interfaces passed 209 tests with one optional skip in 13.05 s.
The second full canonical attempt terminated after 1397.677 s without a public
result: frame 277 at 22.16 s generated a near-collinear clipped triangle, not
demonstrated flooding. Its completed field/failure geometry are retained for
bounded replay; no additional body solve is needed for geometry diagnosis. See
[failure/lifecycle correction](bemio-handle-lifecycle-2026-10-25.md).

WEC-WORKFLOW-001 added `src/topwater/wecsim_backend.py`, shared
`_wecsim_setup.py`, packaged MATLAB resources under `data/wecsim`, and
`docs/cli/solve.md`. The source supplies generic full equilibrium, fixed
duration-independent incoming train and full-hull startup assessment, optional
dependency/fork/shared-engine readiness, native model construction and bus
audits, output-pose conversion, and calm-restoring warning aggregation. New
test source covers equilibrium, fail-before-physics readiness, pose embedding,
warning aggregation, and cleanup exception paths. Root's released lifecycle
checkpoint passed 13 backend tests, all four MATLAB style checks, and two real
normal/forced-failure cycles in the existing hidden shared engine. Exact
cwd/path, diary, warning/thread/figure settings and cache/codegen folders were
restored; only the unique owned model was closed. Distinct temporary class
directories were removed between cycles. The later first canonical run completed
native body integration but did not publish a complete surface result. See
[versioned lifecycle evidence](workflow-engine-lifecycle-2026-10-25.md).

The coordinator accepted the selected-case provisional total-surface policy:
24 rad/s radiation cutoff, spacing <=0.025 rad/s, full native history, zero direct
terms, exact packet sideband diffraction, and explicit exterior-trace continuation.
The Capytaine production API matches the backend's four packed surface arrays.
Source includes the represented-total wet-set 10% warning and initial local
total-elevation check. Current output preserves the unremeshed full field grid;
hull-interior fluid is implicitly invalid, with a separate physical-waterline
routine and temporary exterior P1 integration polygons for diagnostics.
Public workflow completion passed; full-complex near-field recovery and total
disturbed-fluid velocity bounds remain unqualified.
The 42.7 mm maximum/12.5 mm RMS signed-lag comparison is a numerical diagnostic,
not an error bound. See [method review](../numerical_methods_researcher/surface-reconstruction-review-2026-10-25.md)
and [solve contract](../../cli/solve.md). The allocated test suite passed actual
example geometry/M/K/CG identity against the frozen selected 17-panel case.
The old incident-only surface is never promoted to total surface.
Native Topwater dynamics remain outside this assignment.

## Historical completed incident-field video

[Full60s topdown video](../../../artifacts/reference/stronger-buoy-20261024/matplotlib-full60-attempt1.mp4)
SHA15510982f6cef3a246f1129222e0525052215a06cfa2e0073af3c53bf803f23b,
H264800square/751frames12.5fps/60.08s, native0..60retained. All751decode and
owner0/30/51.44/60stillsPASS; ROOT0/peakPASS, USERlinked. Exactterminal/frozen
inputs/ownerleaseFINALbacd1e. [Finaldata/video evidence](full60-scene-video-2026-10-25.md).

Original751hull/time/central65x65waterXYZ BITEXACT, nativeqHUD exact0:4:3000.
ONLY1550same-train borderpoints/frame added to75x77, allfixedviewcorners covered;
no96m expansion or sea/physics changes. Reusablefullscene884667cb (~104MBwater).
Prep40.657s/CLIwall148.792s/libraryrender145.75s/wholejob197.911s. Generic2Dstyle
and renderer unchanged duringfulljob; no bespokeWECdrawing.

Historical workflow `uv run topwater render INPUT.npz OUTPUT.mp4`: two requiredpaths,
wholeNPZ/fixed800/savedcadence/MP4only. GenericMatplotliblibrary/CLI owned by
visualization; WECownsdataadapter/nativeinterface. Active3Drender-only code was
USERapproveddeleted afterimportaudit; physical3Dgeometry/hydro/WECfork untouched.
Historicalmovies/frozenarchives/failures retained. No furthercleanup/jobs queued.

## Physical reference and practical limits

Actual46086/2023-12-23T00:40Z Hm01.65747/Tm027.23657, seed42/unscaled384carriers,
physicaltrain tc40/plateau15/support25; numerical-local-calm/restv0/memory0/ramp0/
viscous0. Originalfixed17H5/CIC10/dt.02/all6free. Native60s implementationPASS:
roll/pitch/yaw6.618565/2.677850/1.325248deg, translations.562/.837/.759m.
Net radiationwork32971.040J/temporaryreturn4074.896J/cumulativemin0;
incident46582.193J, not completefluidenergyclosure/generalpassivity proof.
Cachedrelativewaterplane maxheel3.10deg/trim2.88deg/heave.09m is a proxy,
not nonlinear/wetting/flooding/mesh/timeaccuracy qualification. Topwater physics
implementation remains HELD; reference visualization never substitutes it.

## Versioned evidence

- [Strongerphysicalcase](stronger-buoy-readiness-2026-10-24.md).
- [Telemetry/dataidentity and5s proof](telemetry-2d-readiness-2026-10-24.md).
- [SharedCLI/PyVista history](pyvista-readiness-2026-10-24.md).
- [Previouspresentation history](presentation-readiness-2026-10-24.md).
- [Officialprovidedbaselines](official-examples-readiness-2026-10-24.md).
- [Mildarrival](calm-buoy-arrival-readiness-2026-10-24.md),
  [mildassessment/body](calm-buoy-assessment-2026-10-24.md).

Historical3Dnotes describe archivedcode/results, not currentlysupportedrenderers.
