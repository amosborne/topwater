# Current state: Topwater expert (solver implementation HELD)

## Issue31 integration acceptance, 2026-10-27

Full-complex WEC-Sim surface postprocessing is enabled after selected-case
whole-region qualification (5.784425mm combined numerical budget versus10mm).
This does not implement native dynamics or qualify nonlinear/model accuracy.
Final root-executed integration cohort passes50tests in24.20s, no skips;
independent receipt verifies all numeric module ASTs unchanged apart from the
readiness Boolean. Exact cohort, terminal proofs and artifact identities are
in [the tracked acceptance record](../numerical_methods_researcher/issue31-acceptance-evidence-2026-10-27.json).
Canonical result/video/logs are regenerated; all motion and inputs are unchanged.
No numerical job is delegated. Root owns final checkpoint/issue closure and
future measured-speed allocation; #32–34 remain separate.

## Current workflow geometry checkpoint (2026-10-25)

Latest user explicitly approved the simpler uncut-surface contract. Preserve the
same total incident+diffraction+radiation field and its original rectangular
water vertices/face indices. Finite stored values inside the posed hull are
physically invalid/implicitly excluded by the hull; no extra NPZ keys or forced
mesh hole are introduced. Compute the waterline separately. A future native
solver supplies its own surface and may evolve it internally.

Implemented/frozen shared interface: `surface_geometry.waterline_segments`
reuses the existing visualization intersection kernel, with physical hull faces
only (or explicit artificial-cap exclusion). Visualization reexports the same
function. `validate_opening_clearance` exposes the unchanged complete-rim-edge
and internal-cap compatibility checks without clipping/remeshing. Rim/cap
contact remains fatal; physical coplanar contact remains fatal in the waterline
kernel. No wave values, hull geometry, tolerances or physical gates changed.

ROOT reports `tests/test_surface_geometry.py tests/test_surface_export.py`:
34 passed in 0.50 s, frozen identities/process absence verified, exact producer
terminal and lease released. Evidence:
[uncut contact test log](../../../artifacts/workflow/uncut-contact-tests-attempt1/output.log).
The old real-three-frame clipping A/B benchmark was NOT run because that export
workflow is superseded. Shared helper stays frozen for the producer consumer.
Producer adaptation, cached uncut replay and complete example/video acceptance
remain pending; no complete-case or independent native-physics claim follows.

Boat/Sea/Result NPZ contracts and convex single-opening hull ingestion/internal
caps remain implemented. Nonplanar partial-facet cap closure has an explicit
backend limitation; no general CAD repair is claimed. Earlier clipping studies,
including the guarded frame-277 diagonal correction and 751-frame geometry
pass, remain historical evidence under their original policies:
[mesh tooling](../../mesh-tooling.md).

## Earlier visualization checkpoint

LATEST USER authorizes a standard reusable visualization workflow in the
Topwater library, with producer-neutral prepared meshes, different convex hulls,
code-first CLI, square water-filled fast previews and no custom graphics engine.
Initial package integration/source is owned by `/root/visualization`; WEC owner
owns only reference-data preparation/adapter. Dedicated visualization role is a
new routing entry awaiting native-client availability; no activation claim here.

COMPLETE: `src/topwater/visualization.py`, the Matplotlib/geometry modules and
`cli.py`. Clean command `uv run topwater render INPUT.npz OUTPUT.mp4`: both paths
required; fixed 800 square, whole input, uniform native cadence, only MP4 persists,
diagnostics stdout. Matplotlib and bundled FFmpeg are core dependencies.
USER approved removal of the 3D renderer: generic PyVista code, dispatcher,
3D-only settings and optional dependencies have been removed. `uv.lock` removes
19 unused dependencies, preserving the remaining versions. Three-dimensional
mesh inputs, convex validation and physical waterline intersection math remain.
Final post-cleanup suite `tests/test_visualization.py tests/test_package.py`:
19 passed in0.15s on the actual cleaned source; no source changes were needed.
Exact testmanifest1e32810f9952b5334b97ded502500e80368424254a79184fe466a57b2d8a429f,
terminal101701/CIM73aa68, frozenPASS531e70/separatec617a1/release6a9c91.

LATEST ACCEPTED VIEW supersedes the earlier oblique presentation: Matplotlib
top-down, continuous blue-white-red elevation fixed about zero, no water grid,
light-gray projected hull edges, black mesh-boundary rim and blue actual physical
hull/water triangle intersection. Header time/roll/pitch/heave comes from named
producer telemetry, never inferred from a reference hull. Clean command:
`uv run topwater render data/scene.npz videos/preview.mp4` (no extra required).

Latest accepted `matplotlib-full60-attempt1.mp4` SHA
15510982f6cef3a246f1129222e0525052215a06cfa2e0073af3c53bf803f23b,
under stronger-buoy-20261024. Large14.5pt DejaVuSansMono fixed-column header,
opaque generic projected off-white hull footprint, muted continuous coolwarm,
dark-blue contact. 17 baseline tests plus2new header/footprint checks passed;
font-only finaliteration needed no repeated tests. Full751 frames decode at800
square/12.5fps/60.08s from physical0..60s; owner viewed0/30/51.44/60 decoded
images andROOTviewed0+peak. Header/opaquehull/continuouswater style unchanged.
Actual full public CLI148.792s including process startup/imports/load/flush;
renderer145.75s, geometry/artist64.066s, Aggdraw68.732s, pipe-send.530s,
encoderflush1.703s. Pipe-send is transfer/backpressure, not pure encodeCPU;
encoding overlaps drawing. Figure/artists reused, blitted Agg buffer streams
directly through public imageio-ffmpeg without a second full-figure draw.
Input is104587856bytes:75x77 nodes from original65x65cache plus1550borderpoints
perframe evaluated from the same frozen incoming field only. Original waterXYZ,
hull/time BITEXACT; HUD nativeqrows0:4:3000 exact. Preparation40.657s; whole
producer+QA197.911s; no simulation/physics or sea-state changes. Visible water
draws10366/11248triangles. Exactterminal6f2f47/CIM6c6e36, frozenPASSb012d4,
idle/separate58a912/FINALreleasebacd1e. No jobs remain.

Prior 3D sources are preserved in Git history and old artifacts remain intact.
The earlier same-CLI synthetic convex box checked another topology. Its oblique
timing/receipts are historical, referenced by the
[visualization contract/source note](visualization-workflow-2026-10-24.md).
Latest full60s video and approved active3D renderer/dependency cleanup are
complete. Separate read-only speed suggestions are proposals only; no follow-up
render/optimization job is allocated. Solver and physics remain HELD.

LATEST user, after mild60s video: use a MORE energetic actual sea and report
diagnostics. Coordinator selects `ndbc-46086-20231223T0040Z`, source row29984,
Hm01.6574679484m/Tm027.236567149s/observed direction resultant.457661. Complete48
bands and source provenance; no response/RAO/phase fit or amplitude scaling.
Existing seed42/8perband/default bow policy retained; actual yaw pending export.
New stronger controls and fresh-path preparer reuse EXACT unchanged exporter22bcf2c0.
Proposal43faa35a independently reviewed/coordinator allocated: EXPORT EXECUTED
PASS1.093s,384carriers/48varianceclosure .1717m2/replay/support/frame checks.
Actual bow182.067194deg, meanprop2.067194deg/R.471121; storedomega/globalphi/IDs
exactly equal mildseed42, measuredamplitudes/MEMdirectionschange. Exact38880/32504
terminal0, frozenverified/idle/separatereadback/finalrelease done. Output
`artifacts/reference/stronger-buoy-20261024/export-attempt1`.
No train/assessment/body/video job in this assignment.
[Stronger executed receipt](stronger-buoy-export-readiness-2026-10-24.md).
Original mild input/body/video remain preserved. Topwater solver remains HELD.

LATEST user regroup2026-10-24: stop bug-chasing and official-example variants;
deliver one LESS energetic buoy-generated incoming sea, true locally calm boatrest,
zero history/no ramp, full60s and centered video with wider physical field.
Coordinator selected actual `ndbc-46086-20230319T0840Z`, Hm0.5351635264m,
Tm027.42692891s, all48 observed bands .030-.495Hz and complete direction moments.
Existing generator seed42/8perband/default bow toward actual realized arriving
mean, no amplitude/phase fitting or rescaling. Neutral export independently SOURCE
reviewed and coordinator-allocated072f1fb5: EXECUTED PASS1.110s, all384carriers,
variance.0179m2/replay/support/frameclosure. Actual default bow201.308143deg,
meanpropagation21.308143deg/resultant.454012. Exactproducer2568/guard29816 terminal0,
allfrozenverified/idle/separatereadback/finalrelease done; no numerical owner/job.
Output `artifacts/reference/calm-buoy-20261024/export-attempt1`.
[Source/executed receipt](calm-buoy-export-readiness-2026-10-24.md).
Topwater solver/tests/optimization/BEM/dynamics remain HELD. WEC/researcher own
physical upstream train/calm-arrival construction and its derived statistics;
stationary source carriers alone are not a calm initial water field.

Latest user steering supersedes old queued follow-ups: finalcase60s, rest/zero
velocity andzero radiation memory, viscousOFF, physics-based incident-band policy.
Old unpruned/history-seeded case remains evidence, not the requested finalstartup.
FullCap radiation+diffraction/FK+diffraction applies; projectlambda2beam filter is
alreadyenforced. No Cap/WEC radiationbug or newlyqualifiedcutoff established.
See currenthandoff; no source role implies numerical execution permission.
Read-only data review finds42,575complete directional observations in the frozen
marine-gyro source library beyond9packaged examples; maxobservedHm06.03052m/
Tm0212.90863s is metadata, not selectedcase. Whole-band subset must preserve exact
original component IDs/a/omega/theta/phi and use distinct retained-energy metadata.
Legacy full-spectrum loader must not silently accept pruned provenance. No source
edit, new realization/export/field/job performed by that review.

Follow the [handoff](../../checkpoint-handoff.md) for scope. User authorized
existing sea-state generator inputs for ONE new credible higher-energy canonical
case, target roll5-10deg and appreciable other motions. This does not authorize
Topwater solver implementation, tests, simulations or optimization.

Generator API/source/catalog review and one bounded neutral export completed.
The existing generator already supplies actual higher-energy catalog observations;
no generator or solver source changed, no amplitude scaling or phase fitting.
Three internal candidate banks use seed42/eight components per48bands and preserve
each observed energy/MEM spreading; heading offsets0/45/90 are selection metadata,
not a canonical simulation suite. Default bow faces realized mean arriving
direction; global and boat-frame conventions/transforms are explicit.

Output: artifacts/reference/generated-sea-candidates-20261024/manifest.json.
All energy/replay/support/frame/source/output checks pass. One generator job,
0.1564761s,384components perrecord; no BEM/dynamics/pytest/solver work. ExactPID1960
terminal proved, process inventories empty; direct Windows Python ExitCode was
not exposed (null retained), while producer completed all checks and shell0.
Ledger idle was recorded under claim before final release; no job remains owned.

[Generator review and executed evidence](sea-generator-review-2026-10-24.md)
records exact source/config/catalog/dependency identities, true higher-energy
Hm0 values2.270-2.278m, directional means/resultants and immutable output hashes.
Canonical response selection, physical startup, native hydrodynamic angular
coverage and wetted/righting validity remain pending other owners/coordinator.
Generated waves are stationary harmonic realizations, not calm-start packets or
consistent initialized nonlinear fluid states. Do not infer roll from Hm0 alone.

One separately allocated incident-input diagnostic completed after geometry's
terminal/released checkpoint. Same3banks and offsets0/45/90, no generation or
solver calls:0..60s/dt0.02,57projected full-hull/CG points, exact eta/slopes/Hessian
and eta_t, localCG plane departures, initial values and time histories.
[Executed input diagnostics](incident-input-proposal-2026-10-24.md) contain tables
and artifact identities. CG rolllike WATER-PLANE slope-proxy peaks range4.24-11.85deg;
they are not boat angles or response qualification. Sampled plane departures
reach0.265-0.780m; initial surface/verticalvelocity/curvature are nonzero. Setting
boat pose to a tangent plane alone cannot initialize fluid/radiation state.
ExactPID7120/Popen.wait0, shell0,2.2994226s; all field/check/hash gates pass,
stderr empty and producer/launcher inventories terminal. Idle ledger recorded
under claim, separate idle readback passed, then lease released. No further
numerical allocation remains for this owner.

Historical selected solver foundation remains the P0 Rankine-panel core from
4e8f511; hydrostatics/instantaneous added mass do not supply coupled wave-response
qualification. [Held source map](held-source-map.md) is historical evidence,
not permission or current liveness. External reference loads/trajectories never
substitute for independently calculated Topwater physics. Preserve generic math
coverage, other worktrees and all old failed/interrupted studies.

Stopping point reached: return actual generator inputs and evidence; await a
new bounded coordinator assignment. No automatic numerical follow-up.
