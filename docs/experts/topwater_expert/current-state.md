# Current state: Topwater expert (solver implementation HELD)

LATEST USER authorizes a standard reusable visualization workflow in the
Topwater library, with producer-neutral prepared meshes, different convex hulls,
code-first CLI, square water-filled fast previews and no custom graphics engine.
Initial package integration/source is owned by `/root/visualization`; WEC owner
owns only reference-data preparation/adapter. Dedicated visualization role is a
new routing entry awaiting native-client availability; no activation claim here.

COMPLETE: `src/topwater/visualization.py`, the Matplotlib/geometry modules and
`cli.py`. Clean command `uv run topwater render INPUT.npz OUTPUT.mp4`: both paths
required; fixed 800 square, whole input, uniform native cadence, only MP4 persists,
diagnostics stdout. Matplotlib and bundled FFmpeg are core dependencies; PyVista
remains optional. `uv.lock` is current with existing package versions preserved.

LATEST ACCEPTED VIEW supersedes the earlier oblique presentation: Matplotlib
top-down, continuous blue-white-red elevation fixed about zero, no water grid,
light-gray projected hull edges, black mesh-boundary rim and blue actual physical
hull/water triangle intersection. Header time/roll/pitch/heave comes from named
producer telemetry, never inferred from a reference hull. Clean command:
`uv run topwater render data/scene.npz videos/preview.mp4` (no extra required).

Latest accepted `matplotlib-clip-attempt3.mp4` SHA
7134d6443cc8d66e3e53d2719f63618bf087cb432127f1edb54eb20711dec7b0,
under stronger-buoy-20261024. Large14.5pt DejaVuSansMono fixed-column header,
opaque generic projected off-white hull footprint, muted continuous coolwarm,
dark-blue contact. 17 baseline tests plus2new header/footprint checks passed;
font-only finaliteration needed no repeated tests. All63 frames decode at800
square/12.5fps/5.04s; ROOT and owner viewed the decoded peak still.
Actual full public CLI14.027s including process startup/imports/load/flush;
renderer12.25s, geometry/artist4.132s, Aggdraw4.773s, pipe-send.048s,
encoderflush2.500s. Pipe-send is transfer/backpressure, not pure encodeCPU;
encoding overlaps drawing. Figure/artists reused, blitted Agg buffer streams
directly through public imageio-ffmpeg without a second full-figure draw.
Visible-water culling draws8064/107648triangles. Geometry/time/connectivity are
bit exact; telemetry rows match original nativeq exactly. No wave preparation
or physical correction. Exactterminal75481c/CIMe51809, frozenPASS62d995,
idle/separateca3566/FINALreleasefc086d. No jobs remain.

3D code/artifacts are preserved pending USER review of 2D; do not delete based
on that conditional alone. The earlier same-CLI synthetic convex box checked
another topology. Its oblique timing/receipts are preserved in the
[visualization contract/source note](visualization-workflow-2026-10-24.md).
No further rendering work is allocated. Solver and physics remain HELD.

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
