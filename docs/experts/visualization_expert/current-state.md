# Visualization expert: current state

## Current delivered plots/video, 2026-10-28

Current result-plot source `0e2418b4` and tests `dabadbee` use cividis with a
LINEAR RMS scale. Original-component bins are 15 degrees by 0.025 Hz and display
`sqrt(0.5 * sum(a_i^2))` in metres. Every zero bin remains unmasked at the dark
palette minimum; no positive bin is omitted, clipped or rescaled. The actual
0.115420427 m maximum rounds outward to 0.12 m using two significant digits,
with nice interior ticks and exact colorbar endpoints. Log coloring was rejected.

The 3150x1170 result.png retains the 12-inch six-motion region and centres a
5.1-inch full-height polar circle in the right column, with its horizontal RMS
colorbar below. No overall title/footer or polar subtitles; x labels are exactly
`Time (s)`. Shared translation/rotation limits retain seven symmetric nice ticks
with zero centred. World coming-from histogram edges are unchanged; the polar
axis rotates the default heading up and marks the actual initial heading with
a dotted line. Only frequency rings remain, with labels pointing down on screen.

The spectrum belongs ONLY in result.png. Production solve passes original Sea
explicitly to plotting/publication; Result contains its original 11 physical
members, with NO sea-carrier arrays. Video rendering remains standalone from
result.npz and has no spectrum or Sea input. Historical calls without Sea retain
the 2160x1170 six-motion-only PNG.

The user approved video preview4 and the full re-render is delivered. Production
renderer `3ce2ebd5` and tests `c0a8214b` passed independent source review. It keeps
the full-bleed 800-square view centred on world origin with the original zoom
and full-hull coverage check. Plain equal-inset edge ticks have bottom/left
numbers, top/right ticks, corner clearance and a smaller inside header. The
`Z (m)`/`0` gauge shares the water's exact coolwarm normalization, rounded outward
to 0.1 m (current +/-1.6 m), world-Y metre spacing and zero row. Gray mesh, black
outline, solid blue moving waterline and fixed dashed reference remain unchanged.

Static graticule/gauge artists are cached once through standard transparent Agg
and FigureImage. Each frame draws water/boat, cached instruments, then telemetry;
no per-frame full figure redraw. New encoder-stub tests check two changing poses,
overlay persistence/transparency/order and unit/zero mapping. The final combined
187-test cohort PASSED in 17.22 s. The [final delivery evidence](../../../artifacts/workflow/approved-layout-delivery1/results/report.json)
confirms all 1626 frames decoded at 800x800, 12.5 fps and 130.08 s encoded.
Ordinary CLI time was 433.312 s; renderer 404.656 s (geometry 199.706 s, drawing
195.175 s); the complete display-only driver took 492.061 s. These are whole-run
timings, not an isolated overlay benchmark. Result `aa050fa5` and its saved motion/
field remain unchanged; no solver or reconstruction reran. Root verified terminal
0, frozen identity match, guard/driver 0, empty Job, every owned host process
absent, cleanup 0 and separate idle/null/no-lease state. The video was opened for
the user without agent inspection. No further visualization job is allocated.

## Historical superseded video spectrum inset, 2026-10-28

The following checkpoint was superseded by moving the spectrum into result.png
and removing carrier arrays from Result. Its old hashes/timings are historical.

Current renderer `561dbd60`, tests `0a961993` and render documentation `4d8ee9f2`
are integrated. Root's 177 integration tests passed in 10.83 s after correcting only the
test fixture's yaw-angle shape; the original failed attempt remains preserved.
No preview or agent visual review was performed. The approved reference image
and script supplied the presentation.

The lower-left inset uses original-carrier amplitudes, angular frequencies and
propagation directions copied into result.npz. Rendering reads that archive
alone; legacy missing-carrier archives fail clearly, while three empty vectors
mean calm water. Coming-from angles add pi, fixed world +x/right and +y/up. Magma
bins show `sqrt(sum(a^2))` at 15 degrees by 0.025 Hz; frequency radius extends beyond
0.5 Hz when needed to include every carrier. Radial labels remain, with no
direction labels, colorbar or time chart. The white-backed monospaced WVHT label
uses `4*sqrt(0.5*sum(a^2))`: this case 1.657467948 m, displayed 1.66 m. Histogram and polar
artists are prepared once; cached inset pixels restore above each moving scene.
Tests check fixed encoded inset pixels under hull translation/yaw and exactly
one full Agg draw/histogram. Existing scene layers and physical geometry remain.

Root's [delivery evidence](../../../artifacts/workflow/sea-inset-delivery1/results/report.json)
and [video log](../../../examples/simple-boat/video.log) confirm 1626 decoded
800x800 frames at 12.5 fps, 130.08 s encoded, zero application warnings. Ordinary
command 332.752 s, CLI 329.157 s, renderer 300.625 s; setup 1.187 s, geometry 206.653 s,
draw 81.823 s. These are whole-render timings, not isolated inset overhead.
Carrier attachment preserved every original NPY member byte; no body solve,
field reconstruction, MATLAB call or calm-analysis rerun occurred for the inset.
The retained field already uses the newly simplified reconstruction.

The delivery guard completed terminal 0/frozen-match/driver 0/Job-empty. Root
verified all owned host processes absent, completed lease cleanup and separately
read idle/null/no-lease state. No further visualization job is allocated. Root
owns final acceptance and publication; historical sections below retain their
original policies and evidence.

## Fixed mean-waterline reference, 2026-10-26

The tested source is integrated into baseline: runtime SHA141046c5, regression
SHA130fd07a, CLI documentation SHA595e14d5. All50 renderer tests passed in
`mean-waterline-preview-attempt3`; its separate four-frame preview was supplied
for the user's review, without agent visual inspection. Attempts1/2 remain as
failed launcher/fixture evidence; the final correction only changed the fixture's
single-axis rotation input shape. No physical solution changed.

The initial posed hull intersects a flat world-z0 surface once to construct a
fixed gray dashed mean-waterline reference. It sits above water/mask and below
the full gray boat mesh, black opening outline and solid blue moving contact.
There is no label or legend, per the user's decision. This geometry is derived
from saved Result inputs without a solver-specific field or new NPZ array.
Separation is a qualitative warning about mean-hull geometry, not an error bound.
The full video job `mean-waterline-full-render-attempt1` passed:751 fully decoded
H264 frames,800x800,12.5fps,60.08s; renderer125.468s, command136.093s. Public
video SHA9e654250c924749188aa184847349d5e98204369519dc41e464d4ef49b1191cf.
Log ends success with zero application warnings; the existing encoder duplicate-r
notice does not affect the verified frame rate. All protected physics inputs are
byte-identical. Exact process/frozen-input/idle checks passed and lease released.

## Completed full wireframe layering correction, 2026-10-26

Latest user explicitly wants transparent dry **faces**, with dry edges retained.
The isolated patch in `artifacts/worktrees/dry-hull-wireframe` is integrated into
the baseline working tree. Its 49 targeted renderer tests passed. Root produced
a separate four-frame preview and the user approved it. User instruction:
previews are for the user's review, not agent visual review; keep them separate
from actual full rendering. No preview was produced during the full render.

The public full render completed and verified all751 frames,800x800,12.5fps,
60.08s. Renderer elapsed134.704s; command149.281s, unchanged240s budget.
`examples/simple-boat/video.mp4` SHA256:
39523c0c617450c3fd5477b82317e692c9fd603f82abb466804bfe3f4dfc3c32.
`video.log` ends success with zero application warnings. The external encoder's
duplicate frame-rate-option notice is benign; actual frame rate/count/duration
were decoded and verified. Original NPZ/PNG/sea/result logs are unchanged.
Both jobs are terminal, frozen identities/processes checked and leases released.

The drawing order is saved water plus existing interior-water contour mask,
then the full gray physical boat mesh, black opening outline, and solid blue
physical waterline on top. No blue halo or dry face fill. The physical contact
kernel/mask, saved field/motion, camera/HUD, budget and publication guards are
unchanged. All physical edges/rim are reused directly; renderer calls skip the
now unnecessary submerged-only edge-classification work. The existing private
diagnostic edge path and its generic geometry tests remain available.

A focused source regression spies on actual Agg draw calls in two saved frames
with an encoder stub: it checks artist order/zorder, complete dry+wetted edge
arrays, gray/black/solid-blue colors, physical-only blue contact, no halo, and a
wet mask narrower than the dry flared rim. It fails if rendering invokes the
submerged-only edge classifier. Evidence: `dry-hull-check-attempt1` and
`dry-hull-full-render-attempt1` under `artifacts/workflow`. Integrated runtime and
test bytes match the tested isolated producer. Numerical physics qualification
remains separate; no solver was rerun for this display correction.

## Current waterline display and chart, 2026-10-25

The canonical video now shows water up to the moving physical hull intersection,
transparent dry hull and visible submerged mesh portions. The physical waterline
is computed once per frame; its oriented segments form one standard compound
Matplotlib PathPatch with one subpath per actual contour. No extra water-triangle
mask polygons or water/hull clipping pass remain. Edge intervals reuse contact
endpoints, existing water-coverage boundary crossings and the validated P1
TriFinder. Camera, HUD convention, saved motion/surface and physics are unchanged.
An overtopped opening can need a separate **display-only** internal-closure trace;
it never changes physical contact or supplies a physical cap/waterline/load.

Root's `wet-hull-contour-tests-attempt2` passed **48/48 in 1.89 s**. The private
four-frame preview passed complete decode: renderer 1.344 s versus 1.765 s for
the prior tile-mask preview; geometry 0.500 s versus 0.814 s. Full public CLI
`wet-hull-contour-full-render-attempt1` then succeeded: **135.422 s renderer,
149.813 s command**, unchanged 240 s budget, all 751 frames decoded at
800x800/H.264/12.5 fps/60.08 s. Root verified frozen identities, exact processes
and FFmpeg terminal status and released each lease. The separate
`wet-hull-keyframes-attempt1` again decoded all 751 frames; root viewed the actual
51.44 s peak-roll frame (6.61 degrees) and preview stills. Current video SHA256:
`74ed47aa4157313e05896cce52d7dbd87501fef6905059a81f3359039eed16e5`.
Canonical NPZ/PNG and protected sea/result logs remained unchanged; rendering
replaced only video.mp4 and its invocation log. See [versioned evidence](waterline-transparency-2026-10-25.md)
and [public video log](../../../examples/simple-boat/video.log).

Chart checkpoint `c361c6f` is pushed: exactly seven nice, symmetric Y ticks,
zero in the middle and limits at the first/last ticks. The canonical shared
translation range is +/-1.5 m; the separate rotation range is +/-7.5 degrees.
Every saved motion sample/startup and existing frame convention are retained.
Chart tests passed 23/23; the private PNG was viewed before publication.

The initial tile-mask full render exceeded 240 s and published no video; its
failed attempt remains preserved. Actual profiling motivated removal of those
extra mask polygons. Progress now logs first/every 50th/last frame and stage
timings; failures report completed frames, phase and accumulated/active timing.
No budget increase or physics rerun was used. Physical accuracy remains
provisional, with the qualifications in [warning triage](../../cli/warning-triage.md).
Top-down single-sheet/contact-contour limits and missing-water checks remain
explicit; display closure does not establish flooding/wetting physics. Root
owns integration/checkpoint scheduling; no further job is allocated to this owner.

## Historical first canonical delivery and shared-axis chart

The following records the earlier delivery before waterline transparency and
seven-tick chart styling. Its artifact hash/timings and test counts remain
historical evidence; current outputs/status are recorded above.

Latest user authorizes the canonical NPZ/directory workflow. The dedicated
visualization_expert is loaded as `/root/result_renderer`. Current public command
is `uv run topwater render [DIRECTORY]`: read `result.npz`, atomically replace
`video.mp4` after full decode/frame/dimension checks, and write `video.log`.
No boat/sea input or hydrodynamic computation belongs to rendering.

The complete canonical public workflow is delivered at checkpoint `bda901d`.
`canonical-pipeline-attempt4` completed the 60 s WEC-Sim solve in 2393.743 s and
independent public render in 226.932 s: 3001 motion samples, 751 surface frames,
and 751 fully decoded H.264 frames at 800x800, 12.5 fps, encoded duration 60.08 s
(endpoint frame hold). Root viewed decoded 0/30/51.44/60 s frames and accepted
coverage and readable HUD. Evidence: `artifacts/workflow/canonical-pipeline-attempt4/`
`report.json`, `terminal.json`, and `output.log`. Exit 0 and frozen identities
were verified; root verified the exact process tree terminal and released the
lease. The handoff and execution ledger own subsequent scheduling/liveness.

`result_render.py` streams packed surface frames directly, supports changing
connectivity, and interpolates only rigid poses (linear position/SLERP). The
saved hull is transformed with CG exactly once. Accepted 800-square top-down
coolwarm/off-white-hull/contact/HUD styling remains. One Agg figure and reused
artists use background blitting; the water artist is rebuilt only on topology
changes. Folded/overlapping projected sheets are explicitly unsupported by this
presentation, without narrowing the result format or native solver physics.

`solve` now also publishes `result.png` through `result_plot.py` and the paired
NPZ/PNG publisher. The final canonical PNG was generated from the completed
Result without rerunning physics: `artifacts/workflow/result-plot-backfill-attempt2/`
`report.json` and `terminal.json`. It includes all 3001 saved motion samples and
startup, uses the prior white/light-grid/teal 2x3 style at 2160x1170, and agrees
with archived native motion curves within 2.63e-14 m or deg. Translations share
one symmetric metre range; rotations share a separate symmetric degree range,
each with 5% headroom. Fixed initial heading axes and initial-attitude removal
preserve the motion convention; relative XYZ Euler labels are unwrapped and
warn at gimbal lock. Root and visualization owner viewed the final shared-axis
PNG. The original seven data/video/log files remained bit-exact; original
`result.log` was not retroactively changed to claim PNG production. Latest
shared-axis/solve-publication checks: 44 passed in 5.07 s, evidence
`artifacts/workflow/motion-plot-style-tests-attempt1/output.log`. Root verified
terminal/frozen state and released both execution leases.

Artifact and visual delivery pass; physical accuracy remains provisional.
The producer preserves the total-water mesh, with hull-interior portions
implicitly excluded, and computes geometric waterline contact separately.
The top-down display requires one upward single-valued projected sheet with
supported boundaries/holes; folded or overlapping layers fail clearly. Opaque
masking/contact annotations do not establish wetting, flooding or nonlinear
fluid physics. Radiation-phase sensitivity, reference mesh/irregular-frequency
risk, reconstructed-field/load discrepancies and linear applicability remain
qualified in [warning triage](../../cli/warning-triage.md) and issues #31-34.
No further plot, render, test or physics job is allocated to this owner.

## Historical pre-delivery projection correction and clipped-export evidence

The following preserves earlier tests and the superseded clipped-export study;
current production preserves the producer's uncut total-water mesh. These
records are not current execution allocations or the final public-run evidence.

The renderer's old normalized XY determinant cutoff `1e-11` imposed a whole-patch
minimum triangle area, rejecting small well-conditioned exported cut faces.
The corrected predicate forms local XY edges before patch normalization and
requires signed determinant greater than `64*eps*|edge1|*|edge2|`. Computed span,
determinant and threshold must be finite. Boundary/view tolerances remain
separate; all winding, topology, crossing, overlap, hole and coverage checks are
unchanged. Independent researcher `/root/startup_domain_policy` approves the
source and mathematics, including nonfinite guards.

Owned correction files are `result_render.py`, `tests/test_result_render.py`, and
the predicate paragraph of `docs/cli/render.md`. New controls cover a conforming
10 m patch containing a 10 micrometre right triangle, clockwise/collinear/
near-vertical/folded projection, and finite-input arithmetic overflow. Root ran
`tests/test_result_render.py tests/test_visualization.py`: **31 passed in1.79s**,
including those controls and one real two-frame variable-topology encode.
Evidence: `artifacts/workflow/renderer-projection-tests-attempt1/output.log` and
`terminal.json`; producer29496/launcher27040, exit0, frozen identities match.
Root independently verified exact process/FFmpeg terminal status and released
the lease. The earlier 228-pass/one-skip suite predates this predicate correction.

Root's cached geometry sweep (session87492) also terminated PASS: all751 saved
total-water frames passed exporter/contact checks in1711.7607s, frozen identities
matched, exact processes were terminal and the lease was released. Evidence:
`artifacts/workflow/cached-geometry-sweep-attempt1/`. This is full cached geometry
acceptance, not a full render or a new body/field computation. The execution
ledger and handoff own subsequent liveness.

At this historical checkpoint the complete public workflow was still pending;
the subsequent successful uncut workflow and final PNG are recorded above.

## Historical incident-only rendering evidence

The following records the prior completed SceneSequence workflow. Its command,
timings and process/lease statements are historical and do not describe current
public CLI or liveness. Detailed evidence is retained in the linked records.

User-approved shared Matplotlib2D workflow, full60s video and read-only speed
review are COMPLETE. Active3D renderer/helpers/PyVista dependencies removed;
physical3D mesh/intersection/solver code and historical evidence remain.
No further job, benchmark, physics change or full-render variant allocated.

Watch `artifacts/reference/stronger-buoy-20261024/matplotlib-full60-attempt1.mp4`,
SHA15510982f6cef3a246f1129222e0525052215a06cfa2e0073af3c53bf803f23b.
751decoded/H264800square/12.5fps/60.08s with endpoint framehold. Owner viewed
0/30/51.44/60stills; ROOT viewed0+peak. Source0..60/fullnativecached trajectory.
Original hull/time/central65waterXYZ BITEXACT, exactnativeq HUD. Same-field
minimalborder only:75x77/1550newnodes/frame/~104MB, actualviewcoverage passed.

CLI: `uv run topwater render data/scene.npz videos/preview.mp4`. Both paths
required positional, relative or absolute; wholeinput/uniformcadence/fixed800,
no --extra/input/output/start/end/size. Only requestedMP4 persists; temporary
encoder scratch cleaned, diagnostics stdout. Optionalsettings onlyframing/timeout.

Style: fixedtopdown/coolwarm symmetricheight about0; offwhiteopaquehull blocks
underlyingheightcolors, graymesh/blackrim/darkblueactualwaterline withwhitehalo.
14.5ptmonospace/fixedsignedfield TIME/ROLL/PITCH/HEAVE spans almost800pixels.
OneAggfigure/reusedartists, backgroundblit anddirectbufferFFmpegstream.

API [visualization.py](../../../src/topwater/visualization.py),
[drawing](../../../src/topwater/_visualization_2d.py),
[intersection](../../../src/topwater/_visualization_geometry.py),
[CLI](../../../src/topwater/cli.py). SceneSequence NPZ holds time_s,
hull_points[T,N,3]/hull_faces[F,3], water_points[T,M,3]/water_faces[K,3], metadata,
telemetry__roll_deg/pitch_deg/heave_m. Metres/right-handed-z-up, fixedconnectivity,
convexhullrequirement validated withboundedstorage. Producer normalizes frames/
units/cadence and supplies telemetry; renderer never guesses pose or fillszeros.
Trianglecontact uses actualphysicalfaces; coplanararea contact is explicitly
unsupported. No new flooding/wetting/fluid-accuracy qualification.

FullCLI148.792s (renderer145.75), one-timeborderprep40.657s, wholejob197.911s.
Geometry/update64.066s, Agg68.732s, pipewait.530s, encoderflush1.703s,
encoderstartup1.937s, hullvalidation5.359s, setup3.391s. Encoderwork overlaps
drawing; pipewait is not encodingCPU. Initial63tests17PASS, newerHUD/footprint
checks2PASS; finalfull751decode/provenance/coverage PASS. Lastjob session42135
terminal/frozen/idle/separate58a912/FINALreleasebacd1e. Ledger idle/no lease.

[Read-only speed review](render-speed-review-2026-10-25.md): exactcoveragefinder
reuse and tightercontact broadphase first; preset/imagefastpath needquality or
interpolation review. No performance changes or experiments duringfulljob.

Initialimplementation owner `/root/visualization` is topwater_expert; WECdata
adapter owner `/root/wecsim_reference`. Dedicated visualization_expert definition
andKB/routing exist and inherit parentsettings; currentnativeclient rejected
boundedspawn `unknown agent_type 'visualization_expert'`. Configured NOTloaded;
recheck afterclientrefresh, transferexplicitly atstablecheckpoint. No fallback
or separateCLIcoordinator. See [activation](../activation.md),
[library evidence](../topwater_expert/visualization-workflow-2026-10-24.md).
