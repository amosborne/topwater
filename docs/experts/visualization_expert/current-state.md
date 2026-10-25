# Visualization expert: current state

## Current canonical renderer, 2026-10-25

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
