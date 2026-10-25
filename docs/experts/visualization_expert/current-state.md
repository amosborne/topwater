# Visualization expert: current state

## Canonical result renderer source checkpoint, 2026-10-25

Latest user authorizes implementation of the canonical NPZ/directory workflow.
Dedicated visualization_expert is loaded for `/root/result_renderer`; source
ownership transferred for new `result_render.py`, its tests and CLI render page.
Previous SceneSequence rendering remains available internally and unchanged.

New source streams packed Result surface frames directly, interpolating only
body poses (linear position/SLERP). No boat/sea files, hydrodynamics, frame-wise
full-sequence hull copies, or wave reconstruction. Physical hull stored once in
result transforms with CG exactly once. Accepted800square/topdown/coolwarm/HUD/
offwhitehull/contact styling retained. Persistent Agg figure/blit; water artist
reused for unchanged topology and rebuilt on connectivity changes. Projection
guards include connected manifold sheet/loops, upward nondegenerate triangles,
crossing boundaries, outer-edge intrusion and uncovered holes outside hull.
Some folded native output is intentionally unsupported by this 2D presentation.

Atomic video replacement occurs only after complete decode/frame/dimension
checks. Temp encoding stays outside example; old video survives failures.
`docs/cli/render.md` documents inputs, timing, defaults, fatal detection, warning
and display limits. `tests/test_result_render.py` adds source tests for pose
mapping, variable topology, interior holes, folds, cadence and old-video
preservation. NO tests/imports/renders/probes have run under this assignment;
root must allocate sole numerical lease before validation. Source acceptance
and performance remain pending. No native/sea/physical-result change claimed.

## Previous completed workflow

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
