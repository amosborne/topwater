# Waterline transparency: measured delivery, 2026-10-25

The canonical top-down video now displays saved water up to the moving hull
intersection. Dry hull, including dry rim/mesh lines, is transparent; submerged
physical edge portions remain visible. Camera, HUD convention, all saved motion
and water geometry, and solver physics are unchanged. This report uses root's
completed serial execution evidence, with UTC execution on 2026-10-26.
Final renderer source SHA256:
`a31179fa92c4593cc65d75157f0350797d4f96b165e026a933e3a401c67d8abb`.

## Final display algorithm

The existing shared physical waterline kernel runs once per normal frame.
SciPy endpoint matching joins its oriented segments into closed contours. One
standard Matplotlib compound PathPatch preserves concavity, disconnected
components and nested winding. No per-water-triangle mask polygons, convex
replacement of the contact, or second water/hull clipping pass remains.

Mesh edges split at those existing 3D contact endpoints and existing 2D water
coverage boundaries. Batched interval midpoints reuse the validated TriFinder
and represented P1 height, so edge visibility includes multiple crossings and
cannot bridge missing water. Uncovered legacy hole portions are omitted rather
than assigned invented elevations; exposed missing water fails clearly.

The posed/sloping regression partly overtops its opening, leaving genuine
physical contact open. Only this specifically open-contact path uses the existing
internal closure and intersection core to complete the **display-volume** trace.
Returned physical contact and physical edge classification stay unchanged.
Normal closed-contact frames perform no closure work or extra intersection.
Caps never become physical waterlines, pressure surfaces or flooding physics.
Branched/crossed/touching/inconsistent winding and unresolved open boundaries
fail clearly. Cap reviewed contour mathematics/topology/winding and identified
the legacy coverage gap; the coverage and display-closure corrections are covered
by root's subsequent 48-test passing suite.

## Evidence and timings

| Evidence under artifacts/workflow/ | Actual result |
|---|---|
| wet-hull-contour-tests-attempt2/output.log | 48 passed in 1.89 s, including legacy hole coverage, exact sloping/posed contact, overtopped display-only closure, one-call normal path, contour winding/concavity/components, atomic output and progress/failure reporting. |
| wet-hull-contour-preview-attempt1/preview/report.json | Four unchanged saved surfaces at 0/20/40/60 s, complete four-frame decode, 1.344 s renderer and 0.500 s geometry. Prior tile-mask preview was 1.765 s renderer/0.814 s geometry. |
| wet-hull-contour-full-render-attempt1/{output.log,terminal.json} | Public CLI success, 135.422 s renderer/149.813 s command, unchanged default 240 s rendering budget; exit 0, frozen identities matched. |
| wet-hull-keyframes-attempt1/report.json | Separate complete 751-frame decode and stills at 0/30/51.44/60 s. Root viewed actual 51.44 s peak roll, 6.61 degrees, and the private preview stills. |
| seven-tick-plot-tests-attempt1/output.log | 23 passed in 13.56 s. Seven-tick chart checkpoint c361c6f pushed; private PNG viewed before publication. |

Final video: H.264/YUV420p, 800x800, 751 frames, 12.5 fps, 60.08 s with endpoint
frame hold. Backend is **Matplotlib Agg CPU blit**, Matplotlib 3.11.2; decode
reports FFmpeg 7.1 essentials Windows build. Video SHA256:
`74ed47aa4157313e05896cce52d7dbd87501fef6905059a81f3359039eed16e5`.
The [public video log](../../../examples/simple-boat/video.log) records progress,
publication, verified frame count and actual timing. Its final frame counters
include validation 56.154 s, mask/contours 21.105 s, path construction 0.140 s,
edges 2.002 s, contact 13.170 s, Agg drawing 37.005 s and encoder send 0.814 s.
Encoder finalization was 1.750 s and full verification 2.328 s. Pipe wait is not
encoder CPU time; these phase totals are measured for this saved case only.

No physics was recomputed. Protected boat/sea/result NPZ, result PNG, sea.log
and result.log stayed unchanged. Video and video.log were deliberately replaced
by the successful public render. Root verified exact producer/child/FFmpeg
terminal status, frozen identities and lease release for completed jobs.

The chart retains all 3001 saved motion samples/startup, with shared +/-1.5 m
translation and +/-7.5 degree rotation ranges. Both rows have exactly seven
nice ticks, zero centered, outer ticks equal to the axis limits.

## Preserved failures and limits

`wet-hull-full-render-attempt1` exceeded the original 240 s budget and published
no new video. Its evidence is preserved. Four-frame cProfile in
`wet-hull-profile-attempt2/profile/report.json` measured 0.3526 s display work,
including 3291 polygon clipping calls/0.1504 s, contact 0.0831 s and edge clipping
0.0404 s, with 772-857 mask polygons per frame. Those extra tiles were removed;
the budget was not raised. `wet-hull-contour-tests-attempt1` also remains intact:
45 passed/one failure exposed the real opening-overtopping display case, which
was fixed without changing that fixture or loosening tolerances.

Full decode and visible geometry are artifact acceptance, not hydrodynamic
accuracy. The renderer still requires a supported upward single-valued water
sheet and unambiguous projected contours; the Result format remains broader.
Radiation phase, mesh/irregular-frequency risk, reconstructed-field/load
consistency and linear applicability remain qualified in
[warning triage](../../cli/warning-triage.md) and their issues. No further test,
render, physics or performance campaign is allocated by this record.
