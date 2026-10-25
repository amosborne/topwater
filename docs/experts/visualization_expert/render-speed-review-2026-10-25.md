# Read-only rendering speed review

User requested this review in parallel with the full60s render. The researcher
inspected source and installed Matplotlib/imageio-ffmpeg implementations only:
no edits, tests, benchmarks, or concurrent numerical jobs. These are candidates,
not demonstrated speedups. Do not modify the frozen full-render producer.

The previous63-frame warm CLI took12.205-14.027s. Geometry/artist updates cost
3.326-4.132s, Agg drawing4.283-4.773s, pipe waiting.048-.110s, encoder finalization
2.187-2.500s. Encoder work overlaps drawing: pipe waiting is not encoding CPU.
The full run's compact75x77 grid replaces the earlier233x233 field, so those
old timings overstate the remaining potential savings.

The completed compact full751-frame run took148.792s at the public CLI:
geometry/update64.066s, Agg68.732s, pipewait.530s, encoderflush1.703s,
encoderstartup1.937s, hullvalidation5.359s, setup3.391s; renderer145.75s.
Input preparation40.657s and totalproducer+QA197.911s were measured separately.
These full-run figures remain the better baseline for later experiments; do not
infer a per-frame speedup merely from fewer input cells.

## Ranked next experiments

1. **Reuse the coverage search.** `_visualization_2d.py` currently builds a new
   Matplotlib Triangulation/TrapezoidMapTriFinder each frame for just four view
   corners. Cache by complete XY coordinates, connectivity and mask. A verified
   rigid translation permits querying one reference finder in translated
   coordinates; deformation needs the existing fresh-search fallback. A bounding
   rectangle alone cannot prove coverage when meshes contain holes. High
   confidence avoidable setup work; benefit remains unmeasured.
2. **Tighten the waterline broadphase.** `_visualization_geometry.py` scans
   already-visible water triangles for each hull face. First select a conservative
   full-hull3D bounding box; structured/stationary XY may permit reusable spatial
   bins. Preserve exact triangle intersection, tolerances and coplanar handling.
   This is contact culling in addition to existing view culling. Data-sensitive
   geometry savings; do not assume the whole update cost disappears.
3. **Test an encoder speed preset.** Try standard x264 `veryfast` in an isolated
   later comparison, keeping800x800, cadence, bitrate and thread count fixed.
   At fixed bitrate, quality may decrease. At fixed quality, files may grow.
   Compare full CLI/send/finalize time and decoded images, not close time alone:
   less lookahead can move work earlier. Imageio0.6.0 close polls every10ms;
   there is no artificial two-second sleep to delete.
4. **Consider an image artist for regular grids.** TriMesh Gouraud drawing gathers
   triangle coordinates and colors on every draw. A regular-grid AxesImage fast
   path could reduce Agg cost while preserving the generic triangular fallback.
   Bilinear image interpolation differs from piecewise-triangle interpolation;
   keep actual waterline geometry unchanged and obtain visual acceptance. This
   has more appearance risk than the exact first two changes.

After the full job is terminal, compare one change at a time on the same saved
scene and phase timers. For exact geometry changes compare contact segments and
peak decoded frames; test interior holes and deforming XY for search reuse.
No new physical-accuracy gate or solver work follows from this graphics review.

Primary references: [Matplotlib collections](https://matplotlib.org/stable/api/collections_api.html),
[Agg blitting](https://matplotlib.org/stable/users/explain/animations/blitting.html),
[imageio-ffmpeg close implementation](https://raw.githubusercontent.com/imageio/imageio-ffmpeg/v0.6.0/imageio_ffmpeg/_io.py).
