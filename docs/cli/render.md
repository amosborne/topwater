# `topwater render`

Render the saved hull motion and total free-surface geometry in `result.npz` to
`video.mp4`. Rendering reads no boat or sea input and performs no hydrodynamic
calculation. Editing `boat.npz` or `sea.npz` does not reinterpret an old result.

## Usage

```powershell
uv run topwater render examples/simple-boat
```

The directory is optional and defaults to the current working directory:

```powershell
cd examples/simple-boat
uv run topwater render
```

The filenames are fixed. Relative and absolute directory paths are supported.
There are no public camera, size, interval, FPS, or sea-reconstruction flags.

## Input

`result.npz` is the numeric, self-contained result archive. Its hull is the full
physical body surface used for this result, which may have been remeshed. It
contains no artificial opening cap. The archive supplies:

- Body times `time_s`, absolute world CG positions `position_m`, and unit active
  body-to-world quaternions `orientation_wxyz` in `(w,x,y,z)` order.
- Body-frame `hull_vertices_m`, `hull_faces`, and body-frame `cg_m`.
- Saved `surface_time_s` and packed world-space `surface_vertices_m` and
  `surface_faces`, with per-frame vertex/face offsets. Face indices are local
  to each frame. Counts and connectivity may change between frames.

Positions and vertices use metres, time uses seconds, and world axes are
right-handed east/north/up. Hull vertices transform as
`world = position + rotation @ (vertex - cg)`.

Body and surface times start at zero and have the same ending time. Surface
times must be uniformly spaced for a constant-cadence movie. Body poses may
have a different cadence: translation is interpolated linearly and orientation
uses shortest-path quaternion SLERP at the saved surface times. Water geometry
is never interpolated or recomputed. NPZ loading disables pickle.

The water mesh need not end at the hull. Portions inside the posed hull are
invalid fluid. A display-only mask hides the surface where its actual elevation
intersects the convex hull volume, so visible water reaches the moving waterline
instead of stopping at the dry gunwale silhouette. The shared physical-hull/
surface intersection routine still checks contact. Rendering never changes or
remeshes the saved water triangles or hull.

## Output and replacement

`video.log` records the latest invocation, rendering/verification stages,
timings, errors and publication status using [the common format](logging.md).
It is replaced when rendering starts, including failed attempts.
Frame progress is printed after the first frame, every 50 frames and the last
frame, including cumulative validation, wet-mask geometry/path, submerged-edge,
shared-contact, drawing and encoder-send times. A rendering failure prints its
phase, completed frame count, elapsed/active-phase time and accumulated stage
timings, including encoder finalization and full video verification. Failed
phase time is also visible when that phase could not complete its accumulator.

`video.mp4` is H.264, 800 by 800 pixels, with one frame per saved surface sample.
Frame rate is the inverse surface time spacing. The last saved state is shown
for one frame interval, so encoded duration is `frame_count * spacing`, one
interval longer than the source endpoint. There is no hidden prerun or cropping.

Encoding occurs in the system temporary area. The finished video is fully
decoded and checked for dimensions, cadence, and frame count before atomic
replacement. A failed or interrupted render preserves an existing video.
Cross-volume publication uses a temporary staging directory outside the example
directory. The named video and execution log are the retained outputs.

The command prints its output path. The Python rendering function additionally
returns preparation, setup, geometry, drawing, pipe, encoder and verification
timings; these do not become fields in the result archive.
The log gives a compact video summary (including distinct source and encoded
durations) and stage timings rounded to milliseconds.

## Assumptions and hard-coded presentation

- Matplotlib Agg CPU rendering and bundled FFmpeg via `imageio-ffmpeg`; no GPU
  rendering is claimed.
- Fixed top-down square view, equal metres, centered over the entire sampled
  hull-motion envelope, with framing factor 1.45.
- Whole saved sequence, default total rendering budget 240 seconds, H.264
  3 Mbit/s, YUV420p, one encoder thread. Encoder finalization is bounded to
  30 seconds (or the smaller configured Python rendering budget).
- Continuous coolwarm water elevation colors, fixed symmetrically about zero
  using the entire saved surface range. No wave or motion amplification.
- Off-white mask only where the saved water lies inside the posed hull volume.
  Dry hull faces are transparent; their mesh and rim lines remain visible.
  Bottom-to-top drawing order is water surface (with its interior-water mask),
  a dashed gray reference waterline, the full gray physical boat mesh, the black
  opening outline, then the solid
  blue physical waterline. All dry and wetted edges are shown. The blue line
  has no halo and uses physical contact only, never an artificial closure trace.
  Water still reaches the actual moving hull intersection. This is a top-down
  schematic, without optical refraction or depth shading.
- The dashed gray reference is the initial posed physical hull intersected with
  flat world z=0, computed once and fixed in world coordinates. It is independent
  of the first saved, possibly wavy surface. For the WEC-Sim workflow it shows
  the mean-hull waterline used by the linearized model; for other producers it
  is simply an initial-pose reference, without assuming their solver is linear.
  Separation from the moving blue line is a qualitative geometry comparison,
  not a force or motion accuracy bound. No label, legend, metric, or HUD change
  is added. If the initial hull does not intersect z=0, the reference is empty.
- TIME/ROLL/PITCH/HEAVE header. Roll/pitch are XYZ Euler departures from the
  initial attitude expressed in its initial horizontal heading frame. Heave is
  change in world CG height. Initial equilibrium trim is not wave-driven pitch.
- One persistent figure with reused hull/header artists and Agg background
  blitting. Water artists are reused when connectivity is unchanged and rebuilt
  when topology changes. A full figure redraw is not performed for each frame.

## Fatal errors and detection

| When | Detection |
|---|---|
| Load | Missing/malformed NPZ members, invalid arrays, nonfinite values, invalid hull or surface indices, invalid quaternions or incompatible time endpoints. |
| Before drawing | Nonuniform surface cadence, undefined initial horizontal heading, invalid output path or unsupported framing/budget. |
| Before drawing | Ambiguous coplanar physical-hull/flat-z=0 intersection prevents a unique reference waterline. |
| At each saved frame | Folded/degenerate XY projection, inconsistent face winding, non-manifold/disconnected sheet, crossing boundaries, overlapping/invalid triangulation, or an outer boundary entering the view. |
| At each saved frame | Uncovered interior water hole outside the opaque hull footprint, or water missing a fixed view corner. A boat-shaped hole is supported. |
| During display masking | A saved water hole extends outside the hull volume at its boundary elevation and would be exposed by transparent dry hull. Missing water is not invented. |
| During display masking | The display boundary remains open/branched after any opening-closure completion; contours cross/touch or have inconsistent nested winding. The renderer fails instead of convexifying the waterline. |
| At contact calculation | Degenerate contact triangles or coplanar hull/water area contact without a unique waterline. |
| During encoding/verification | Time budget exceeded, missing encoder/dependency, filesystem/pipe failure, unreadable video, wrong dimensions/cadence or incomplete frame count. |

All failures occur before replacing the named video. Triangle projection checks
require positive signed XY doubled area exceeding
`64 * machine_epsilon * projected_edge_length_1 * projected_edge_length_2`.
Edges are formed before whole-patch normalization, so a small, well-conditioned
cut triangle is not rejected merely for being small relative to the patch.
Nonfinite computed spans, determinants or area thresholds are rejected explicitly.
Boundary-distance/crossing checks retain the span-normalized `1e-11` tolerance;
contact annotations use the shared intersection helper's `1e-9` times hull span.
These are numerical geometry tolerances, not physical-accuracy thresholds.

## Warnings

If a saved orientation reaches Euler gimbal lock, the command warns that HUD
angles are nonunique. Geometry still uses the saved quaternions. The renderer
does not invent linear-range, flooding, radiation, or wetted-area warnings:
those physical diagnostics belong to `solve`.

## Detection limits

The NPZ format can represent more general surfaces than this top-down renderer.
This presentation supports one connected, upward, single-valued projected
water sheet with simple boundary loops and hull-hidden holes. Folded or
overlapping layers are rejected, rather than flattened into a misleading view.
That is a visualization limit, not a restriction on a future native solver.

Waterline segments and opaque hull masking are presentation geometry. They do
not establish correct flooding, wetting, wave scattering, nonlinear contact or
physical accuracy. Full video decode verifies the artifact, not the physics.

The shared physical waterline is computed once. Its oriented endpoints are
matched within `1e-9 * hull_diagonal` using SciPy, then assembled into closed
contours. One standard Matplotlib compound PathPatch draws those contours,
preserving concavity, disconnected components and nested holes. No additional
water-triangle/hull clipping or per-water-triangle mask polygons are generated.
The full physical mesh and opening outline are drawn without wet/dry edge
clipping; restoring dry wireframe visibility adds no edge/water intersection
work. Hull closure planes reject exposed legacy water holes. When actual physical contact ends at an
overtopped opening, the existing internal closure supplies a separate trace to
close the display-volume mask. This fallback intersects only closure faces; it
never changes returned physical contact, draws a physical cap, or claims flooding
physics. Closed physical contact takes no extra closure/intersection call.
The rendering budget remains 240 seconds.
For a saved interior-water hole fully hidden inside the hull, the hole remains
masked. Drawing the full physical wireframe does not invent water elevations in
that hole. The canonical full-grid result does not have these holes. Transparent
hull faces require actual saved water in any newly exposed region.
