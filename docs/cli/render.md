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

## Output and replacement

`video.mp4` is H.264, 800 by 800 pixels, with one frame per saved surface sample.
Frame rate is the inverse surface time spacing. The last saved state is shown
for one frame interval, so encoded duration is `frame_count * spacing`, one
interval longer than the source endpoint. There is no hidden prerun or cropping.

Encoding occurs in the system temporary area. The finished video is fully
decoded and checked for dimensions, cadence, and frame count before atomic
replacement. A failed or interrupted render preserves an existing video.
Cross-volume publication uses a temporary staging directory outside the example
directory. No scene cache, images, reports, or sidecar files are retained.

The command prints its output path. The Python rendering function additionally
returns preparation, setup, geometry, drawing, pipe, encoder and verification
timings; these do not become fields in the result archive.

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
- Opaque off-white hull footprint, gray physical mesh, black opening rim, and
  dark-blue geometric hull/water intersection with a white halo.
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
| At each saved frame | Folded/degenerate XY projection, inconsistent face winding, non-manifold/disconnected sheet, crossing boundaries, overlapping/invalid triangulation, or an outer boundary entering the view. |
| At each saved frame | Uncovered interior water hole outside the opaque hull footprint, or water missing a fixed view corner. A boat-shaped hole is supported. |
| At contact calculation | Degenerate contact triangles or coplanar hull/water area contact without a unique waterline. |
| During encoding/verification | Time budget exceeded, missing encoder/dependency, filesystem/pipe failure, unreadable video, wrong dimensions/cadence or incomplete frame count. |

All failures occur before replacing the named video. Projection checks use
span-normalized geometric tolerance `1e-11`; contact annotations use the shared
intersection helper's tolerance `1e-9` times hull span. These are numerical
geometry tolerances, not physical-accuracy acceptance thresholds.

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
