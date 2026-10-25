# Standard visualization workflow

## Current standard: Matplotlib top-down

### Completed full 60-second reference

Watch `artifacts/reference/stronger-buoy-20261024/matplotlib-full60-attempt1.mp4`,
SHA256 `15510982f6cef3a246f1129222e0525052215a06cfa2e0073af3c53bf803f23b`.
All751 frames decode as H264800x800/12.5fps/60.08s, retaining the physical60s
endpoint. Owner viewed decoded0/30/51.44/60s; ROOT reviewed0+peak.

Full public CLI wall time was148.792s; renderer145.750s. Same-field border/input
preparation took40.657s and producer+CLI+QA197.911s. Recorded renderer phases:
input load/schema validation.157s, convex validation5.359s, imports/setup3.391s,
geometry/artist updates64.066s, Agg draw68.732s, pipe send.530s,
encoder startup1.937s, close/final flush1.703s. Pipe sends are transfer and
backpressure, not isolated encoder CPU. Encoding overlaps drawing.

The full input is104587856bytes (~104.6MB decimal), with751frames on75x77 water
nodes. Original65x65 waterXYZ, hull poses and time are bit identical to the
accepted reference cache. Only1550missing border nodes per frame were evaluated
from the same immutable train, rather than generating a96m field for every time.
The final fixed camera halfspan is6.100702m; water halfwidths6.9375/7.125m
cover its corners at every frame. Native HUD rows are exactly0:4:3000.
Rendering culls to10366visible of11248waterfaces.

Prepared `scene-full60-attempt1.npz` SHA256
`884667cb0a43f55ec86303c38fae303a76c58a8869c547b1d3290f70b4f330bc`.
Reproduce without input preparation:

```powershell
uv run topwater render artifacts/reference/stronger-buoy-20261024/scene-full60-attempt1.npz videos/reference-full60.mp4
```

Evidence is under `matplotlib-full60-attempt1-source/proposal.json` SHA256
`f03211a6f1fd83d7ce34b76f753cfc5f26a159e19c31ef0a9284e05b1072febe`
and `matplotlib-full60-qa-attempt1/verification.json`. Exact process/children
terminal6f2f47/CIM6c6e36, frozen inputs/dependencies verifiedb012d4, separate idle
readback58a912, final lease releasebacd1e. No numerical work remains allocated.
The read-only speed review did not edit or benchmark this frozen run.

Final integration check after removing 3D API/dependencies ran exactly
`tests/test_visualization.py` and `tests/test_package.py`: **19 passed in0.15s**
on the current cleaned source, with no fixes or rendering repeats needed.
Manifest `render-final-tests-attempt1-source/proposal.json` SHA256
`1e32810f9952b5334b97ded502500e80368424254a79184fe466a57b2d8a429f`;
terminal101701/CIM73aa68, frozen check531e70, separate idlec617a1,
final lease release6a9c91. Its completion receipt is beside the proposal.

### Accepted short preview and styling history

Latest accepted styling is opaque off-white projected hull fill, gray mesh,
black rim, dark-blue geometric waterline, muted standard `coolwarm` elevation
palette about zero, and14.5pt DejaVu Sans Mono header with fixed-width signed
fields. The font is constant throughout the file; values do not shift columns.
The fill comes from the input hull's convex XY footprint, never a boat fixture.

Latest movie under the stronger reference root: `matplotlib-clip-attempt3.mp4`,
SHA256 `7134d6443cc8d66e3e53d2719f63618bf087cb432127f1edb54eb20711dec7b0`.
ROOT and owner reviewed the decoded peak image; full63-frame decode passed.
The earlier17 contract/geometry/CLI tests and2new footprint/header tests passed;
the final font-only rerender did not repeat unrelated checks.

### Complete-process timing

Actual clean public CLI wall time, measured externally from launch to exit, was
**14.027s**. Renderer wall time was12.250s; CLI `main` time was12.375s, excluding
interpreter/package imports before main. The wrapper including QA took14.770s.

| Recorded stage | Seconds |
| --- | ---: |
| Input load and schema validation |0.110 |
| Convex mesh validation |0.079 |
| Renderer imports/setup/initial canvas |0.671 |
| Geometry and artist updates |4.132 |
| Agg drawing/blitting |4.773 |
| Encoder pipe sends |0.048 |
| Encoder startup |0.047 |
| Encoder close/final flush |2.500 |

Encoding runs concurrently with plotting. Pipe-send time measures transfer and
backpressure, not independent encoder CPU usage; these stages do not constitute
a serial encoder profile. Process overhead outside CLI main is also included
in the external total. No whole60-second rendering cost is established here.

Practical next candidates, not implemented or benchmarked: cache contact/coverage
geometry for repeated style iterations (the4.132s includes other updates too),
use a faster standard encoder preset to reduce the flush tail with a quality
tradeoff, or introduce a coarser display mesh/structured-grid drawing path to
reduce Agg cost while retaining fixed800 output. The0.048s pipe wait does not
justify another Python encoder thread. No optimization campaign is allocated.

Exact final manifest `matplotlib-clip-attempt3-source/proposal.json` SHA
`9a77d47fe82c07b54763bd51b6958854c7187eb2290a5ca80075b1349f1b4e87`;
`matplotlib-qa-attempt3/verification.json` contains external and stage timings.
Same prepared telemetry scene and all geometry/time values remain unchanged.

### Original top-down baseline and input contract

```powershell
uv run topwater render data/scene.npz videos/preview.mp4
```

Both paths are required, relative or absolute. The whole prepared file renders
at its uniform saved cadence to a single800x800 MP4. Matplotlib and the bundled
FFmpeg encoder are core dependencies; no extra flag, global FFmpeg installation,
GPU or PyVista is needed for the standard view. Diagnostic information returns
on stdout. The renderer creates no permanent sidecars.

SceneSequence uses world metres in right-handed z-up coordinates:
`time_s[T]`, `hull_points[T,N,3]`, `hull_faces[F,3]`, `water_points[T,M,3]`,
`water_faces[K,3]`, JSON metadata and an optional `telemetry` dictionary of finite arrays with one value per timestamp.
The standard header requires `roll_deg`, `pitch_deg`, `heave_m`; time comes from
`time_s`. NPZ serialization stores these as `telemetry__roll_deg`,
`telemetry__pitch_deg`, `telemetry__heave_m`. Older geometry-only scenes load,
but must receive producer-defined telemetry before the default header can render.
No mesh-based pose inference or first-frame rezeroing occurs in the renderer.
Connectivity is fixed within a sequence; the hull shape and resolution are data,
not renderer constants. All physical hull faces must lie on a convex boundary;
nonconvex/internal faces fail validation. Open gunwales are permitted. Unit/frame
conversion and any time resampling belong to the producer adapter.

The original top-down baseline used Gouraud-shaded standard Matplotlib TriMesh with continuous `bwr`
colors; the latest style uses muted `coolwarm` with the same symmetric
normalization about zero for the full input. Values
are surface elevation in metres, not bathymetric depth in this infinite-depth
reference. Hull edges project in light gray, open boundary edges form the black
rim, and blue contact lines intersect actual physical hull and water triangles.
No convex closing-cap waterline is invented. Hull wireframes are annotations
projected over the water field, not a simulated interior flood surface.

Contact geometry cuts each candidate triangle pair by the other's normalized
plane, then intersects their collinear one-dimensional intervals. Bounding-box
culling bounds candidate storage per physical face. Point contacts have zero
line length; coplanar area contact has no unique waterline and fails explicitly.
This also handles a hull crossing a large water triangle without the earlier
point-distance mask limitation.

The camera has fixed world XY orientation and fixed scale/center over the input;
translation remains visible. All-time triangle bounds cull water outside that
view before plotting. A single figure, TriMesh, LineCollections and header are
reused. A cached Agg background is restored, animated artists are drawn, and the
RGBA buffer streams directly to the public `imageio_ffmpeg.write_frames` encoder.
There is no second full-figure redraw after blitting and no per-frame image files.

Accepted movie `matplotlib-clip-attempt1.mp4` under the stronger reference root,
SHA256 `954343c52b6bc5690fcf421bc0ac6804c2a4643b2571c26ffc34682f44e87eed`:
63 frames,800square,12.5fps,5.04seconds from50.00–54.96s. All17 focused schema/
CLI/geometry/package tests pass; full decode and ROOT/owner peak-still review pass.
Geometry/time/connectivity are bit identical to the accepted scene; producer
telemetry matches native WEC-Sim q rows exactly.

Actual CPU render11.25s includes validation.094s/setup.765s; frameaverage.120s.
The frame loop totals 7.563 s: geometry/artist updates 3.499 s (0.0555 s/frame)
and combined Agg drawing/encoder streaming 4.064 s (0.0645 s/frame). Encoder
startup/flush and initial canvas work are included in the total, not separately
profiled. These are recorded timings, not a new profiling experiment.
Only8064/107648watertriangles are drawn. Combined tests, telemetry augmentation,
CLI and QA took17.797s. Wider water preparation remains the earlier one-time
100.75s cache operation and was not repeated. Evidence:
`matplotlib-clip-attempt1-source/proposal.json` SHA
`cc6ba729f6104f635f453cd58dd70a44fc9fc26ae9bc8f3fd1fa0ae890a2b589`,
`matplotlib-qa-attempt1/verification.json` and `augmentation.json`.


### Retired 3D renderer

The user approved removal of active 3D rendering after reviewing the 2D result.
The PyVista/VTK renderer, its settings and optional dependencies are removed.
Git checkpoints8450c2c/662a148 preserve their historical source and evidence;
previous artifacts remain unchanged. Three-dimensional mesh inputs, physical
geometry and actual triangle-intersection math remain part of the 2D workflow.
No physical solver or hull definition was removed.
