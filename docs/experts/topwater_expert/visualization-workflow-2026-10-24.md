# Standard visualization workflow

## Current standard: Matplotlib top-down

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

SceneSequence retains the neutral mesh contract documented below and adds an
optional `telemetry` dictionary of finite arrays with one value per timestamp.
The standard header requires `roll_deg`, `pitch_deg`, `heave_m`; time comes from
`time_s`. NPZ serialization stores these as `telemetry__roll_deg`,
`telemetry__pitch_deg`, `telemetry__heave_m`. Older geometry-only scenes load,
but must receive producer-defined telemetry before the default header can render.
No mesh-based pose inference or first-frame rezeroing occurs in the renderer.

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

The earlier optional3D implementation and artifacts below remain preserved
until the user reviews2D and explicitly elects removal. They are not the default.

## Earlier optional 3D workflow and preserved evidence

User-authorized library visualization only. Solver implementation and physical
model qualification remain held. Initial standard workflow is executed and
accepted; source and exact run evidence below distinguish tuning attempts.

`topwater.visualization` accepts prepared `SceneSequence` world meshes in metres,
right-handed z-up coordinates: increasing `time_s`, `hull_points[T,N,3]`,
`hull_faces[F,3]`, `water_points[T,M,3]`, `water_faces[K,3]`, and JSON metadata.
Each sequence has fixed triangular connectivity; hull geometry and resolution
are data. Producer adapters convert their poses/surfaces into this contract.
The library has no WEC-Sim references, wave model or current-boat dimensions.
`save`/`load` use pickle-free NPZ so camera/material changes reuse prepared data.

Hull convexity is an explicit requirement. All supplied physical faces must lie
on the convex boundary. Open gunwales are supported; the derived closed envelope
masks water for display only and does not model flooding or physical wetting.
Unsupported internal/nonconvex faces fail explicitly. Unit/frame conversion and
time resampling belong to the adapter, never silent renderer adjustments.
Water mesh resolution must resolve the hull contact region: standard VTK
point-distance clipping interpolates within cells and cannot discover a small
hull wholly inside a coarse triangle. Prepared WEC data retains its fine local
water grid. This is a display-resolution requirement, not a fluid-accuracy claim.

`RenderSettings` controls camera/framing, materials and light/shadows. Resolution
is fixed at 800 x 800 and cadence is inferred from uniform saved timestamps.
`render_sequence` lazily imports optional PyVista,
reuses one VTK context and actors, updates mesh data, uses standard VTK clipping,
normals, boundary extraction, lighting/shadows and transparency. Standard VTK
cell picking checks actual water coverage at all four camera corners per frame.
Only the explicitly named MP4 persists. Backend capability text and separate
setup/frame timings are returned to the caller/stdout; temporary encoding files
are cleaned automatically. Diagnostic capture belongs to an explicit tuning job.
No bespoke rasterizer, shader, fake background water or vertical amplification.

Install graphics with `pip install 'topwater[visualization]'`. Ordinary package
imports do not require or initialize graphics. Actual backend and timing were
measured in the allocated job below; installed GPU names alone are not proof
of the backend used.

Code-first CLI after a producer prepares its scene:

```powershell
uv run --extra visualization topwater render data/scene.npz videos/preview.mp4
```

Both input and output paths are required positionals. `--settings` is optional
and takes a JSON `RenderSettings` object for camera/material tuning. Both relative
and absolute paths are supported. There are
no input/output flags, time-range flags, size or retiming options. The entire
prepared input file is rendered. `topwater render --help` opens no graphics
context. At least two uniformly spaced timestamps are required; irregular
cadence is rejected, with resampling a producer-adapter responsibility.

## Executed preview and limits

- 13 focused schema/package/CLI tests passed, including required positional paths,
exact output filename, full input retention, convex boundary rejection and
irregular-cadence rejection. A distinct convex box used the same public CLI.
- Initial attempt1 stopped before rendering because the existing global pytest
temporary directory was inaccessible. Runner-local `--basetemp` fixed it.
- Attempt2 produced valid meshes/video but its shadows, boundary lines and
three-bin listed colormap were visually rejected. Artifacts remain preserved.
- Accepted attempt3 uses standard continuous `Blues_r`, readable ordinary VTK
material lighting, default shadows/boundary lines off and perspective. No new
custom graphics code or physical changes were introduced.

Accepted movie:
`artifacts/reference/stronger-buoy-20261024/pyvista-clip-attempt3.mp4`, SHA256
`c80e14001a9e3818291e428cea71378b3e158b44831d084afcc34b67d1ec2b0f`.
All 63 frames decode at 800 square / 12.5 fps / 5.04 seconds; source time 50.00–54.96 s.
All four camera corners hit actual water geometry at every frame. ROOT and owner
viewed the decoded peak frame. This is visual acceptance, not fluid accuracy.

Actual backend is Intel Arc Pro OpenGL 4.5; the NVIDIA device was not used.
Rendering took 31.515 s, including convex validation 0.078 s and context setup
2.437 s; mean frame cost 0.417 s. Cached CLI plus full decode/poster QA took
35.829 s. Separate wider-water display preparation took 100.75 s once; the
central accepted water is bit identical.
Prepared scene SHA256
`1079ae90ec78ea3171204be2996a101caf86015e44771b509b9af82352a3ef4f`
is reused without resampling waves for camera/material iteration.

To reproduce this example's producer-heading camera, use the accepted settings:

```powershell
uv run --extra visualization topwater render artifacts/reference/stronger-buoy-20261024/pyvista-scene-attempt2/scene.npz videos/preview.mp4 --settings reference/wecsim/pyvista-cached-settings.json
```

The earlier prepared scene's `render-settings.json` belongs to the rejected
shadow-enabled style and is preserved only as historical evidence.

Evidence: `pyvista-clip-attempt3-source/proposal.json` SHA256
`5b434785e9daf58e0213e7fbc86f8f4d0d6742bdb1d6587f380e2dc7975d7259`,
`pyvista-qa-attempt3/verification.json` and the prepared scene's
`pyvista-scene-attempt2/preparation.json` under the same stronger reference root.
No real-time or future large-mesh cost is established by this five-second case.
