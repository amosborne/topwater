# Standard visualization workflow

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
