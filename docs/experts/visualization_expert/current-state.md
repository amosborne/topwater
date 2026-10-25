# Visualization expert: current state

User requested a dedicated PyVista/visualization expert and a shared Topwater
render workflow accepting standard inputs from either WEC-Sim or Topwater.
Different convex hulls are supplied through data; convexity is a validated mesh
requirement. Standard PyVista/VTK handles graphics. Physics remains producer-owned.

Initial workflow and preview COMPLETE: fixed square800x800, source50.00-54.96s,
63frames/12.5fps/5.04s, continuousblue actualwater covering the fullviewport,
litwarm/ivoryhull, perspective, unchangedphysicalmotion. ROOT viewed peakposter.
Watch `artifacts/reference/stronger-buoy-20261024/pyvista-clip-attempt3.mp4`,
SHAc80e14001a9e3818291e428cea71378b3e158b44831d084afcc34b67d1ec2b0f.
Full63decodePASS. No full60s encode or follow-on work is allocated.

Final CLI: `uv run --extra visualization topwater render INPUT_PATH.npz OUTPUT_PATH.mp4`.
Both paths are required positional arguments, absolute or relative. Always render
the whole input at its uniform timestamp cadence; no start/end/size/output flags.
Only the exact requested MP4 survives normal rendering. Camera/light/material
settings may be supplied as JSON; diagnostics return on stdout and scratch is
cleaned. Examples prefer relative paths, e.g `data/scene.npz videos/preview.mp4`.
Temporary reference QA is allowed while tuning.

Implementation: [shared module](../../../src/topwater/visualization.py),
[CLI](../../../src/topwater/cli.py). `SceneSequence` stores time_s,
hull_points[T,N,3]/hull_faces[F,3], water_points[T,M,3]/water_faces[K,3], metadata.
Metres/right-handed-z-up, fixed triangular connectivity per sequence, uniform
timestamps; producers normalize geometry. Prepared NPZ is a named-array archive.
OnePlotter/reusedactors; standardVTKfilters/materials, optional/lazygraphicsdeps.
Convex validation uses bounded storage rather than a dense face-by-plane matrix.
Water cells must resolve hull intersection for standardVTK clipping to mask it.

Initial library/schema implementation is owned by `/root/visualization`
(topwater_expert); WEC-Sim adapter by `/root/wecsim_reference`. Transfer ownership
explicitly at a stable checkpoint; do not duplicate their source.
This new role may first review the contract and rendering approach read-only.

Isolated visualization environment installed under ignored workspace artifacts:
`artifacts/dependencies/pyvista-venv/Scripts/python.exe`, PyVista 0.49.1,
VTK9.7.1/SciPy1.18.1. Actualbackend IntelArcPro/OpenGL4.5 hardware. Finalcached
render31.515s (validation.078s/context2.437s/frameavg.417s); wholeCLI+QA35.829s.
Initialsame-field widerDISPLAYdataprep100.75s is cached separately and not repeated
for styling. 13focusedcontract/package/CLItestsPASS and a secondconvexshape2frame
CLI smoke PASS inattempt2. Failedtemp-fixtureattempt1 and visuallyrejected
shadow/waterline/discretepaletteattempt2 remain preserved underoriginalpolicies.
Defaultshadows/waterline nowoff, standardcontinuousBlues_r/ambientlighting.
[Implementation evidence](../topwater_expert/visualization-workflow-2026-10-24.md).
No numerical/render job allocated to this role. Role TOML schema validation
passed; model/settings inherit the parent. Native spawn attempt for a bounded
read-only contract review returned `unknown agent_type 'visualization_expert'`.
This session does not expose the newly added role. No fallback agent or separate
CLI coordinator was started; existing Topwater implementation owner continues.
