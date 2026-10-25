# WEC adapter to shared PyVista rendering: completed preview

Updated2026-10-24. ROOT accepted the final poster. WEC owner supplies only
accepted-data preparation/reference tooling; generic rendering/API/CLI lives in
Topwater library and is owned by visualization. No further WEC jobs/code edits.
Executionowner finalizes terminal/lease evidence in the shared ledger.

## Result and reusable command

[Five-second square preview](../../../artifacts/reference/stronger-buoy-20261024/pyvista-clip-attempt3.mp4)
SHA c80e14001a9e3818291e428cea71378b3e158b44831d084afcc34b67d1ec2b0f.
H264800x800,63frames at12.5fps,5.04s. Saved physical times50.00..54.96s,
includingpeakroll51.44s; full63frame ffprobe decode and ROOTposter review PASS.
ReferenceQA: `pyvista-qa-attempt3/verification.json` and `poster.png` under same
stronger artifactroot. Actual OpenGL device IntelArcPro, cachedrender31.515s
(35.829s reference runner plusQA). Initial displaydata preparation100.75s.
No GPU speed claim is inferred from installedhardware alone.

Public generic CLI, relativepaths preferred:

```powershell
uv run topwater render examples/scene.npz videos/preview.mp4
```

Both paths are REQUIRED POSITIONALS. Render the WHOLE prepared NPZ at its saved
uniform cadence to exactly the named MP4, fixed800square. No start/end/size/fps
flags and no implicit outputfilename. Optional `--settings camera.json` is for
styling. The default public renderer writes onlyMP4; explicit reference tooling
captures stdoutJSON/timing/backend and decodedposters separately in referenceQA.

On this machine the isolatedenv is `artifacts/dependencies/pyvista-venv`:

```powershell
$env:UV_PROJECT_ENVIRONMENT = (Resolve-Path artifacts/dependencies/pyvista-venv).Path
$env:UV_CACHE_DIR = "$PWD/artifacts/tools/uv-cache"
uv run --no-sync topwater render artifacts/reference/stronger-buoy-20261024/pyvista-scene-attempt2/scene.npz videos/wec-preview.mp4 --settings reference/wecsim/pyvista-cached-settings.json
```

Use a fresh outputfilename; scripts preserve oldoutputs. Renderjobs still require
bounded coordinatorauthorization/frozenidentity/serialclaim procedure.

## Producer adapter and cached scene

`prepare_pyvista_scene.py` selects63times BEFORE serialization, reads accepted
GLOBAL hullposes/faces and the same frozen incident coefficientarchive. Wider
96m DISPLAYONLY wave sampling uses the existing train.profiles API for those
exacttimes; original65x65nearwater nodes/heights are BITIDENTICAL to accepted
cache. No new sea, force, motion, hull dimensions or simulation is generated.
Shared scene arrays: time_s, hull_points(T,N,3), hull_faces(F,3),
water_points(T,M,3), water_faces(K,3), units/worldframe metadata.
Preparedscene reusable acrossstyleiterations:
`pyvista-scene-attempt2/scene.npz`, SHA1079ae90ec78ea3171204be2996a101caf86015e44771b509b9af82352a3ef4f.
Camera/materials/clip/datareferences are inputs/exampleconfig, not hull-specific
constants in the genericrenderer. A differentconvexsyntheticbox exercised the
same publicCLI with differenttopology/aspect, without a physicalsimulation.

StandardPyVista/VTK perspective/lighting/transparency/surfaceclipping only;
no bespoke rasterizer/shader. FinalcontinuousBlues_r/standardambient configuration
has shadowsFalse/waterlineFalse because those optionalVTKpasses caused literal
artifacts in attempt2. All camera corner rays hit actualwater in everyframe.
Convexphysical/render hullrequirement is validated; an openrim's convexclosure
is a visual watermask, not a flooding/wet-body model. Water is prescribedincident
field, not scattered/radiated fluid. Underlyinglinear/mesh/geometry applicability
remains provisional; visualization does not change that qualification.

## Tools and preserved attempts

Workspace-only Python3.12.8/PyVista0.49.1/VTK9.7.1/NumPy2.5.4/SciPy1.18.1;
imageio2.38.1/imageio-ffmpeg0.6.0. Fullversions in scopedrequirements record.
Simulationvenv/globalPython unchanged. EditableTopwater exposes publicCLI;
workspaceUV_CACHE_DIR avoids sandboxusercache permissions. No ParaView install.

- Attempt1:9testsPASS/4pytestuser-tempERROR, beforefield/GPU. Preserved.
  Exactterminal/releaseb4ac49; minimalworkspace--basetemp runner correction only.
- Attempt2: prepared validcache100.75s, actualIntelArcPro clip33.89s and distinct
  syntheticmeshroute; ROOTliteralREJECTshadow/waterline/blackinterior artifacts.
  Exactterminal/release514ac2; no physics defect inferred.
- Attempt3: cached-only standardCLI rerender, proposal5b434785e9daf58e0213e7fbc86f8f4d0d6742bdb1d6587f380e2dc7975d7259,
  no fieldpreparation/tests/install/simulation. 31.515s cachedrender, full63decode
  and ROOTposterPASS. No more WEC source/job changes authorized.
