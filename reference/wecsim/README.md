# Historical WEC-Sim reference workflows

These scripts, commands and receipts describe the workflow at commit
[`7570f94`](https://github.com/amosborne/topwater/tree/7570f94/reference/wecsim).
The scene-cache/two-path render interface below is retired. In particular,
`run_full_reference.py` expects that former CLI and its JSON stdout; it is
incompatible with the current directory-based CLI. The scripts and evidence
remain preserved for their original checkpoints.

For the current workflow, use `examples/simple-boat` and the
[generate-sea](../../docs/cli/generate-sea.md), [solve](../../docs/cli/solve.md)
and [render](../../docs/cli/render.md) commands. Current rendering reads the
self-contained `result.npz`; this historical scene cache is not that format.

## Recorded 2D result at the historical checkpoint

[Watch the full60s topdown video](../../artifacts/reference/stronger-buoy-20261024/matplotlib-full60-attempt1.mp4).
SHA `15510982f6cef3a246f1129222e0525052215a06cfa2e0073af3c53bf803f23b`.
H264800square/751frames12.5fps/60.08s; native0..60inclusive. Full751decode and
0/30/51.44/60stills passed, ROOT0/peak accepted. Same acceptedstrongerbuoy WECmotion,
not a newsimulation. [Finalreceipt](../../docs/experts/wecsim_expert/full60-scene-video-2026-10-25.md).

```powershell
uv run topwater render examples/scene.npz videos/full60.mp4
```

INPUT andOUTPUT are two REQUIREDPOSITIONALpaths. Render WHOLEpreparedNPZ atits
savedcadence, fixed800square, onlynamedMP4defaultoutput. No --extra/start/end/
size/fps/output flags. GenericTopwaterMatplotlib code knows no WEC/hull constants.

Actualreusableinput: `artifacts/reference/stronger-buoy-20261024/scene-full60-attempt1.npz`,
SHA `884667cb0a43f55ec86303c38fae303a76c58a8869c547b1d3290f70b4f330bc`.
Original751hull/time/central65x65waterXYZ are BITEXACT; only1550missingbordernodes
perframe evaluated fromthe SAMEincidenttrain (75x77,~104MBwater). NativeHUDqexact.
Prep40.657s/CLIwall148.792s/render145.75s/wholejob197.911s. Futurestyleiterations
reuseinput, not waves. Peak/cadence/physics unchanged; no motion amplification.

The producer is `prepare_full_reference_scene.py`; `run_full_reference.py` wraps
producer, cleanpublicCLI andexplicitreferenceQA. Existingoutput paths failclosed;
use a fresh reviewedconfig/sourcefreeze and ordinary serialclaim procedure for
any rerun. Read [handoff](../../docs/checkpoint-handoff.md),
[contract](../../docs/experts/shared-contract.md) and
[coordination](../../docs/experts/coordination.md) first. No automaticjob isqueued.

USERapprovedactive3Drendering cleanup is complete: obsoletecustom/PyVista render
helpers/configs removed afterimportaudit; active2D has no dependency on them.
HistoricalMP4s/frozenproducers/evidence remain. Physical3Dhydrodynamicgeometry,
solver code and nativeWEC-Sim fork visualization are unchanged. Currentpipeline
is2Donly. Older3Devidence links describe archivedresults, not currentcommands.

Physical reference: collected46086/2023-12-23T00:40Z Hm01.657m/Tm027.237s,
unscaledseed42bank, physicalincomingtrain; rest/zero radiationmemory/no ramp/
viscous0, fullnative60s sixDOF. Original17H5/CIC10/dt.02. Peakroll6.619deg.
Implementation passed; linear/mesh/time/curvedwaterwetting/flooding accuracy
remains provisional. Displayprescribedwater/contact is not a floodingmodel.
[Caseevidence](../../docs/experts/wecsim_expert/stronger-buoy-readiness-2026-10-24.md).

## Dependencies and engine

This machine used project `.venv` Python/Capytaine3.0.0, MATLAB R2025b with
Simulink/Simscape/Multibody, NumPy/SciPy/xarray/threadpoolctl and Matplotlib, plus
FFmpeg/ffprobe on PATH. Exact versions/source hashes are in job manifests.
Use the maintained WEC-Sim project branch pinned in
[wecsim-fork.json](../../docs/experts/wecsim-fork.json), currently
`f962103a348d63101ce669c2611f7d7c109755c3`. Its reviewed fixes are already applied;
`fix_*.m` files are recorded migration/regression tools, not startup steps.

`run_matlab.py` attaches only to the existing `CAPYTAINE_REFERENCE` shared engine.
It never starts another MATLAB. Reconcile actual process/session identity first;
restricted engine-cache visibility is not proof the engine is absent. The tested
hidden bootstrap uses this inline PowerShell command after reconciliation:

```powershell
$matlabProcess = Start-Process -FilePath matlab -ArgumentList '-nosplash','-nodesktop','-r',"matlab.engine.shareEngine('CAPYTAINE_REFERENCE');set(0,'DefaultFigureVisible','off');" -WindowStyle Hidden -PassThru
```

`start_reference_engine.ps1` records that body. This host blocks its direct
`powershell -File` invocation by execution policy; the inline command above was
the working method. Record launcher/actual engine PID, start time and host under
the lease. `engine_query.py` is a read-only inventory, not numerical validation.

## Recorded pipeline

From `C:/Users/alex/topwater`, with execution allocated separately for each job:

```powershell
.venv/Scripts/python.exe reference/capytaine/export_case.py --output artifacts/reference/minimum-20261024/native-attempt2
.venv/Scripts/python.exe reference/capytaine/correct_hydrostatics.py --source artifacts/reference/minimum-20261024/native-attempt2 --output artifacts/reference/minimum-20261024/native-gl2
.venv/Scripts/python.exe reference/capytaine/check_hydrostatics.py --corrected artifacts/reference/minimum-20261024/native-gl2
.venv/Scripts/python.exe reference/wecsim/run_matlab.py reference/wecsim/prepare_bemio.m --log artifacts/reference-bemio.log
.venv/Scripts/python.exe reference/wecsim/incident_packet.py --output artifacts/reference/minimum-20261024/packet-attempt3
.venv/Scripts/python.exe reference/wecsim/extend_packet_region.py --output artifacts/reference/minimum-20261024/packet-full-planform
```

The producer exports native hydrodynamics unchanged on17 triangles. The separate
native GL2 correction restores exact degree-two hydrostatic moments and provides
the importer-compatible `body_name` alias. Official BEMIO imports/writes H5;
exact native CG/CB/volume/K/Ainf supersede legacy text rounding with reviewed
normalization and one phase conversion. Excitation uses+-20s history/future
support; radiation memory is10s. Full dry inertia is supplied from native MAT.

The frozen packet preserves canonical carriers0.1/0.25/0.5Hz and20/10/5mm,
headings0/45/90degrees, and phases. Focus times100/11.5/200s place it upstream.
Input calm/force/arrival gates and full physical-hull planform coverage passed
independent review. The historical fixed-order `screen_packets.py` study is
preserved; reproduction of the accepted input does not select a new packet.

Serialize the accepted retained-force array for native runtime comparison
(format conversion only, no altered samples), then run the cold case:

```powershell
.venv/Scripts/python.exe -c "import numpy as np;from scipy.io import savemat;p=np.load('artifacts/reference/minimum-20261024/packet-attempt3/startup-forces.npz');savemat('artifacts/reference/minimum-20261024/expected-retained-force.mat',{'time':p['time'],'force':p['force_convolution']})"
.venv/Scripts/python.exe reference/wecsim/run_matlab.py reference/wecsim/run_reference_case.m --log artifacts/reference/minimum-20261024/case-attempt1-run.log
```

Cold construction creates a native Global Reference Frame template, executes
official preprocessing/variant/bus setup and radiation-persistent reset, then
instantiates the native free6DOF joint/body and reloads/updates the model. There
is no simulated warmup. The mean-world STL is shifted by mean CG once for the
local File Solid; the equilibrium quaternion and world inertia are not reapplied.
Native angular velocities are world omega, not Euler rates. Native orientation
columns are radians with `Rz(yaw)*Ry(pitch)*Rx(roll)` reconstruction.

Runtime checks read actual H5 metadata, reconstruct compiled full dry+added mass,
retain absolute world CG and raw60-column bus, verify initial displacement/
velocity/radiation force, and compare all501x6 signed applied excitation samples
to the accepted retained history. Native outputs are exported directly; this
case does not use the stock scalar-wave postprocessor. No PTO, mooring, empirical
damping, prescribed trajectory or restricted motion is added.

## Rendering evidence and preserved history

Only shared2D publicCLI above is active. Priorcustom3Drendercommands were removed;
exacthistoricalproducer snapshots remain with artifactmanifests and git history.
Preserved5sproof: `matplotlib-clip-attempt1.mp4`, data/HUDidentity in
[2Devidence](../../docs/experts/wecsim_expert/telemetry-2d-readiness-2026-10-24.md).
Prior3Dhistory: [PyVista](../../docs/experts/wecsim_expert/pyvista-readiness-2026-10-24.md)
and [presentation](../../docs/experts/wecsim_expert/presentation-readiness-2026-10-24.md).
No upstreammessages, physicsrefinement or Topwater solver work follows automatically.
