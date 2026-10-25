# Native WEC-Sim references

Latest shared-library preview:
[5s square PyVista clip](../../artifacts/reference/stronger-buoy-20261024/pyvista-clip-attempt3.mp4).
H264800square/63frames12.5fps/5.04s; physicaltimes50.00..54.96, same accepted
stronger motion/train. SHA `c80e14001a9e3818291e428cea71378b3e158b44831d084afcc34b67d1ec2b0f`.
Full63decode and ROOTposter review PASS. ActualIntelArcPro cachedrender31.515s;
initial reusable widerDISPLAYdata preparation100.75s.

WEC adapter prepares producer-neutral worldmeshes/time metadata; all drawing is
shared `topwater.visualization`, with no currentboat dimensions/heading embedded.
The reusable command takes two required paths and renders the wholeinput:

```powershell
uv run topwater render examples/scene.npz videos/preview.mp4
```

Only the namedMP4 is the publicdefaultoutput, fixed800square and cadencefromNPZ;
no start/end/size/fps flags. Optional `--settings camera.json` controls style.
The WEC adapter selects the desiredinterval beforewriting its reusablecache.
Prepared actualdata: `artifacts/reference/stronger-buoy-20261024/pyvista-scene-attempt2/scene.npz`.
Use fresh outputpaths and the ordinary serial jobprocedure. [Recorded isolatedenv,
relativepath example and limits](../../docs/experts/wecsim_expert/pyvista-readiness-2026-10-24.md).
Convexhulls are an explicit validatedinput requirement; visualwatermask is not
flooding physics. Finalstandardlighting has optionalshadow/waterlinepasses off
afterliteral artifacts; no customrenderer or changedwaves/bodytrajectory.

The stronger buoy-derived case is complete:60s/3001native samples, source record
46086/2023-12-23T00:40Z, Hm01.657m/Tm027.237s, seed42 and bow182.067deg.
It starts from rest in numerically locally calm water with zero radiation memory,
no temporal ramp and no added viscous damping. Source amplitudes and directional
spreading come from the observation, not a scaled mild realization.

Latest clearer presentation movie:
[reference-presentation.mp4](../../artifacts/reference/stronger-buoy-20261024/presentation-encode-attempt2/reference-presentation.mp4).
H2641280x900,751frames/12.5fps,60.08s; MP4 SHA
`7a9debb31a440b6967fd3340f2fcd6f89072699dd7ecc2e36136d3bca6079aa8`.
Full751decode and actualencoded0/30/peak51.44/60frame review PASS. Selected24deg
camera/10m span uses equalmetres, solidhull, actualcached-contact and noaxes/grid.
Distantwater is a320x225 display pass; hull/contact/neighborhood remain1280.
Visible finitepatch edges/coarsening are rendering limits, not altered physics.
Actual fullmovie render1064s ononeCPU; [receipt and preserved iterations](../../docs/experts/wecsim_expert/presentation-readiness-2026-10-24.md).

The previous12m isometric movie remains available unchanged:
`artifacts/reference/stronger-buoy-20261024/render-encode-attempt1/reference-isometric.mp4`,
SHA `42c8ca1c80e961784200d076faff7e9a34a67c782120019e606bbf553ec004e2`.
Both use the same native motion and prescribed incoming field.

Peaks surge/sway/heave[.562,.837,.759]m and roll/pitch/yaw[6.619,2.678,1.325]deg
meet the requested stronger-motion target. Net radiation work32971.040J is
positive throughout cumulatively; returned instantaneous work4074.896J and
incident work46582.193J are recorded separately. No complete energy closure or
general passivity proof follows. [Case evidence](../../docs/experts/wecsim_expert/stronger-buoy-readiness-2026-10-24.md).

This is a provisional fixed-mean linear reference, not mesh/time/nonlinear
accuracy qualification. A cached local-water-plane proxy gives~3.10deg maximum
relative heel,~2.88deg trim and0.090m relative vertical offset. Existing3deg
flatwater curves show~7.3% righting-secant and15.3% tangent-stiffness departures;
combined curved-water wetting/flooding and small geometry error are unconfirmed.
[Applicability context](../../docs/experts/wecsim_expert/stronger-relative-water-proxy-2026-10-24.md).

The current mild buoy-derived case is complete:60s/3001native samples, locally
calm incoming-wave startup, zero boat velocity/radiation memory, no temporal ramp
or added viscous damping. Source observation46086/2023-03-19T08:40Z has
Hm0.535m/Tm027.427s; its derived finite train and both physical-coordinate frames
are retained separately. Mean bow201.308deg faces the dominant incoming energy.

Local video: `artifacts/reference/calm-buoy-20261024/render-encode-attempt1/reference-isometric.mp4`.
H2641280x900,751frames/12.5fps,60.08s including the native60s endpoint hold.
The close12m view uses equal physical metres and the same prescribed incident
train as the force input. Decoded0/30/60frames and all751frame metadata pass review.
MP4 SHA `205400a829eef3ffb9b1fdf83a359299dc05e4edf15f1f051f3765fb5ca7e8f9`.

Peak surge/sway/heave are[.238,.274,.240]m and roll/pitch/yaw[1.481,.663,.617]deg.
All3001actual excitation samples match the prepared forces exactly. Radiation
memory absorbs net3271.532J (temporary returned work309.916J, cumulative minimum0),
against3492.514J incident work. This is not a complete fluid-energy budget or
mesh/time/nonlinear accuracy qualification. [Case evidence](../../docs/experts/wecsim_expert/calm-buoy-arrival-readiness-2026-10-24.md).

The stronger case and clearer presentation are complete; the mild case remains
immutable. Current scope and actual execution state are in the handoff/shared ledger.

The accepted minimum17-panel case has501 native six-DOF samples covering0..10s
from rest at dt0.02s. The same frozen dispersive incoming field drives the
excitation and visualization. The centered isometric video and literal encoded
frames passed independent review, alongside the native motion. This establishes a minimum-case
workflow, not temporal or mesh motion accuracy.

Read the [handoff](../../docs/checkpoint-handoff.md),
[physical contract](../../docs/experts/shared-contract.md) and
[serial execution procedure](../../docs/experts/coordination.md) before running.
Each numerical or rendering job needs a bounded assignment, frozen manifest,
one owner and the shared atomic lease. These scripts do not claim/release it.
Preserve accepted outputs and failed attempts; many scripts intentionally reject
existing output directories/logs. The commands below reproduce the recorded
pipeline with fresh artifact paths. A new run must use a separately assigned
fresh output root or reviewed path adaptation, not overwrite these artifacts.

## Seeded strong-motion prototype: superseded startup

Latest user requirements for the final60s validation case are TRUE locally calm
water with waves initially upstream, body rest/zero velocity, zero radiation
memory, viscous damping off and a physics/accuracy-justified
incident band. The nonzero-history prototype below is historical evidence,
not that final requested startup. No final new band has been selected.

The preserved REF-002 prototype uses observed record `ndbc-46086-20250105T1140Z`, generator seed42 and
all384 directional components without amplitude or phase fitting. Hm0 is2.270m
and Tm02 is9.296s. The initial mean bow is253.218deg CCW from global east,
90deg from the incoming-aligned default to excite roll. This is a physical
heading/frame change, applied once in the force inputs and display.

The native60s case has3001samples at dt0.02 and40s CIC memory. Peak motions are:

| Surge | Sway | Heave | Roll | Pitch | Yaw |
| ---: | ---: | ---: | ---: | ---: | ---: |
|0.368m|1.549m|1.520m|5.127deg|1.723deg|2.200deg|

Roll exceeds5deg only at the final60.00s sample; this is a peak criterion,
not sustained5-10deg statistics. Full startup and all six histories are retained.
Pose/velocity and physical radiation prehistory come from a deterministic
mean-linear particular solution. All post0 motion is native and free. Native
acceleration-feedback startup is unchanged; its nonzero initial residual is
reported, rather than cancelled. Independent native fixture and complete
actual-state/force/history/wrench reviews pass.

This remains a **provisional17-panel model trajectory**. Its incident-water geometry assessment was intentionally stopped under the new
pruning/startup steering. Radiation passivity, mesh/time accuracy
and nonlinear equivalence are unqualified. The rejected held-from-rest sea run
and calm-water control are preserved as diagnostics, not additional canonical
sea cases. Topwater solver implementation remains held.

Source entry: [run_seeded_jan90_case.m](run_seeded_jan90_case.m).
[Case, initialization, commands and evidence](../../docs/experts/wecsim_expert/seeded-native-case-proposal-2026-10-24.md)
records the frozen generator/native/BEMIO dependencies and fixture prerequisites.
Actual outputs are in
`artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1`.
The generated renderer preserves the full bank and native motion in an equal-metre
10m XY tracking view. Previews pass; encoding was never allocated. Do not finalize this prototype as
the newly requested zero-memory/pruned validation case.

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

## Rendering and evidence

Preview before encoding, under separately allocated artifact jobs:

```powershell
.venv/Scripts/python.exe reference/wecsim/render_reference.py --case artifacts/reference/minimum-20261024/case-attempt1 --output artifacts/reference/minimum-20261024/case-attempt1/render-attempt3 --preview-only
.venv/Scripts/python.exe reference/wecsim/render_reference.py --case artifacts/reference/minimum-20261024/case-attempt1 --output artifacts/reference/minimum-20261024/case-attempt1/render-attempt3
```

The fixed orthographic isometric camera uses a12m square centered on the full
mean hull planform. Actual metres and native motion are unchanged; mm/degree
telemetry makes small motion visible without amplification. The transparent
surface is the exact prescribed incident spectral field, explicitly labelled;
it is not a computed scattered/radiated fluid surface.251frames at25fps encode
10.04s because the native t10 endpoint is held for one frame. Original draft
videos remain at their original paths; final delivery uses its unique directory.

Evidence:

- [Accepted input/inertia](../../docs/experts/wecsim_expert/accepted-input-2026-10-24.md)
- [Imported history/interface](../../docs/experts/wecsim_expert/imported-excitation-2026-10-24.md)
- [Independent native motion review](../../docs/experts/capytaine_expert/native-motion-review-2026-10-24.md)
- [Native run and artifact record](../../docs/experts/wecsim_expert/minimum-reference-2026-10-24.md)

The stock1e-7s residual-acceleration Transport Delay explains initial apparent
full-mass/acceleration imbalance and fixed-step forward extrapolation. Independent
all501 native Solid Newton-Euler/force-bus checks passed; no correction was
justified. This numerical consistency is not a mesh/time convergence result.
Topwater implementation, optimization and automatic refinement remain held.

## Separately bounded follow-up

Selecting a sea with larger visible motion requires nonlinear confirmation;
scaling this linear reference alone is insufficient. The selected WEC-Sim7.1
path rejects `elevationImport` with `nonlinearHydro=1`, so preserving physical
upstream arrival in that follow-up needs a feasibility decision. The accepted
minimum input and result remain frozen; this follow-up does not unhold Topwater.
