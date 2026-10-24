# Minimum-mesh native reference

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
