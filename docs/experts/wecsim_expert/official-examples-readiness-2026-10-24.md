# Official provided examples: executed receipts

2026-10-24. Coordinator bounded assignment: prepare TWO existing official linear
baselines before further custom-sea work. No MATLAB connection, model compilation,
simulation, dependency download or source-library change in this preparation.
All old studies remain preserved. Future canonical case still requires true local
calm arrival, zero body velocity/radiation memory, viscous off and full60s.
Official examples below retain THEIR documented physics and ramp/PTO settings.

## Exact selected versions and assets

Maintained checkout `artifacts/dependencies/topwater-wec-sim` HEAD
`f962103a348d63101ce669c2611f7d7c109755c3`, clean. Official v7.1.0 base
`a23d2993807410fe78cc63546c1ba9d81912024f`. Each copied official asset has the
same Git blob as this base and byte-identical SHA256 to the local checkout;
the preparer accounts for Git checkout newline conversion without editing assets.
Every individual model/input/STL/H5/golden hash and size is in each proposal's
`official_asset_proof`; actual source/library dependencies are frozen separately.

| Case | Official input/model | Native settings | Golden |
|---|---|---|---|
| standard | `examples/RM3/{wecSimInputFile.m,RM3.slx}` | regular H2.5m/T8s,400s,dt.1,ode4,ramp100s | No supplied trajectory golden for this400s example |
| cic | `tests/RegressionTests/RegularWaves/regularCIC/{wecSimInputFile.m,RM3.slx}` | regularCIC H2.5m/T8s,200s,dt.1,defaultode4,ramp100s,CIC20s | Supplied `regularCIC_org.mat`; official full-history heave tolerance1e-5m for each body and relative heave |

Both use supplied `examples/RM3/hydroData/rm3.h5` (4,211,144 bytes), equilibrium
masses, official dry inertias, float/plate STL and PTO damping1,200,000Ns/m,
PTO stiffness0. Default linear hydro (`nonlinearHydro=0`), quadratic/linear drag
zero, Morison off, stateSpace0 and b2b0 are retained. They are NOT closed body-only
energy systems: incident forcing, PTO energy extraction and mechanical constraints
are essential. No removing damping/PTO or substituting our boat coefficients.

Model source links distinguish the standard Floating(3DOF) constraint from the
CIC regression's Translational constraint (heave-only bodies). The latter model
uses old `WECSim_Lib/...` links; the pinned checkout contains that official library
alongside modular libraries. No model migration/edit was made here. ActualR2025b
load/update compatibility is still unexecuted, so report any failure rather than
silently replacing its model/constraint. The original regression runner is preserved
but not invoked: it rebases heave only for comparison and deletes outputs; our
wrapper retains all native output/workspace and repeats its comparisons directly.

## Startup and energy interpretation

Official `initializeWecSim.m` clears convolution persistent functions, reads the
official input, then initializes the model. Native CIC histories start zero;
body initial displacement/velocity defaults remain unchanged. Both inputs use
the native100s wave ramp, not a physically travelling calm-to-arrival packet.
The standard regular mode uses constant-frequency A/B, not CIC. The second case
exercises actual native convolution. None of these facts changes our future
canonical startup requirements.

Harness records actual initial position/velocity/radiation, all six full histories,
native stored applied total/added-mass forces, and native force-sum residuals.
`postProcessWecSim.m` replaces reported added-mass/total force after restoring mass;
the original applied forces remain in `body.hydroForce.hf1.storage`, which we save
separately. CIC resisting power uses worldv/worldomega dot resisting Frad; negative
instantaneous power or component work is NOT itself an error. Cumulative work is
a diagnostic, not a positivity or full fluid-energy acceptance test. Supplied
golden agreement checks implementation/regression, not our BEM accuracy/passivity.

## Archived preparation and commands (subsequently executed below)

Preparer: `reference/wecsim/prepare_official_examples.py`. It already copied fresh
unmodified assets into `artifacts/reference/official-examples-20261024-ready2/{standard,cic}/official`.
Directory layout preserves all official relative paths. Do not rerun against this
existing destination; later attempts require a new explicit directory and archive.
No old result directory is overwritten.

Source helper: `reference/wecsim/run_official_example.m`. Native entry is exactly
`wecSim` from the appropriate official working directory. Sole display changes:
`defaultFigureVisible=off` and official native MCR input override
`mcr.header={'simu.explorer'};mcr.cases={'off'};imcr=1`. They suppress windows;
neither changes ramp, time, memory, constraints, force coefficients nor damping.
Original input remains byte-identical. Official postprocessing is retained.

For EACH case separately, after coordinator allocation and authoritative engine
readiness, use its `proposal.json`: acquire lease, guard `record-claim`, SEPARATE
guard `verify`, then guard `run`. Exact prospective producer commands are:

```
.venv/Scripts/python.exe reference/wecsim/run_matlab.py artifacts/reference/official-examples-20261024-ready2/standard/entry.m --log artifacts/reference/official-examples-20261024-ready2/standard/engine.log
.venv/Scripts/python.exe reference/wecsim/run_matlab.py artifacts/reference/official-examples-20261024-ready2/cic/entry.m --log artifacts/reference/official-examples-20261024-ready2/cic/engine.log
```

One engine/oneCPU, max600s EACH, estimated30-240s each (not measured). No parallel
examples, automatic retry, full test suite or stock checkout/download. Exact terminal,
frozen identity verification, idle-under-lease/separate readback/final release before
the next allocation. Expected grids4001/2001 samples including0; retain native
output workspace plus `audit.mat`/`audit.json` in each official working directory.
Golden comparison saved before its tolerance assertion. Old standard golden is
absent, explicitly: success is native completion/grid/startup checks and diagnostics,
not an invented reference trajectory. Source syntax lint passed; no runtime PASS.

Stopping point: source/readiness reviewed, return proposal paths to coordinator,
then STOP for numerical allocation. No custom sea/packet/pruning/seed work resumes.

Independent non-author SOURCE PASS received from researcher: harness
`72713aa792cd57c0d4095a5917033aea253d2e348ba49f996c52ea4a2ea2a56a`;
standard proposal `c7ae2279f33242aa3d326d99186cc193466a8b2b1da9f47bb74b1371dd4c3b43`;
CIC proposal `4c5aee5d281bd337da9e0d938695ede2cfc5a2a1594c10b2a074fc49ebe600d0`.
Reviewer independently checked all134/138 frozen identities and official-base
asset proofs. Native path/stale-model/complete-grid/selected-wave/finite-force
guards and original golden logic passed source review. This is NOT executed PASS.

## Executed two-job receipt

Coordinator subsequently allocated exactly these TWO SERIAL jobs, unchanged
source/settings. Both completed on sharedMATLAB9168/R2025b, oneCPU, no retries.
Authoritative engine queries proved RM3 absent before each run; each wrapper closed
only its owned RM3 after saving results. Maintained fork remains clean.

| Result | Standard | Official CIC regression |
|---|---:|---:|
| Native elapsed including postprocessing |31.2566057s|28.2803306s|
| Complete grid |0:.1:400,4001samples|0:.1:200,2001samples|
| Initial velocities and radiation forces |exact0, both bodies|exact0, both bodies|
| Applied native force-sum residual |exact0, both bodies|exact0, both bodies|
| Finite positions/velocities/forces |PASS|PASS|
| Combined resisting radiation work |64,757,239.5686495J|23,520,328.3683494J|
| Minimum cumulative work, each body |0J|0J|
| Official golden max heave error |not supplied|[4.156675e-13,1.065814e-13,3.126388e-13]m|

CIC golden errors are float/spar/relative heave over the full200s and pass the
official1e-5m tolerance. Instantaneous CIC resisting power is negative at some
times (minima -6765.693W and -1.473974W): this is not automatically an energy bug.
The standard case uses constant-coefficient radiation, not CIC. PTO dissipation
and constraints remain native; these are open forced systems, not closed-fluid
energy-conservation or all-frequency passivity demonstrations. Neither qualifies
our17-panel boat BEM coefficients or removes the existing native-pressure witness.

All outputs stay in each proposal's original `output` directory: official full
`output/RM3_matlabWorkspace.mat`, `audit.mat`, `audit.json`, native simulation.log
and engine.log. Main audit SHA256:

- Standard MAT `0a9273c4c9593e6f731646802d3804a58701a9d24f6a4bd5917a5e1c2dcb1ad2`;
  JSON `ccfb68bf6934ff0dd4334123c52dc7f80bd876a9c91cbe3f02ff9d264e42f341`.
- CIC MAT `57b4179a177117b11888c79e1273d13b346817648bd27c7a6bf47a6cfbe996fb`;
  JSON `e8386d8b2eedf138c681f5aa7395b2f286b0b5f175f6a8eb61a99741cbe613b9`.

Independent researcher recomputed both full grids, initial values, applied-force
identities/work and CIC supplied-golden comparisons; all134/138 frozen inputs
match. No runtime physical edits/model replacement/fixture recovery needed.

Terminal bookkeeping:

- Standard exactclient33180/guard38440, session67753 terminal0/chunk78b1b2;
  authoritative CIM6b3d0e proves both absent, only known engine docsserver33628.
  Allfrozen/idle-under-lease e1d7f9; SEPARATEreadback51485b; FINALrelease d93e01.
- CIC exactclient34992/guard25872, session70311 terminal0/chunke1c9a2;
  authoritative CIM68b7b4 proves both absent, only docsserver33628.
  Allfrozen/idle-under-lease de7823; SEPARATEreadback3e74ad; FINALrelease5abf98.

No ledger writes after either release. Both numerical jobs terminal; no remaining
allocation. STOP: report results, no automatic custom-case/packet/pruning work.
