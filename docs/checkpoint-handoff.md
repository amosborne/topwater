# Checkpoint handoff

## Active goal: canonical Capytaine -> BEMIO -> WEC-Sim,10s movie and mesh variation

The USER replaced the old goal and expressly excludes ALL Topwater-specific
code work. Read docs/capytaine-wecsim-workflow.md and current AGENTS override.
Goal ACTIVE/INCOMPLETE. No10s boat simulation or movie exists yet.

User requirements:
- Minimum hull mesh supplies geometry. Solvers may refine that same surface.
- Boat initially at equilibrium/rest in locally calm water, rho1025,g9.81,
  infinite depth. Full dry inertia/products/CG must be preserved.
- Waves START ROLLING IN from upstream. NO temporal amplitude ramp, global
  wave appearance under hull, or suddenly switched harmonic loads as surrogate.
- Canonical mild directional sea. Initial waves may be upstream; excitation
  and the rendered incident field must share the same actual arrival input.
- Capytaine frequency-domain outputs -> official BEMIO -> WEC-Sim six-DOF
  time-domain dynamics. Do not build another custom reference integrator.
- First produce the working minimum-mesh10s result/isometric video. Then compare
  actual six-motion histories as solver meshes refine. Intermediate added-mass,
  damping or coupling coefficients are diagnostics, NOT standalone gates.
  Preserve old failed studies under their original policies; never relabel them.

Verified setup:
- MATLABR2025b25.2, Simulink/Simscape/SimscapeMultibody installed and all3license
  tests returned1. Inventorysession5520terminalexit0 (default sandbox startup
  initially lacked writable preferences; normal approved MATLAB run succeeded).
- OfficialWEC-Simv7.1.0 cloned at artifacts/dependencies/wec-sim-v7.1.0,
  commita23d2993807410fe78cc63546c1ba9d81912024f. Clone40953terminalexit0.
- MATLABsession53356terminalexit0 loads WECSim_Lib, instantiates simulationClass
  and waveClass(elevationImport), resolves official readCAPYTAINE. Source path
  added in batch only; no global MATLAB path was saved. Setup evidence:
  docs/capytaine-wecsim-environment.json. This is NOT a boat simulation.
- Capytaine3.0.0 installed; FFmpeg9.0.1available. Matplotlib3.11.2 installed in
  the analysis venv only, no projectdependency or Topwaterruntime changes.
- Neutral geometry/data: reference/capytaine/canonical-boat.json, full minimum
  hull plus exact17wet triangles, full dry M, world CG/pose and sea. Exported
  from exact existing frozen input/NPZ (hashes included), no solver imports.
  Arrival mapping remains explicitly prospective: freeze it before actual run.
- The brief custom Python dynamics experiment has been set aside under
  reference/capytaine/history and is NOTused. Four analytic math checks passed;
  they establish no boat accuracy. Do not extend that prototype.
- Old b64cb63 quadrature job session36762 was USER-INTERRUPTED terminalexit1 after
  G24, with G20 reused/G28 unfinished. Local partial artifacts remain at
  artifacts/validation/fixed-field-2772-continuation-001.*, unchanged. Not a
  numerical failure. All older coefficient studies remain historical.
- No live numerical jobs. No src/topwater or existing tests were modified.

Important WEC-Sim source findings:
- examples/BEMIO/CAPYTAINE/sphere/main.py and bemio.m show native netCDF plus
  legacy hydrostatics export, readCAPYTAINE/radiationIRF/excitationIRF/writeBEMIOH5.
- waveClass elevationImport interpolates an imported elevation history; with
  rampTime0 it imposes no ramp. bodyClass.userDefinedExcitation applies the
  official NONCAUSAL excitation-IRF convolution relative to the wave reference
  point, choosing direction from native hydroData. Do not mistake zero local
  elevation alone for zero prearrival body forcing. Check padding and phase.
- Imported-wave markers are disabled; built-in imported surface FFT reconstruction
  warns qualitative. Use a consistent incident-wave visualization and label it;
  do not pretend it is a total scattered/radiated CFD field.
- bodyClass supports inertiaProducts[Ixy Ixz Iyz] as well as diagonal inertia.
  Preserve product signs/conventions through the full body model.
- Independent setup review passed; no checkpoint blocker. Before launch: use
  per-direction imported components (one scalar with one heading loses the sea),
  explicitly set simu.rampTime=0 (default100), verify missing hydrostatic fields
  are not replaced with importer zero CG/CB/volume defaults. CSV padding beyond
  simu.time is discarded before convolution: pad ACTUAL timeline. Preserve the
  first10s from the true rest initial condition, never crop to an already-moving
  later window. See workflow doc for these resolved prospective requirements.

**Next bounded action:** FIRST end-to-end minimum-mesh boat run and10s isometric
movie using Capytaine/BEMIO/WEC-Sim. Scope: generate minimum-mesh native Capytaine
netCDF/hydrostatics; officialBEMIO H5; single freely-floating six-DOF Simulink
model, no PTO/mooring/invented damping; freeze incoming spatial arrival input
and its force/rendering mapping; run10s from rest, retain all startup. Archive
actual code/inputs/dependency hashes BEFORE each numerical launch; declare
frequency band/spacing/IRF/time step/resources prospectively, serial jobs,
immutable producers until exact process terminal. Review startup, full inertia,
phase/sign/units and surface interpretation before promotion. Deliver actual
raw native data, histories, plots/movie, reproducible commands and honest limits,
reviewed evidence, checkpoint commit/push, issue30 and fresh handoff; THENSTOP.
No coefficient or extra-quadrature qualification campaign. No automatic large
mesh sweep: after first working result, characterize motion variation over at
least3 actual panel resolutions and relevant time/frequency/IRF refinements,
all judged in six motion units. The full video/mesh-study goal remains active
until those deliverables exist and are checked. No Topwater work.
