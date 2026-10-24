# Canonical Capytaine -> BEMIO -> WEC-Sim workflow

The user replaced the previous goal with: establish the end-to-end canonical
workflow, deliver an isometric video of a10s simulation beginning with a boat
at rest in calm water, and characterize the motion's variation with mesh.
Capytaine and WEC-Sim are the reference toolchain. No Topwater solver work.

## Physical case and startup

Use the existing24ft boat, full dry mass/inertia including products, equilibrium
CG/pose,1025kg/m3 water,9.81m/s2 gravity, infinite depth and zero forward speed.
Neutral input reference/capytaine/canonical-boat.json includes the minimum FULL
hull and its exact frozen17-triangle wet surface. Refinement changes numerical
panels, not physical geometry, mass or inertia. Hull details and equilibrium are
shared data; importing or modifying the Topwater runtime is excluded.

Keep the mild three-component sea (amplitudes20/10/5mm, frequencies0.1/0.25/0.5Hz,
propagation headings0/45/90deg and phases0/pi2/pi) as the carrier specification.
The incident wave TRAIN must begin upstream and travel toward the boat. Initially
the boat is at equilibrium with zero linear/angular velocity and the local
water is calm. No waves suddenly appearing underneath the boat, no global force
step, no arbitrary temporal ramp. Startup transients are part of the10s output;
settling within10s is not required. Arrival geometry/delays must be frozen before
the actual run and identical in excitation and rendered water. A leading edge
adds frequency content beyond the three carrier frequencies: use the transient
input/excitation path rather than just switching harmonic RAO forces.

WEC-Sim7.1.0 supports elevationImport and a non-causal excitation-IRF convolution
relative to its wave reference point; spatial propagation and prearrival padding
must therefore be checked explicitly. An elevation trace by itself is not a
validated traveling field. Body excitation includes diffraction plus incident
pressure; preserve heading and phase conventions through BEMIO. Built-in imported
wave visualization is documented in source as qualitative. Use a consistent
incident surface for the movie, label it incident-wave visualization, and do not
present it as a resolved total scattered/radiated fluid surface.

## Workflow and deliverables

1. Capytaine on the minimum solver mesh: native hydrodynamics and hydrostatics,
   source/input/dependency identity, no intermediate coefficient accuracy gate.
2. Official BEMIO readCAPYTAINE, radiationIRF/excitationIRF and writeBEMIOH5.
   Preserve full six-mode matrices, density/depth/CG and dry inertia. Read raw
   outputs if something is wrong; do not resume blanket coefficient convergence.
3. WEC-Sim single freely floating six-DOF body with no PTO, mooring, invented
   damping or motion constraints. Start at rest; incoming-wave elevation input,
   causal radiation memory, sufficient excitation prehistory/future padding.
4. Save all six10s motion histories, incident signal/arrival metadata, raw native
   outputs, inputs/versions, physical interpretation and an isometric video.
5. AFTER the first working workflow, rerun the same case on progressively finer
   solver meshes and compare full time histories. At least three actual meshes
   characterize a refinement sequence; do not call two-level agreement convergence.
   Report maximum per-axis motion differences in metres/degrees, including
   startup, with time-step/frequency-band/spacing/IRF sensitivity in the SAME
   motion units. Do not gate individual coefficient or coupling errors.

Final deliverable: working reproducible10s video and numerical six-motion
mesh-variation characterization. The previously stated1cm/0.1deg targets inform
reference motion uncertainty; they do not imply a new coefficient threshold.
State uncertainty and limits; a visualization alone does not prove motion
convergence. Small sea amplitude does not by itself establish the linear regime.

## Current bounded checkpoint: establish the approved toolchain

Scope: neutral input export, verified MATLAB prerequisites, pinned official
WEC-Sim release/library/Capytaine importer, inspect the transient input path,
record the new workflow and exact next action. No hydrodynamic job, boat
simulation, mesh study, video or Topwater code changes at this setup checkpoint.
Deliver the setup evidence, reviewed handoff, commit/push and issue30 update;
then stop this bounded setup. Overall goal remains active/incomplete.

Verified MATLAB R2025b with installed/licensed Simulink, Simscape and Simscape
Multibody; official WEC-Sim v7.1.0 commit a23d2993807410fe78cc63546c1ba9d81912024f
is downloaded under artifacts/dependencies/wec-sim-v7.1.0. MATLAB successfully
loaded WECSim_Lib and instantiated simulationClass and waveClass(elevationImport).
It resolves readCAPYTAINE from that pinned release. Setup checks are not a boat
simulation. Capytaine3.0.0 is installed. FFmpeg9.0.1 is available. Matplotlib
3.11.2 was installed only in the analysis environment, not project dependencies.

The abandoned custom Python time-domain prototype is preserved under
reference/capytaine/history, with four passing mathematical checks. It is NOT
the approved time-domain runtime. Use WEC-Sim. Old b64cb63 G24/28 coefficient
continuation was interrupted by the user (exactsession36762 terminalexit1 after
G24); preserve its local partial evidence, never call it a numerical failure.

Next bounded work: produce the FIRST minimum-mesh canonical Capytaine netCDF,
convert it with official BEMIO, build a freely floating WEC-Sim model and freeze
the spatial arrival input after inspecting its native transient handling. Run
the actual10s case and make the first isometric video. Before native launches,
archive actual producer/physical inputs and declare frequency/IRF/time/resource
controls. Jobs serial; no edits to their producers until their exact handle is
terminal. Review physical startup, signs/phase/units, full inertia and arrival
visualization before promotion. Stop at first reviewed end-to-end result and
fresh handoff; THEN begin the separately bounded motion-based mesh study.

Independent setup review passed. Resolve these concrete native-code integration
details BEFORE the first boat launch:

- Use per-direction imported components. waveClass elevationImport generates a
  scalar elevation without direction, while bodyClass selects a directional
  excitation IRF; a single summed trace with one heading loses the canonical sea.
- Padding must exist in the ACTUAL convolution timeline. waveElevUser first
  interpolates onto simu.time, so extra CSV rows alone do not pad its subsequent
  conv(...,'same'). Include required future padding in the simulated timeline
  and render the first10s from the actual rest initial condition. Never replace
  the requested rest-start with a later cropped, already-moving window.
- Set simu.rampTime=0 explicitly (the default is100), and verify the imported
  signal is not inadvertently ramped.
- Verify hydrostatic metadata through native netCDF/legacy files/BEMIO/H5.
  readCAPYTAINE can fall back to zero CG/CB/volume if missing; defaults cannot
  silently replace the canonical boat. Preserve full inertia products and frames.

These are workflow/input correctness checks, not intermediate coefficient
accuracy gates. No first-run result is claimed at this setup checkpoint.

Official references:

- https://wec-sim.github.io/WEC-Sim/main/user/getting_started.html
- https://wec-sim.github.io/WEC-Sim/main/user/workflow.html
- https://wec-sim.github.io/WEC-Sim/main/theory/theory.html
- https://github.com/WEC-Sim/WEC-Sim/releases/tag/v7.1.0
