# `topwater solve`

Solve a boat/sea example and replace its complete `result.npz` only after the
motion and total free-surface output have passed validation.

```powershell
uv run topwater solve [DIRECTORY] DURATION_SECONDS [--solver=native|wecsim]
uv run topwater solve examples/simple-boat 60 --solver=wecsim
```

DIRECTORY defaults to the current directory. The input/output filenames are
fixed. This checkpoint implements the interfaces but keeps total-surface
publication closed pending the documented motion-driven/contact pilot.

The native backend is the intended default. Native coupled dynamics are not yet
available; omitting `--solver=wecsim` must report that limitation without silently
using another solver.

## Inputs

The example directory supplies `boat.npz` and `sea.npz`. The shared file-format
documentation defines their numeric arrays, units, and coordinate conventions.
The requested duration is a positive number of seconds. WEC-Sim currently
requires a duration that is an integer multiple of its fixed 0.02 s time step.

Boat inputs include the complete physical convex hull shell, mass, body CG, and
full inertia tensor about CG. The single opening is inferred from the uncapped
shell; its artificial topology cap is an internal product and carries no pressure.
Sea inputs supply explicit component amplitudes, angular frequencies, global
propagation directions, and phases. Solving does not regenerate the sea or change
the input files.

## Outputs and replacement

`result.npz` contains computed boat motion and total water-air surface geometry.
It is a visualization/comparison result, not a restart file. The caller publishes
the complete validated file atomically. A failure preserves an existing result.
The adapter itself never writes a partial public result.

The solve does not render or remove `video.mp4`. An existing video may therefore
belong to an earlier result; tracking that association is the user's responsibility.

Capytaine exports, BEMIO H5/MAT files, STL, Simulink models/caches, and native
audit arrays live in operating-system temporary storage and are removed after
execution. No public persistent cache is produced. Diagnostics go to the terminal.

## WEC-Sim dependencies

The reviewed path uses Capytaine 3.0.0, MATLAB R2025b with matlabengine 25.2.2,
and licensed Simulink, Simscape, and Simscape Multibody. These are optional
dependencies of the WEC-Sim workflow, not requirements for rendering.

The maintained fork is pinned to
`f962103a348d63101ce669c2611f7d7c109755c3`. Set `TOPWATER_WECSIM_PATH` to its
checkout directory, or use the repository-local fallback
`artifacts/dependencies/topwater-wec-sim`. A different or unverifiable revision
fails readiness rather than running unreviewed engine code.

The adapter attaches to an existing shared MATLAB engine named
`CAPYTAINE_REFERENCE`. It does not start or quit MATLAB. Engine visibility and
required licenses are checked before the native body run. Project automation
must follow the serial numerical execution procedure before calling the workflow.

MATLAB stays minimized/headless: the adapter never opens the desktop, model
editor, or Mechanics Explorer and temporarily disables new figure visibility.
Before any model construction or initialization, a guard snapshots the original
working directory, search path, diary filename/on-off state, warning settings,
thread limit, figure default, and Simulink cache/code-generation settings. It
restores those settings on success or failure and closes only the unique model
owned by this run. Restoration must finish before scratch deletion. If restoration
cannot be verified, the solve fails and retains its system-temporary directory
with a recovery location in stderr; the previous public result remains intact.

This lifecycle passed repeated normal/forced-failure checks in the existing
hidden shared engine, including removal of each temporary class directory before
the next cycle. [Verification and limits](../experts/wecsim_expert/workflow-engine-lifecycle-2026-10-25.md).
That pass does not qualify the body/surface physics pipeline.

## Assumptions and hard-coded policies

- Open water, infinite depth, zero forward speed, density 1025 kg/m3 and gravity
  9.81 m/s2.
- Boat horizontal CG starts at the origin. Bow faces the energy-weighted dominant
  incoming direction; an energetic sea with no resolved dominant direction is
  rejected. Calm water uses zero heading.
- Calm-water equilibrium balances heave, roll and pitch before simulation time
  zero. The heave/roll/pitch branch must be stably restoring. Initial motion and
  radiation memory are zero; no simulated warmup is hidden in the duration.
- The incoming train uses the accepted fixed initial spatial window: group-time
  center 40 s, plateau half-width 15 s, support half-width 25 s, and quintic
  transitions. Signed-kappa dispersion propagates that field. Requested duration
  changes the stopping time, not the startup field. There is no temporal force ramp.
- Sideband angular-frequency spacing is 0.005 rad/s. The project support limit
  is `f <= min(0.6 Hz, sqrt(g/(4*pi*beam)))`. This is an operating-envelope
  assumption, not a mesh or time-convergence certificate. A finite packet's
  spectrum differs from the stationary measured spectrum.
- Body excitation preserves the existing periodic heading interpolation and
  declared low-frequency bridge: below 0.1 rad/s,
  `H = H0 + (omega/0.1)^2 * (H(0.1,heading)-H0)`, with the exact constant-pressure
  limit `H0 = K[:,heave]`. This approximation applies to body loads; it is not
  substituted for missing spatial surface fields.
- Native WEC-Sim dynamics retain fixed-step ode4, 0.02 s body/CIC sampling,
  10 s body-radiation memory, full mass/inertia coupling, and all six free modes.
  Viscous/linear damping, PTO, mooring, and prescribed motion remain off.
- Equilibrium orientation is included once in the mean geometry and inertia.
  Global heading is an embedding applied once to output geometry and body pose.

## Fatal errors and detection stage

| Stage | Detected fatal conditions |
|---|---|
| Readiness, before physical preparation | Missing/wrong dependency versions, missing maintained fork, unreviewed fork revision, unavailable total-surface capability. |
| Physical setup | Invalid inputs, unsupported incident support, unresolved default heading, missing force/moment balance, unstable supported equilibrium, initial rim immersion or lack of buoyant contact. |
| Incoming-field assessment, before body integration | Full-hull initial elevation/fluid-speed bound fails or remains inconclusive; initial force/moment bounds fail; incoming-wave envelope does not arrive during the requested interval. |
| Native preprocessing | H5/physical CG, CB, volume, stiffness or infinity added mass differ; phase-convention mapping or compiled full inertia fails; imported force does not bind to exactly one positive total-force branch. |
| Native integration/output audit | Solver failure, nonfinite/incomplete output, wrong grid/endpoints/initial state, nonzero initial radiation memory, applied excitation differs from prepared input, inconsistent force-component sum, Newton-Euler/CIC/added-mass-delay audit failure. |
| Required surface/export | Missing or invalid total surface, unsupported contact/export geometry, incomplete coverage, or invalid common result arrays. Surface publication is unavailable until the physical reconstruction pilot has passed. |

Calm-start engineering limits are `max|eta| <= .001*sum(a)` and
`max|u| <= .001*sum(a*omega)` throughout the initial full-hull bounding region.
Initial force/moment norms must not exceed `.001*m*g` and `.001*m*g*L`, with
`L` the maximum hull dimension, and every native component must also be within
`.001` of its carrier-load scale. The current combined-field bound uses order-six
Taylor boxes, initial widths 1 m/1 m/0.2 m, at most three refinement levels, and a
480 s assessment budget. An inconclusive bound is reported distinctly from an
observed violation. At least 25% combined incoming envelope must arrive.

## Warnings

Warnings indicate reduced physical applicability while supported computation
continues. They do not become extra NPZ fields. Print the first event per category
and summarize its sample count, first/last time, and maximum departure.

The implemented restoring diagnostic compares exact **calm-water** restoring at
each computed pose with the linear fixed-mean `-K*q` wrench. Force and moment are
assessed separately. A warning starts above 10% relative difference, normalized
by the larger exact/linear norm with floors `.001*m*g` and `.001*m*g*L` to avoid
percentages of numerical zero. This is an engineering reporting convention,
not a nonlinear dynamic error bound. If flat-reference geometry reaches the rim
or loses contact, the restoring diagnostic reports that it is unavailable there.

The proposed 10% actual-versus-assumed wetted-set warning requires the accepted
total-surface/contact reconstruction. It is **not currently implemented**. The
workflow must not call an incident-only or flat-water proxy the actual wetted set.

## Current support and detection limits

The generic source is being integrated. The total-surface readiness gate is
intentionally closed until endpoint/tail treatment, causal radiation transfer,
and near-hull export have passed independent verification. Presence of field
evaluation functions alone does not open that gate. No successful motion-only
result is substituted for the mandatory total-surface output.

Linear WEC-Sim hydrodynamics act about a fixed mean wetted hull. Reconstructing
additional surface output does not add feedback or alter its governing dynamics.
Geometric trimming against a moving hull does not establish nonlinear moving
contact, flooding, slamming, breaking, or viscosity. Those events require explicit
supported models/detectors before the command can claim to assess them.

Finite results and passing internal balances do not establish mesh/time/domain
convergence, accuracy, or passivity. Radiation can temporarily return stored
energy; its instantaneous work sign is not an automatic fatal-error criterion.
