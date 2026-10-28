# `topwater solve`

Solve a boat/sea example and write its complete `result.npz` and six-DOF
motion plot `result.png` after motion and total free-surface validation.

```powershell
uv run topwater solve [DIRECTORY] DWELL [--rest SECONDS] [--ramp SECONDS] [--solver=native|wecsim]
uv run topwater solve examples/simple-boat 100 --rest 10 --ramp 5 --solver=wecsim
```

DIRECTORY defaults to the current directory. The input/output filenames are
fixed. Full-complex surface reconstruction is enabled following the selected
60-second simple-boat qualification: maximum elevation difference plus reference
uncertainty is 5.78 mm against the chosen linearized model, below the 10 mm target.
This is conditional selected-case evidence, not a guarantee for other inputs.
See [the qualification and its assumptions](../experts/numerical_methods_researcher/issue31-selected-case-acceptance-2026-10-27.md).
The new dwell-based 130-second packet integration is awaiting its own production
validation; the earlier 60-second trajectory is historical evidence.

The native backend is the intended default. Native coupled dynamics are not yet
available; omitting `--solver=wecsim` must report that limitation without silently
using another solver.

## Inputs

The example directory supplies `boat.npz` and `sea.npz`. The shared file-format
documentation defines their numeric arrays, units, and coordinate conventions.
DWELL is the positive full-amplitude target interval in seconds. Rest defaults
to 10 s on each end; ramp defaults to 5 s on each side. Total simulation time is
`2*rest + 2*ramp + dwell`: the example above runs 130 s, including both endpoints.
Each timing input must align to the fixed 0.01 s integration/packet grid. Ramp
must be positive and rest nonnegative. The end time is fixed; settling is
reported through the motion, not used to extend this run automatically.

Boat inputs include the complete physical convex hull shell, mass, body CG, and
full inertia tensor about CG. The single opening is inferred from the uncapped
shell; its artificial topology cap is an internal product and carries no pressure.
Sea inputs supply explicit component amplitudes, angular frequencies, global
propagation directions, and phases. Solving does not regenerate the sea or change
the input files. Their phase numbers refer directly to the midpoint of the
requested dwell, 65 s in the example, using
`cos(k dot x - omega*(t-midpoint) + phase)`. Phases are not retimed inside the
sea archive. The target at the spatial origin is a carrier multiplied by the
quintic rest/ramp/dwell/ramp/rest envelope.

## Outputs and replacement

`result.npz` contains computed boat motion and total water-air surface geometry.
It is a visualization/comparison result, not a restart file.

`result.png` contains the complete six motion histories: surge, sway and heave
in metres, with roll, pitch and yaw perturbations in degrees. It is derived from
all saved motion samples, including startup; plotting never recomputes physics.
The axes and angle conventions are described below.

Both files are prepared privately before publication. Individual replacements
are atomic; handled failures restore the previous pair. If the operating system
also blocks rollback, the error identifies retained recovery copies. Abrupt
process termination between replacements cannot guarantee pair consistency;
an invocation without a successful log END is incomplete.

`result.log` records the latest invocation, stages, backend output, warnings,
fatal errors, timings and publication status. It is replaced when solving starts,
including failed attempts. See [the common log format](logging.md) and
[the current warning triage](warning-triage.md).

The solve does not render or remove `video.mp4`. An existing video may therefore
belong to an earlier result; tracking that association is the user's responsibility.

Capytaine exports, BEMIO H5/MAT files, STL, Simulink models/caches, and native
audit arrays live in operating-system temporary storage and are removed after
execution. No public persistent cache is produced. Diagnostics go to the terminal.

A failed surface export may retain private field/replay data in its system-temp
directory so export can be diagnosed without repeating the body/field solve.
The private cache records its computational mean-water frame and vertical world
offset. Public replay restores that offset once to both motion and surface;
older caches default to zero offset. No extra public NPZ metadata is required.
The original stage error is preserved and the recovery location is printed.
If temporary-file cleanup also fails, that secondary failure is reported without
replacing the original error. A cleanup failure after otherwise successful work
is fatal and prevents publication.

The motion plot uses a fixed inertial frame: the initial bow projected
horizontally, port, and world up. Translations are CG displacement from the
initial position in those axes. Rotations are initial-attitude departures,
expressed in that same frame as unwrapped extrinsic XYZ Euler angles; initial
trim is removed. All native motion timestamps are used, with no smoothing.
The PNG is 2160 by 1170 pixels (12 by 6.5 inches at 180 dpi).
Translations share one symmetric Y range in metres; rotations share a separate
symmetric Y range in degrees. Each range includes all three histories in its row.
Each row has exactly seven equally spaced major Y ticks at nice decimal values,
with zero always in the middle. The first and last ticks are the axis limits;
the common row range retains at least 5% headroom over all saved motions.

An initial bow with no resolved horizontal projection is a fatal plot limitation
before either new output is published. Euler gimbal lock emits a presentation
warning; the saved quaternions remain authoritative. Plot encoding/filesystem
failures are fatal, with the previous output pair preserved when rollback works.

## WEC-Sim dependencies

The reviewed path uses Capytaine 3.0.0, MATLAB R2025b with matlabengine 25.2.2,
and licensed Simulink, Simscape, and Simscape Multibody. These are optional
dependencies of the WEC-Sim workflow, not requirements for rendering.

The maintained fork is pinned to
`6b774228e1f6b31dee473e626204b0ef65e5ec38`. Set `TOPWATER_WECSIM_PATH` to its
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
thread limit, figure default, the newer-model load policy (`ErrorIfLoadNewModel`),
and Simulink cache/code-generation settings. It
restores those settings on success or failure and closes only the unique model
owned by this run. Restoration must finish before scratch deletion. If restoration
cannot be verified, the solve fails and retains its system-temporary directory
with a recovery location in stderr; the previous public result remains intact.
If native work and engine restoration both fail, the error retains
`EngineCleanupError` classification and recovery scratch, while its displayed
message includes the primary MATLAB failure/transcript and cleanup failure.

This lifecycle passed repeated normal/forced-failure checks in the existing
hidden shared engine, including removal of each temporary class directory before
the next cycle. [Verification and limits](../experts/wecsim_expert/workflow-engine-lifecycle-2026-10-25.md).
That pass does not qualify the body/surface physics pipeline.
The two additional current-path lifecycle regression sources are documented in
[source checkpoint](../experts/wecsim_expert/current-path-lifecycle-corrections-2026-10-25.md);
their focused regression checks passed within the linked checkpoint's scope.

## Assumptions and hard-coded policies

The standard public-interface spectral handoff passed the selected canonical60s
force, startup, applied-force, rest/radiation-memory and motion comparisons at the
pinned6b774228 reference. Fresh total-surface reconstruction, six-DOF PNG and all
751 frames through the ordinary render CLI also passed independent review.
These checks establish compatibility for that case; the physical limitations below
remain unqualified. [Evidence and scope](../experts/wecsim_expert/current-state.md).

- Open water, infinite depth, zero forward speed, density 1025 kg/m3 and gravity
  9.81 m/s2.
- Boat horizontal CG starts at the origin. Bow faces the energy-weighted dominant
  incoming direction; an energetic sea with no resolved dominant direction is
  rejected. Calm water uses zero heading.
- Calm-water equilibrium balances heave, roll and pitch before simulation time
  zero. The heave/roll/pitch branch must be stably restoring. Initial physical
  pose is that equilibrium, velocity and radiation memory are zero.
- Each original carrier gets its own real quintic time packet at the origin.
  One FFT supplies DC and all contiguous positive harmonics up to the first
  certified passing cutoff, at frequency spacing `1/total_time`. Negative
  harmonics are the real signal's conjugates, not additional opposite-heading
  waves. Every retained positive harmonic propagates in its original carrier
  heading with its own deep-water dispersion. There is no force multiplier.
- The sole packet acceptance criterion is that absolute elevation never exceeds
  its target amplitude envelope by more than its allocated tolerance. Allocate
  the aggregate 0.01 m budget proportional to original carrier amplitude, so the
  sum of conservative per-carrier bounds stays within 0.01 m at the origin.
  Fixed 10 ms samples plus a curvature/roundoff allowance certify that bound
  between samples. Dwell waveform mismatch is a diagnostic, not another gate.
  Calm conditions over the actual full hull require a separate startup check.
- Original carrier eligibility remains
  `f <= min(0.6 Hz, sqrt(g/(4*pi*beam)))`. Generated packet sidebands above this
  operating assumption are retained and warned about, with cutoff, count and
  amplitude L1 contribution. All require genuine coefficient coverage. This
  rule is not a mesh/time convergence certificate; finite packets differ from
  the original stationary sea spectrum.
- A nonzero packet DC is preserved as uniform static mean water level `d`.
  BEM geometry, CG and restoring stay in the original computational mean-zero
  frame. Native initial heave is `-d`; native restoring therefore supplies the
  full static `d*K[:,heave]` load, including coupled moments. Exported body and
  water positions receive `+d` once. No zero-frequency diffraction/radiation
  solve, artificial epsilon-frequency wave or applied force history is added.
- WEC-Sim computes harmonic excitation from one full-directional spectrum with
  supplied phases, then integrates body motion. Topwater encodes the unchanged
  coherent physical components as spectrum density/spreading/phase inputs; no
  applied force history is supplied. Signed-kappa branches retain their separate
  headings, phases, absolute time and spatial origin, with no additional ramp.
  Equal frequency/heading coefficients are summed coherently; exactly zero
  direction rows are omitted. A single-heading sea has an exact-zero second
  column solely to satisfy native FullDir's directional-grid requirement.
- Topwater samples coefficient inputs at every actual component frequency and
  heading using the retained complex bilinear interpolation. Below 0.1 rad/s,
  `H = H0 + (omega/0.1)^2 * (H(0.1,heading)-H0)`, where `H0 = K[:,heave]`.
  This is a declared coefficient approximation, not a new WEC-Sim capability or
  a substitute for spatial free-surface fields. Ordinary native spectral
  preparation consumes these coefficient nodes and computes the force sum.
- Official BEMIO import/normalization and `rebuildHydroStruct` prepare the existing
  public `bodyClass(hydroStruct)` input. Its ordinary excitation axes contain
  derived nodes. The consumed radiation IRF/time and direct Ainf remain from
  a genuine 0.1--12 rad/s operator, extended if packet coverage requires it.
  Unused raw A/B and period metadata
  are omitted from this fixed FullDir+CIC struct; no raw coefficient is invented.
  The genuine original producer arrays remain unchanged. Spectral forcing needs
  no excitation IRF or H5 export; the legacy imported-elevation study path retains
  ordinary excitation-IRF/H5 preparation under its original policies.
- Native WEC-Sim dynamics retain fixed-step ode4, 0.01 s body/CIC sampling,
  10 s body-radiation memory, full mass/inertia coupling, and all six free modes.
  Viscous/linear damping, PTO, mooring, and prescribed motion remain off.
- Equilibrium orientation is included once in the mean geometry and inertia.
  Global heading is an embedding applied once to output geometry and body pose.
- Surface reconstruction adds no feedback to body dynamics. It combines the
  represented incoming packet, direct sideband diffraction, and radiation driven
  by all six saved fixed-mean native modal velocities. Angular velocities are
  not replaced by Euler-angle rates.
- The radiation method preserves the full complex frequency response.
  Its default band is 0--48 rad/s at 0.025 rad/s spacing, with a 1024 s FFT
  period and FFT timestep half the saved motion timestep. Every saved velocity
  knot is retained. Velocity is zero before startup; a prescribed 20 s cubic
  Hermite continuation takes final velocity/acceleration smoothly to rest after
  the saved interval. This continuation is numerical padding, not predicted
  motion. Higher frequencies are set to zero, leaving tail uncertainty to be
  evaluated. Raw initial radiation is preserved; no separate acceleration
  response is added. Spatial chunks contain at most 128 evaluation points.
  These settings passed whole-domain validation for the selected simple-boat
  case against the 1 cm linear-model elevation requirement, including reference
  uncertainty. Other inputs/settings need their own accuracy assessment. Offline synthesis
  does not establish exact causal rest-start physics for arbitrary inputs.
- Disturbance fields near/inside the mean hull footprint use an explicit exterior
  trace at offset `min(0.025 m, 0.01*beam)`. The incident contribution is evaluated
  at the actual target locations. This continuation and geometric intersection
  do not solve nonlinear moving-contact physics.
- The water output uses a piecewise-linear triangular graph sampled at 0.1875 m
  grid spacing. Export retains the sampled vertices and original connectivity;
  it does not cut triangles or remesh at the posed hull. Surface portions inside
  that hull are invalid fluid. A separate shared intersection routine computes
  the waterline from the total surface and physical hull faces. Its finite patch covers the
  trajectory/default render view; its edge is an export boundary, not a wall.
  Surface timestamps target 0.08 s spacing, uniformly adjusted to include both
  endpoints. CG interpolation is linear and orientation interpolation is Slerp.

## Fatal errors and detection stage

| Stage | Detected fatal conditions |
|---|---|
| Readiness, before physical preparation | Missing/wrong dependency versions, missing maintained fork, unreviewed revision or modified fork source, unavailable total-surface capability. |
| Physical setup | Invalid inputs, unsupported incident support, unresolved default heading, missing force/moment balance, unstable supported equilibrium, initial rim immersion or lack of buoyant contact. |
| Incoming-field assessment, before body integration | Full-hull initial elevation/fluid-speed bound fails or remains inconclusive; incoming-wave envelope does not arrive during the requested interval. |
| Native preprocessing | Prepared physical CG, CB, volume, stiffness or infinity added mass differ; phase convention, finite coefficient/kernel data, full-radiation preservation or parameter-level full inertia identity fails; native initial force/moment exceeds calm limits; native FullDir setup loses physical components or coefficient-node identity. |
| Native integration/output audit | Solver failure, nonfinite/incomplete output, wrong grid/endpoints/initial state, nonzero initial radiation memory, applied excitation differs from WEC-Sim's prepared total, inconsistent force-component sum, Newton-Euler/CIC/added-mass-delay audit failure. |
| Required surface/export | Missing or invalid total surface, unsupported contact/export geometry, initial represented total elevation exceeds the local calm limit, incomplete coverage, or invalid common result arrays. |

Calm-start engineering limits are physical `max|eta| <= 0.01 m` and
`max|u| <= .001*sum(a*omega)` throughout the initial full-hull bounding region.
Initial net excitation-plus-restoring force/moment norms must not exceed `.001*m*g` and `.001*m*g*L`, with
`L` the maximum hull dimension, and every native component must also be within
`.001` of its carrier-load scale. The current combined incoming-field bound uses order-six
Taylor boxes, initial widths 1 m/1 m/0.2 m, at most three refinement levels, and a
480 s assessment budget. An inconclusive bound is reported distinctly from an
observed violation. At least 25% combined incoming envelope must arrive.

## Warnings

Warnings report physical applicability, numerical limitations or operational
conditions while supported computation continues. They do not become extra NPZ
fields. Restoring diagnostics print the
first event per category and summarize its sample count, first/last time, and
maximum departure. Wet-set departures are detected at each saved surface frame;
their first/last exceedance and worst fraction/time are emitted when surface
export completes.

The implemented restoring diagnostic compares exact **mean-sea-reference calm-water** restoring at
each computed pose with the linear fixed-mean `-K*q` wrench. Force and moment are
assessed separately. A warning starts above 10% relative difference, normalized
by the larger exact/linear norm with floors `.001*m*g` and `.001*m*g*L` to avoid
percentages of numerical zero. This is an engineering reporting convention,
not a nonlinear dynamic error bound. If flat-reference geometry reaches the rim
or loses contact, the restoring diagnostic reports that it is unavailable there.

The surface producer warns when the symmetric-difference area between the
represented total-water wetted set and the fixed mean wetted set exceeds 10% of
the initial hydrostatic wetted area. It reports first/last exceedance and the
worst fraction/time. This uses the reconstructed piecewise-linear total graph;
it is not a nonlinear-physics truth estimate.

The method logs its band, FFT period/timestep, input interpolation,
future continuation, pre-response and raw startup value in `SURFACE_NUMERICS`
records. These are numerical diagnostics, not automatic accuracy certificates.
The earlier cosine projection and failed causal-only cut remain investigation
evidence; their old discrepancies are not reused as new-run error estimates.
The producer also reports direct sideband
field-versus-prepared-body load differences, without modifying the accepted
body loads. Total-minus-incident surface metrics describe the produced geometry;
they are not an accuracy certificate.

Delta moments exclude hull interiors. They integrate the represented piecewise-
linear field over valid exterior support using temporary polygons for numerical
integration, without changing the exported mesh. The one-beam summary selects
source triangles by their centroid's horizontal distance to the separately
computed total-surface waterline, then integrates their exterior portions. This
is an explicit discrete approximation to the band, not an exact geometric
distance contour. Full-patch and band summaries report area/time weighting,
RMS, incident-relative RMS when defined, and the resolved peak/time/location.

Capytaine mesh-resolution and irregular-frequency warnings remain visible. Their
first occurrence is retained and repeated families are summarized; they are not
silently suppressed or converted into numerical-accuracy claims.

## Current support and detection limits

The former cosine projection is historical. The full-complex method passed
issue #31's selected-case qualification at all 5929 targets and 751 saved frames
through 60 s, including exterior-trace continuation. Maximum runtime/reference
difference was 0.341 mm; the largest pointwise difference plus reference
uncertainty was 5.784 mm. The 1 cm requirement concerns numerical elevation error
against the chosen linearized model, not nonlinear calculations or real water.
The omitted-frequency allowance assumes that each regular response stays below
twice its maximum calibrated from 36--48 and 64/96 rad/s, with no unseen larger
resonance. This is an explicit engineering assumption, not a proven universal
bound or a qualification of the underlying BEM mesh. Runtime input,
initial-state, topology, finite-data, and output gates remain mandatory. A
motion-only result is never substituted for the required total-surface output.

The full-hull initial incident elevation/velocity bounds remain in force. The
producer additionally checks initial **represented total elevation** on local
water triangles touching the boat's initial horizontal bounding region. No
full-hull total disturbed-fluid velocity bound has been established.

Linear WEC-Sim hydrodynamics act about a fixed mean wetted hull. Reconstructing
additional surface output does not add feedback or alter its governing dynamics.
Geometric intersection with a moving hull does not establish nonlinear moving
contact, flooding, slamming, breaking, or viscosity. Those events require explicit
supported models/detectors before the command can claim to assess them.

Finite results and passing internal balances do not establish mesh/time/domain
convergence, accuracy, or passivity. Radiation can temporarily return stored
energy; its instantaneous work sign is not an automatic fatal-error criterion.
