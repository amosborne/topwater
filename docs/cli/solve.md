# `topwater solve`

Solve a boat/sea example and write its complete `result.npz` and combined
motion/sea-spectrum plot `result.png` after motion and total free-surface validation.

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
The dwell-based 130-second packet integration completed its native audits. Its
saved motion was audited and continued after an engine-cleanup failure; this
does not claim a successful original solve invocation. The current simplified
surface path retains the uncut grid and basic data validation; moving-waterline,
wetted-area, disturbance-statistics and opening-clearance work is absent from
normal solve. The earlier 60-second surface-error qualification does not qualify
the changed packet. See the [historical integration record](../experts/packet-integration-run-2026-10-27.json).

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

`result.npz` contains computed boat motion, the full display hull, and total
water-air surface geometry. It stores no sea-spectrum arrays or plot metadata.
The archive remains sufficient for standalone video rendering and is a
visualization/comparison result, not a restart file.

`result.png` contains the complete six motion histories and the original sea's
polar spectrum. Surge, sway and heave are in metres; roll, pitch and yaw
perturbations are in degrees. All saved motion samples, including startup, are
used. The original `sea.npz` data already loaded by solve are passed directly to
plot generation; neither the plot nor video rendering recomputes physics.

The spectrum uses 15-degree by 0.025 Hz bins, colored by RMS elevation:
`sqrt(0.5 * sum(amplitude**2))` for the original carriers in each bin. Radius is
frequency; the colorbar is labeled **Bin RMS elevation (m)**. Coming-from directions add pi to
the original propagation angles. The polar chart is rotated so the default
boat heading points up: the default is the incoming direction of the
amplitude-squared-weighted mean propagation vector. A dotted line marks the
actual initial boat heading relative to that default. In the current workflow
they coincide. Angular spokes and compass ticks are omitted; frequency labels
point down and the dotted heading line remains. The enlarged polar panel fills
the right column and has a horizontal colorbar below it. The upper color bound
is rounded outward to two significant digits, with nice interior ticks and exact
endpoints (0 to 0.12 m for this case); this never changes or clips bin values. Neither the spectrum nor the six-motion figure has
an overall title/subtitle or footer; motion x-axis labels are **Time (s)**.
Significant wave height is `4*sqrt(sum(amplitude**2)/2)`, but is not drawn as a heading.
The plot uses the perceptually uniform cividis amplitude scale, with zero bins
colored at its lower end; empty bins are not masked white. Labels remain black.
These are original sea carriers, not generated packet sidebands.
No spectrum is added to the video or the result archive.

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
The combined PNG is 3150 by 1170 pixels (17.5 by 6.5 inches at 180 dpi). Direct
plot API calls without a sea retain the original 2160 by 1170 motion-only PNG.
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

If automation overrides TMP/TEMP, configure both initial Simulink cache and
code-generation folders as siblings of that temporary folder before invoking
the adapter. The temporary folder must not sit beneath either initial folder.
Otherwise MATLAB can reject restoration after otherwise completed physics.
The adapter preserves the original session settings rather than replacing them.

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
- Radiation preserves the full complex frequency response. The default uses
  exactly the saved native 10 ms velocity samples: one FFT and inverse FFT,
  with no oversampling, zero-padding or invented future motion. Odd sample
  counts (including 13001 samples over 130 s) are retained; the DFT period is
  sample count times timestep (130.01 s in that example). The transfer bank
  remains 0--48 rad/s at 0.025 rad/s spacing, with linear real/imaginary
  interpolation and zero response above that band. Raw initial radiation is
  preserved and no separate acceleration response is added. Spatial chunks
  contain at most 128 evaluation points. Explicit private settings can retain
  the historical padded policy; its earlier qualification does not qualify
  these new defaults. Saved-field calm analysis is separate from solve and
  checks elevation against +/- 0.01 m without changing the reconstructed field.
- Interior target locations use diffraction/radiation evaluated at the closest
  point on the mean-waterline boundary, with zero outward offset. Exterior
  targets retain their actual locations. Only explicit boundary-trace queries
  are admitted; strictly interior Green-function queries remain forbidden.
  Incident water is always evaluated at the actual target. The small boundary
  evaluation screen found finite values and exterior approaches toward them;
  its original 0.1% comparison screen failed and is retained as evidence, not
  relabeled a pass. Finite offsets include real spatial variation, particularly
  at high frequency, so that failure does not establish a bad boundary value.
  This engineering choice does not certify a 1 cm elevation-error bound or
  solve nonlinear moving-contact physics. Existing saved fields retain their
  original offset policies until explicitly regenerated.
- The water output uses a piecewise-linear triangular graph sampled at 0.1875 m
  grid spacing. Export retains the sampled vertices and original connectivity;
  it does not cut triangles or remesh at the posed hull. Surface portions inside
  that hull are invalid fluid. Rendering alone computes the moving waterline
  from the saved total surface and physical hull faces. Its finite patch covers the
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
| Required surface/export | Missing/nonfinite total surface, malformed clocks/connectivity, incomplete coverage, or invalid common result arrays. Opening overtopping is retained as data; surface export does not apply a moving-opening clearance gate. |

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
maximum departure. Normal surface export does not compute moving waterlines,
wetted-area comparisons, disturbance statistics or opening-clearance checks.

The implemented restoring diagnostic compares exact **mean-sea-reference calm-water** restoring at
each computed pose with the linear fixed-mean `-K*q` wrench. Force and moment are
assessed separately. A warning starts above 10% relative difference, normalized
by the larger exact/linear norm with floors `.001*m*g` and `.001*m*g*L` to avoid
percentages of numerical zero. This is an engineering reporting convention,
not a nonlinear dynamic error bound. If flat-reference geometry reaches the rim
or loses contact, the restoring diagnostic reports that it is unavailable there.

`SURFACE_NUMERICS` records describe the frequency band, native FFT length and
timestep, transfer interpolation, unchanged initial sample and lack of padding.
These are numerical diagnostics, not automatic accuracy certificates. Calm-water
analysis of saved fields is a separate diagnostic. A height above 1 cm does not
by itself identify temporal wraparound rather than a physical wave or body wake.

Historical clipping-based wet-set and total-minus-incident statistics remain
available only through explicit private cache analysis (`analyze=True`). They
are excluded from normal solve and default cache replay. This historical path
retains its original geometry checks and does not modify cached fields.

Capytaine mesh-resolution and irregular-frequency warnings remain visible. Their
first occurrence is retained and repeated families are summarized; they are not
silently suppressed or converted into numerical-accuracy claims.

## Current support and detection limits

The former cosine projection is historical. The historical padded full-complex policy passed
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

The earlier padded-policy evidence above does not qualify the native unpadded
FFT or exact mean-waterline trace defaults. The full-hull initial incident
elevation/velocity bounds remain in force. A separate saved-field calm analysis
uses the +/- 1 cm elevation requirement; normal surface packing performs no clipped
startup-elevation or disturbed-fluid velocity bound.

Linear WEC-Sim hydrodynamics act about a fixed mean wetted hull. Reconstructing
additional surface output does not add feedback or alter its governing dynamics.
Geometric intersection with a moving hull does not establish nonlinear moving
contact, flooding, slamming, breaking, or viscosity. Those events require explicit
supported models/detectors before the command can claim to assess them.

Finite results and passing internal balances do not establish mesh/time/domain
convergence, accuracy, or passivity. Radiation can temporarily return stored
energy; its instantaneous work sign is not an automatic fatal-error criterion.
