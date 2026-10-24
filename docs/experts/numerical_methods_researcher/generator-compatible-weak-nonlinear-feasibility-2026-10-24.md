# Generator-compatible weakly nonlinear follow-up, 2026-10-24

**Superseded task:** the user subsequently clarified that the physical linear
regime is an a priori engineering assumption from wetted-surface/buoyancy curves
relative to WEC-Sim's mean-geometry assumptions. No weakly nonlinear trajectory,
native nonlinear extension or adapter is required. Preserve this source review
as historical/unexecuted feasibility; do not implement it. Current scope is ONE
buoy-generated canonical case targeting 5-10-degree roll and appreciable other
motions, supported by the geometry envelope reviewed in
[geometry-operating-envelope-review-2026-10-24.md](geometry-operating-envelope-review-2026-10-24.md).

Bounded read-only method proposal; no implementation/numerics. Baseline video
closure takes priority. Accepted artifacts remain immutable; Topwater solver
work stays held. Selected WEC-Sim source inspected at
`f962103a348d63101ce669c2611f7d7c109755c3`.

## Controlling input and target

Use the actual project sea-state generator to generate a higher-energy sea.
Preserve its spectrum, directional spread and actual seeded component realization
(frequencies, headings, amplitudes and phases). Orient the box along the generated
dominant direction using the project's established propagation convention.
Do not replace this with a beam sea, fitted period bank, rescaled old realization
or another random generator carrying the same numerical seed.

Desired roll is 5-10 degrees, subject to the user's retained physical-linearity
intent. Blind scaling of the 0.0191-degree baseline would require roughly
261-523 times its wave amplitudes; it is not a suitable default. Higher energy
must come from generator inputs. Existing-coefficient response predictions and
exact nonlinear righting/wetting diagnostics can reject unsuitable generated
cases; they cannot certify a nonlinear trajectory. No arbitrary percentage from
the earlier proposal is a new requirement.

## Verified native API boundary

Selected `waveClass.m` and model functions establish:

- `regular`/`regularCIC` call `waveElevReg`, which uses `cos(omega*t)` and does not
  consume arbitrary component phase. A regular-wave-object phase bank is not an
  existing faithful arbitrary-phase API.
- `spectrumImport` is type 21 and accepts frequency Hz, density m^2/Hz and an
  explicit third phase column. Source `setWavePhase`, `irregWaveSpectrum` and
  `waveElevIrreg` preserve those phases. For one heading, complex generator
  coefficient `c` under `Re(c*exp(i*k.s-i*omega*t))` maps to native phase
  `-arg(c)` and density `pi*abs(c)^2/dOmega`, using the ACTUAL native frequency
  widths. This gives native amplitude `sqrt(2*S_rad*dOmega)=abs(c)`.
- A type-21 file supplies one phase per frequency, not an arbitrary
  frequency-by-direction phase matrix. Multiple directions, multiple wave
  objects and independently seeded directional phases cannot be assumed
  equivalent without an exercised roundtrip.
- `spectrumImportFullDir` is type 35 and its setup currently generates random
  phases. Native `waveElevation` and `nonFKForce` implement spectral evaluation
  only for types below 30. Type 35 or imported elevation type 40 is therefore not
  a verified nonlinear replacement. The type-40 option-2 loophole remains invalid.

Primary documentation confirms the
[three-column spectrum import](https://wec-sim.github.io/WEC-Sim/main/user/code_structure.html#spectrumimport)
and distinguishes
[weakly nonlinear buoyancy/FK integration](https://wec-sim.github.io/WEC-Sim/main/user/advanced_features.html#nonlinear-buoyancy-and-froude-krylov-excitation).
The selected local source, rather than live documentation, determines these
version-specific limits.

## First choice and cheapest rejection

First obtain the generator's exact component-bank API from its owner. Attempt
the smallest native type-21 representation only if it can preserve that entire
bank. Before any boat run, compare independently evaluated generator eta and
incident pressure with native fields at several noncollinear surface/submerged
points and times, and compare all six linear excitation histories. This exposes
heading, phase, unit/energy and directional double-count errors. The comparison
uses an identity/roundoff or explicitly justified interpolation tolerance, not
an invented physical-linearity percentage. Reject the native mapping on any
unexplained mismatch; do not simplify the generator to make it fit.

For a general directional bank that the stock API cannot preserve, prefer a
narrow, independently reviewed shared-component-bank field/pressure adapter into
the existing native weakly nonlinear buoyancy/FK path while retaining its linear
diffraction/radiation dynamics. The same bank drives linear loads, nonlinear
incident pressure and rendering. This extends the consumer interface, not the
Topwater solver or a new fluid solver. It is a proposed extension, not an already
supported feature. A new incoming-startup construction must be explicitly
reconciled with the generator bank; no global temporal ramp or hidden phase/seed
change substitutes for that requirement.

Two preflight rejection checks precede even a short nonlinear trajectory:

1. Native `waveElevation.m` evaluates `t/rampT` when `t<=rampT`. At `t=0` and
   `rampTime=0`, source predicts NaN. Independently reproduce this finite-input
   endpoint before relying on a minimally reviewed zero-ramp correction. This
   is a source-derived lead, not an executed defect reproduction here.
2. Native nonlinear buoyancy/FK uses panel-centroid wet/dry classification. Refine
   the full nonlinear integration geometry separately from the frequency-domain
   BEM mesh; prove zero-wave equilibrium and stable righting/pressure integrals.
   The minimum display STL is not nonlinear-convergence evidence.

Then WEC-Sim owner runs the same generated case with linear and weakly nonlinear
hydrodynamics from actual rest, comparing all six full histories and relative
water-to-moving-hull wetting. At 5-10-degree roll, inspect actual righting,
chine/rim interactions and model limitations rather than assuming small angles.
Passing native weakly nonlinear checks addresses buoyancy/FK geometry only;
unmodeled nonlinear diffraction/radiation, viscosity, slamming and flooding do
not become certified. Coordinator and user retain acceptance decisions.

Stopping point: this source-backed feasibility/proposed rejection plan. No
higher-energy generated case or 5-10-degree response is claimed. Implementation
and serial numerical allocation belong to the respective owners after baseline
closure and the coordinator's bounded assignment.
