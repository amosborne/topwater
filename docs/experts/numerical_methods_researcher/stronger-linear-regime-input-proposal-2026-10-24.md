# Stronger useful reference input: bounded proposal, 2026-10-24

**Input-selection prescription superseded by later user steering:** use the
project sea-state generator at higher energy, preserve its spectrum, directional
spread and seeded realization, and point the box along its generated dominant
direction. Do not scale the old sea or prescribe a bespoke beam-sea/period bank.
The later hard target is appreciable roll of 5-10 degrees while retaining the
intent to test the physical linear regime. Current bounded backend feasibility
is in [generator-compatible-weak-nonlinear-feasibility-2026-10-24.md](generator-compatible-weak-nonlinear-feasibility-2026-10-24.md).
The amplitude ladder and numerical thresholds below are an earlier unexecuted
proposal, not current instructions or requirements. The user subsequently also
cancelled the nonlinear-dynamics/adapter prescription: current support is an
a priori geometry envelope for ONE generated canonical case, reviewed in
[geometry-operating-envelope-review-2026-10-24.md](geometry-operating-envelope-review-2026-10-24.md).

Read-only follow-up proposal requested while the minimum-case video finishes.
No numerical jobs, new solver, new sea qualification or Topwater unhold. Accepted
packet/data remain immutable. Selected WEC-Sim fork/source checked at
`f962103a348d63101ce669c2611f7d7c109755c3`; Capytaine 3.0.0.

## Decision and useful scale

Increase the response until physical nonlinear corrections become material,
then retain the strongest passing input. Start with a controlled amplitude ladder
`alpha=[1,2,4,8]`, preserving frequencies, directions, phases and envelope shape
for the candidate prediction. This isolates amplitude effects before considering
a different sea spectrum. Do not require every axis to exceed a historical
absolute tolerance by increasing waves regardless of physical consequences.

Native baseline peaks from `case-attempt1/summary.json`:

| Amplitude factor | Translation peaks, mm | Rotation peaks, deg |
|---|---|---|
| 1, executed | 5.222 / 4.970 / 5.935 | 0.0191 / 0.0306 / 0.0204 |
| 4, linear forecast only | 20.889 / 19.880 / 23.742 | 0.0765 / 0.1224 / 0.0817 |
| 8, linear forecast only | 41.778 / 39.759 / 47.484 | 0.1529 / 0.2447 / 0.1634 |

The earlier 1 cm / 0.1 degree comparison budgets exceed current response peaks;
absolute-only checks could therefore accept nearly zero output. About 5.23 times
the baseline is needed merely to exceed all those historical budgets, with roll
limiting. This is motivation for stronger input, not a new acceptance requirement
or evidence that factor eight is physically linear.

Compare all six full histories using `q_alpha/alpha` against `q_1`, with per-axis
peak/RMS normalizations and weak axes explicitly reported. Preserve physical
rotation conventions, using relative rotation matrices/logarithms where needed.
Exact scaling in a linear hydrodynamic model verifies implementation consistency;
it cannot establish the nonlinear boat's linear regime. Nonlinear rigid-body
kinematics alone also cannot test changing wetted hydrodynamics.

## Cheapest physical rejection experiment

Capytaine owner first evaluates exact geometric buoyancy/restoring on the full
physical hull at the candidate scaled poses, with no new frequency-domain BEM
solves. Transform the full hull, clip at mean water level, retain/reapply reviewed
degree-two GL2 integration, and compute volume/CB and gravity-plus-buoyancy wrench
about instantaneous CG. Compare with `-K*q_alpha`. Start at at most 24 frozen
times selected by extrema of heave, roll, pitch and relative water-to-hull
excursion; inspect all 501 samples only for a promising candidate. This is a
genuine nonlinear geometry check, independent of linear self-scaling, and cheap
enough to reject an excessive amplitude before another dynamics run.

Also inspect the incident surface relative to the MOVING hull:
`eta(X_hull(t),t)-z_hull(t)`, rather than wave amplitude or CG heave alone.
Report waterplane area/second-moment changes, wetting at chines/knuckles, rim
clearance, and wavelength/steepness. Ordinary waterline motion is expected; a
chine crossing is a diagnostic, not an automatic failure. Deck flooding,
slamming/breaking or loss of contact reject the supported regime.

Proposed screening rule for coordinator approval BEFORE execution: reject a
candidate if nonlinear restoring differs by more than 5 percent of the baseline
linear restoring amplitude scaled by alpha in an informative heave/roll/pitch
axis, or a prohibited geometric event occurs. Report near-zero restoring axes
without denominator floors. Passing only means 'not rejected by this screen';
flat-water restoring omits incident-pressure, diffraction and radiation changes.

## Required nonlinear confirmation and feasible routes

The selected WEC-Sim source is decisive:

- `initializeWecSim.m:330-335` rejects `elevationImport` with nonlinearHydro 1.
- `waveElevation.m` and `nonFKForce.m` implement their spectral surface/pressure
  branches only for `typeNum<30`; imported type 40 falls through. Setting option
  2 is not a supported loophole.
- `nonLinearBuoyancy.m` and `nonFKForce.m` integrate at triangle centroids with
  wet/dry classification. The minimum display STL cannot independently establish
  nonlinear convergence. Nonlinear integration needs a suitably refined full
  geometry, separate from the unchanged frequency-domain solver mesh.

Official [nonlinear settings](https://wec-sim.github.io/WEC-Sim/main/user/advanced_features.html#nonlinear-buoyancy-and-froude-krylov-excitation)
exclude user-defined irregular waves. The documented method changes buoyancy and
Froude-Krylov pressure on instantaneous geometry while retaining linear wave
theory; it is weakly nonlinear, not a complete nonlinear fluid solution.
[Theory](https://wec-sim.github.io/WEC-Sim/main/theory/theory.html#nonlinear-buoyancy-and-froude-krylov-wave-excitation).

Two bounded routes, both needing consumer ownership and independent review:

1. Re-express the identical deterministic incident spectral bank through a
   supported native irregular-wave path with explicit amplitudes/phases, then
   verify eta/pressure/excitation and actual-rest startup equivalence before
   paired linear/nonlinear runs. Existing directional/nonlinear paths must be
   exercised; do not assume representation alone makes them correct.
2. Add a narrowly reviewed imported-incident field/pressure interface to native
   nonlinear buoyancy/FK evaluation, using the same spectral bank everywhere.
   This is a maintained-fork extension, not a Topwater implementation. Geometry
   clipping/integration and recovery of baseline linear limits require regression.

For either route, a proposed final operational criterion is all six nonlinear
versus linear histories within 2 percent of each informative scaled response
peak AND within 1 mm / 0.01 degree, with nonlinear spatial/temporal integration
uncertainty smaller than those differences. These are proposed working limits,
not established physics constants or already authorized project requirements.
They reserve roughly one tenth of the historical comparison budgets. Weak axes
need explicit absolute-error reporting. A weakly nonlinear pass supports only
the nonlinear mechanisms modeled; it does not certify viscous or fully nonlinear
diffraction/radiation effects. Existing reference mesh/time uncertainty remains.

## Startup, ownership and stopping point

Scaling the accepted packet is not an accepted new startup. Even factor two
would make its direct initial heave about 0.1486 N and native pitch about
0.1123 Nm, exceeding the existing 0.1 cap. Any actual new input requires separately
frozen upstream placement/timing, full-region calm and force-tail checks, actual
zero body initial state, physical band and useful arrival. Preserve old evidence.

Recommended next bounded assignment after minimum-case closure: Capytaine owns
the geometric rejection screen; WEC-Sim owns incident equivalence and eventual
paired weakly nonlinear dynamics. Coordinator freezes criteria/resources and
serial execution lease. Stop at the first rejected amplitude or strongest
screen-surviving candidate and its explicit limitations; no new sea is qualified
until the chosen nonlinear confirmation exists. Topwater remains held.
