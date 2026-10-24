# Mild periodic case: prospective reference protocol

Scope: freeze the controlled sea and one-second comparison, implement/review
periodic motion reconstruction and its spatial uncertainty using complete raw
Capytaine responses, and specify matching body/fluid initialization. No old
failed study is relabeled; no Topwater runtime substitute or speed campaign.
The overarching practical time-domain goal remains incomplete.

The case is declared in validation/mild-periodic-case.json BEFORE computing its
acceptance. Three waves span the existing 0.10/0.25/0.50 Hz test frequencies and
0/45/90 degree headings, with amplitudes 20/10/5 mm and fixed phases 0/pi/2/pi.
This is a controlled directional sea, not a renamed observed spectrum. These
frequencies exercise the supplied finite-wave responses without interpolating
a broadband radiation kernel. Phases/amplitudes will not be fitted to Topwater.

## Reference scope and acceptance

The historical raw study is an immutable INPUT with a failed coefficient-policy
decision. A new output-specific protocol is justified by the user's explicit
0.1 degree / 1 cm motion target. It does not qualify every coefficient for every
possible forcing. Retain all 36 raw matrix entries and their uncertainties.

For each selected wave, reconstruct Re(H*a*exp(i*phase-i*omega*t)); derivatives
follow analytically. Sum all three waves without fitting phase or amplitude.
For each physical DOF, sum the absolute complex mesh-increment envelopes over
waves. This bounds the difference of the represented periodic histories for all
times, avoiding cancellation in a short sample window. It is empirical mesh
uncertainty, not a certified continuum bound. Require both finest increments
within 10% of the user's motion limit and the non-increasing component trend
within the declared roundoff floors. The 10% allocation is an engineering reserve
chosen prospectively for the reference; final error plus ALL relevant uncertainty
must still fit the full user limit. Keep original coefficient gates unchanged.

Report exact wave frequencies/directions/amplitudes/phases, equivalent Hm0 from
the component variance, all six motion ranges and velocities, and which motions
are meaningfully excited relative to the target. A negligible-forcing or quiet
startup interval cannot establish a useful baseline. A constant-initial-pose
negative control must exceed the user's tolerance in at least one translation
and one rotation after subtracting twice the reference spatial/sampling envelope,
before accepting the case as informative. This does not claim
that the sea alone validates each weakly excited mode. Sample-step refinement of
the analytic reference is a comparison/sampling check, not an ODE accuracy claim.

## Full initialization is mandatory

The reference is periodic with prehistory, not a quiescent release. Initialize
the body with its reconstructed pose/velocity and the fluid with the total
incident, diffraction and body-radiation periodic elevation/potential at t=0.
World rotation coordinates map to physical orientation by the exponential map,
with angular velocity from its consistent derivative. A body-only match is
insufficient. Linear mean-plane fields must be sampled/pulled onto the nonlinear
surface and projected to the current body boundary condition; quantify that
representation/compatibility error, including wetting changes near the chine.
Periodic prehistory is exact for the linear reference. The projected nonlinear
initial state can launch adjustment transients, which count against the same
budget. The rotation chart is frozen as theta=Log(R*R_equilibrium^T) in mean
world axes; compare those components, not a differently ordered Euler chart.
The nonlinear grid-to-reference geometry map remains a blocking field-export
task: reject points lying inside the mean hull or outside the reference fluid
domain instead of silently using an interior continuation of the potential.
Validate the map/projection in motion units, not just a small boundary residual.

Future reference body trajectories, pressures or forces cannot drive Topwater.
Only prescribed analytic incident waves and explicit physical initial-state data
may enter its own time-domain fluid solver. The historical Cummins array math
remains useful for later causal tests; this discrete periodic case does not
require an invented broadband kernel, frequency tails or zero-memory startup.
Analytic acceleration and future phasor/motion arrays are reference diagnostics;
the runtime packet carries only initial pose/velocity and sampled eta0/phi0 plus
incident waves. Topwater must calculate its own pressure and acceleration.

## Deliverables and stopping checkpoint

Archive the case, raw data/producer identity, reconstruction code, histories,
complex component amplitudes/derivatives, empirical uncertainty and separate
qualification flags. Initial periodic fluid sources/field data and their
uncertainty are not yet exported: that missing evidence MUST block complete
reference readiness even if the motion-only spatial gate passes.

Stop this checkpoint after the case/motion artifact and independent review,
with an exact next action for full periodic field initialization and minimum
Topwater TD physics/cost screening. Commit/push, update issue #30 and hand off.
No passing reference or actual simulation claim follows from a plan or traces
alone. Jobs remain serial with archived immutable producers and fresh outputs.

## Candidate 001 result

Source checkpoint a5ebf7e. The declared array-only candidate completed; overall
exit 1 / passed=false is required because periodic fluid fields, field uncertainty,
nonlinear compatibility and linear-limit checks are missing. No Topwater TD
simulation or pressure solve ran.

The controlled sea's equivalent component Hm0 is 0.0648074 m. The spatial+sampling
gate passed with non-increasing component increments, using at most 44.10% of
the allocated 10% reference budget. This is not full reference qualification.
The conservative frozen-pose negative control fails the user's tolerances in
surge, heave, roll and pitch, so the one-second interval is informative. Sway
and yaw are weaker; this sea alone does not establish every mode's accuracy.
Six-mode/coupling checks remain required before baseline promotion.

[Complete candidate/producer/input bundle](../validation/evidence/mild-periodic-001.zip)
includes reference histories/derivatives, all raw matrix entries/increments,
the body/wave-only initial packet and nested complete raw-reference evidence.
[Summary and SHA256](../validation/evidence/mild-periodic-001-summary.json).
The original study 003 decision is unchanged and failed; it is retained as raw
input to this separate prospective, forcing-specific output assessment.

The focused validation suite now passes 60 tests, including harmonic signs and
derivatives, coupled-component uncertainty without cancellation, reference
sampling allocation, exact rotation-map/body-rate consistency, conservative
signal checks and rejection of changed case/readiness conventions.

Independent review rebuilt the phasors, derivatives, uncertainty sums, sampling
bound and conservative signal checks, and verified body angular velocity by
independent rotation finite differences. Provenance and the nested raw dataset
were verified; all four missing-readiness flags and overall passed=false remain
correct. No blocking finding remains for this limited stage.

Next bounded action: freeze the minimum Topwater TD fluid grid/domain and a safe
reference-to-nonlinear geometry map; export the required total periodic eta0/phi0
fields with field/mapping uncertainty. Reject invalid reference-fluid points.
Only after a complete compatible initial-state packet exists, run a bounded
preparation/actual coupled-RHS cost screen before a one-second simulation. No
Capytaine trajectory/pressure/acceleration may drive the runtime.
