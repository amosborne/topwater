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
