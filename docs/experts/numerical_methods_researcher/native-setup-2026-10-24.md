# Native setup and focused reference question, 2026-10-24

Native agent identity: `/root/researcher_setup`; exact registry role:
`numerical_methods_researcher`. The role's developer instructions were delivered
on activation, independently of reading its TOML: own bounded evidence-backed
method proposals; coordinator selects experiments; reference experts own their
implementations; no numerical job without a lease; no formulation/optimization
campaign or solver edits; independent review for physical corrections; inherited
parent settings; durable context in files rather than permanent model learning.

Scope/acceptance: create reusable proposal template and one focused reference
question, record source identity and stop. No literature search, numerical
experiment, MATLAB invocation, test, solver/model job, code change or lease
request was made. No physical result is established by this setup.

## Exact paths read

- `.codex/agents/numerical_methods_researcher.toml`
- `docs/experts/numerical_methods_researcher/README.md`
- `docs/experts/numerical_methods_researcher/current-state.md`
- `docs/experts/shared-contract.md`
- `docs/experts/coordination.md`
- `docs/experts/capytaine_expert/native-setup-2026-10-24.md`
- `artifacts/dependencies/wec-sim-v7.1.0/source/objects/bodyClass.m`
  (targeted search and `userDefinedExcitation`, lines 1086-1107)
- `artifacts/dependencies/wec-sim-v7.1.0/source/objects/waveClass.m`
  (targeted search and `waveElevUser`, lines 874-888)

Read-only `git rev-parse HEAD` for that official checkout returned
`a23d2993807410fe78cc63546c1ba9d81912024f` (selected v7.1.0).
Only this note, `proposal-template.md` and own `current-state.md` were changed.

## Active question for the reference owners

What exact wave-reference location, directional incident histories and retained
convolution timeline allow the first ten seconds to start from rest in locally
calm water, while the actual excitation and displayed upstream incident front
share the same arrival, with no imposed temporal amplitude ramp?

Existing official selected-source anchors:

- [bodyClass.m at selected commit](https://github.com/WEC-Sim/WEC-Sim/blob/a23d2993807410fe78cc63546c1ba9d81912024f/source/objects/bodyClass.m#L1086-L1105):
  `userDefinedExcitation` explicitly uses a noncausal excitation convolution;
  it samples the kernel on `min(kt):dt:max(kt)` and uses `conv(...,'same')*dt`.
- [waveClass.m at selected commit](https://github.com/WEC-Sim/WEC-Sim/blob/a23d2993807410fe78cc63546c1ba9d81912024f/source/objects/waveClass.m#L874-L888):
  imported elevations are interpolated onto the supplied simulation timeline
  and multiplied by its ramp function. Extra input-file rows alone therefore
  do not establish retained convolution prehistory/future support.

These are source observations, not an independently reproduced defect. The
mapping between kernel time zero and discrete `same` alignment remains to be
established by the WEC-Sim owner, especially for asymmetric kernel support.

Analytical question to settle before an execution assignment: under the stated
convention F(t)=integral K(tau)*eta_ref(t-tau) d tau, a retained kernel on
[tau_min,tau_max] requires eta_ref on [-tau_max,10-tau_min] for t in [0,10].
The consumer must establish the actual sign/reference/alignment before using
this interval. Negative-time support can make force precede arrival at the
chosen reference point; local zero elevation alone does not prove zero load or
compatible body rest. Finite-band tails and their physical significance are
unproven project questions, not setup acceptance thresholds.

Requested bounded answer from WEC-Sim/BEMIO owner, jointly checked by Capytaine
owner: kernel support/alignment and interpolation contract, per-direction
reference histories, actual retained padding, earliest excitation versus local
surface arrival, and initial radiation/body-memory conditions. Preserve the
actual rest-time origin; do not crop an already moving run. Agree `rampTime=0`
through the native consumer path and use the same prescribed incident field in
rendering. A field display is not a computed scattered-fluid solution.

Cheapest subsequent rejection plan (not executed or authorized here): first
derive an asymmetric synthetic-kernel convolution with independent index math;
then compare native preprocessing under a later sole execution lease, including
retained timeline edges and sign/direction cases. Physical hydrodynamic runs
follow only a separate bounded assignment. Consequential corrections require
independent review. This question does not introduce a coefficient gate or a
whole-library audit.

Next action: wait for a focused reference-owner question; numerical work remains
unassigned and unleased.
