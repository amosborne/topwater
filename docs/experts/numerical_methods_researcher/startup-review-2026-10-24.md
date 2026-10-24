# Minimum-case startup review, 2026-10-24

Assignment: independent read-only startup/padding review for the authorized
minimum-mesh reference. No solver edits or numerical jobs. This note records
review reasoning and proposed independent expectations, not executed validation.

## Source identity and observations

Read-only HEAD and status checks of
`artifacts/dependencies/topwater-wec-sim` returned clean official v7.1.0 commit
`a23d2993807410fe78cc63546c1ba9d81912024f` before owner implementation.

- `source/functions/initializeWecSim.m:255`: wave setup receives `simu.time`.
- `source/objects/waveClass.m:874-888`: `waveElevUser` interpolates the imported
  data onto only that timeline and applies the native ramp.
- `source/objects/bodyClass.m:1086-1105`: `userDefinedExcitation` interpolates the
  excitation kernel and uses `conv(...,'same')*dt` on that trimmed wave array.
- `source/functions/BEMIO/excitationIRF.m`: BEMIO creates a symmetric time grid
  `linspace(-tEnd,tEnd,nDt)` and integrates the complex excitation coefficients.

Thus extra CSV rows outside `simu.time` do not reach the convolution. This is a
source-proven missing-input path; an executed independent reproduction remains
required before relying on a correction. `same` alignment is suitable only when
the kernel zero lag is the central sample; use odd, symmetric, commensurate grids
or explicitly account for physical lag indexing.

The selected source's theory also states that radiation memory assumes no body
motion for negative time. Excitation history and body/radiation history are
different: padding the former must not initialize motion by cropping an already
moving body trajectory.

Primary explanatory documentation:
[WEC-Sim convolution theory](https://wec-sim.github.io/WEC-Sim/main/theory/theory.html#convolution-integral-formulation),
[imported elevation](https://wec-sim.github.io/WEC-Sim/main/user/code_structure.html#elevationimport).
Those live pages explain the method; selected local source above determines the
version-specific implementation observation.

## Physical boundary and acceptance

Mathematical inference: a nonzero finite-band packet cannot have exactly compact
spatial support. Setting one elevation sample to zero does not set the velocity
potential or fluid velocity to zero. A compact initial elevation with one-way
deep-water dispersion has nonlocal potential/velocity; a sharp truncation also
introduces broad frequency content. No spectral or excitation tail may be hidden
by overwriting initial samples or by a force switch.

The coordinator accepts numerical calm with explicit quantified tails. A
defensible minimum case must therefore retain the canonical carrier amplitudes,
frequencies and directions, declare its spatial envelope and resulting bandwidth,
and report initial hull-region elevation and fluid velocity plus all six initial
excitation forces/torques. Exact body displacement, velocity and radiation memory
are zero at the actual start. Meaningful incident arrival and motion must remain
inside the first ten seconds. No temporal ramp is introduced.

For the declared convention
`F(t) = integral K(tau) eta_ref(t-tau) d tau`, output `[0,10]` and retained kernel
`[tau_min,tau_max]` require input `[-tau_max,10-tau_min]`. For symmetric `[-T,T]`
this is `[-T,10+T]`. Padding is incident input, not body warmup. Values must come
from the same dispersive incident field used by the renderer, at its stated
wave-reference origin; an independently moving visual wavefront is insufficient.
Noncausal force before origin elevation arrival can be physical for an extended
body and must not be clipped.

Prospective bounds adopted by the coordinator before packet diagnostics:

- At time zero over the hull region, `max(abs(eta)) <= 1e-3 sum(a_j)`:
  `3.5e-5 m` for the canonical amplitudes.
- At time zero over the hull region, fluid speed magnitude is at most
  `1e-3 sum(a_j omega_j) = 4.39823e-5 m/s`. This carrier-derived incident-velocity
  scale is a declared reference scale; actual broad-band velocities must be
  evaluated, not estimated from envelope elevation alone.
- Per translational/rotational DOF, `abs(Fexc_i(0)) <= 1e-3 S_i` and also at most
  `0.1 N` / `0.1 Nm`, where
  `S_i = sum_j a_j abs(H_i(omega_j,theta_j))` uses dimensional native excitation
  at the canonical carriers. A zero scale requires an explicit zero-force check,
  not an invented denominator floor.
- At least one incoming packet must reach 25 percent of its canonical amplitude
  at the hull within ten seconds, with nontrivial resulting body response.
  Exact body displacement/velocity and radiation memory remain zero initially.

These quantify the initial-calm approximation, not reference motion accuracy.
They prevent choosing tolerances after observing a particular packet. Packet
geometry may be adjusted before freezing the physical forcing; tails may not be
clipped and motion histories may not be fitted.

Feasibility is unproven pending incident-only diagnostics. The supported
wavelength minimum `2*beam = 5.4864 m` implies approximately `0.533 Hz` in deep
water, stricter than the separate `0.6 Hz` ceiling. The 0.5 Hz carrier therefore
has little available envelope bandwidth and can require a long incident packet.
Keeping its peak upstream beyond ten seconds is consistent with the stated
at-least-one-arrival gate, while preserving its canonical carrier amplitude.
The 0.25 Hz component has more bandwidth for the first visible arrival. A smooth
spectrum bounded away from zero avoids introducing a poorly controlled low-k
tail; localization and velocity still require actual diagnostics. No assertion
of a feasible packet or simulated motion follows from this reasoning.

## Independent cheapest rejection expectations

No execution performed. Test with `dt=1`, physical lags `tau=-2:2` and
`K=[1,2,4,8,16]`; use direct indexed convolution as the independent expectation.
Output times are `0:10`.

1. Unit elevation impulse at time 11: `F(9)=1`, `F(10)=2`. The unmodified selected
   code drops that future input and returns zeros, exposing the right-edge issue.
2. Unit elevation impulse at time -1: `F(0)=8`, `F(1)=16`, exposing missing past
   support.
3. Unit elevation impulse at time 5: `F(3:7)=[1,2,4,8,16]`, confirming time/sign
   alignment away from edges.
4. Preserve output wave times/elevations and existing ramp behavior; padding
   should affect only excitation preparation, with each direction kept separate.
5. Reject insufficient history rather than silently infer zeros, unless an
   explicitly declared zero-extension policy is independently justified.

Subsequent packet check: compare imported elevation with the renderer's field at
the reference point, and excitation with an independent retained-lag sum on the
same times. Report initial tails and onset from the unmodified signals.

Stopping point: owner proposal and ensuing padding correction review; no broad
formulation research, mesh study or implementation ownership.

## Executed stock evidence and correction source review

Independently inspected
`artifacts/reference-imported-excitation-stock/result.json`: the three actual and
expected vectors agree exactly with the independent endpoint predictions above,
with maximum errors `[2,16,0]`. This confirms the selected scalar native crop
defect, rather than merely a source suspicion.

Read-only review of the owner's initial `bodyClass.m` correction finds the
explicit query `output_time - physical_lag` and `sum(K*eta)*dt` correct in sign,
index and scaling. It retains raw file history for excitation, leaves the output
wave array unchanged, rejects missing support/nonfinite interpolated kernels,
and preserves the declared ramp on the queried physical timeline. The draft
`check_imported_excitation.m` uses independent nested sums for time steps 1 and
0.5, future/past/interior samples, short support, ramp and output-grid invariants.
This is source review, pending executed fixed-regression evidence.

The separate draft `reproduce_multiple_imports.m` uses native `hydroForcePre`
with radiation alone stubbed, two input amplitudes 1 and 2, and expected `3*K`.
The stock first-wave-only path predicts missing `2*K`, maximum error 32. It is
adequate to test omitted component handling, but its equal headings do not test
direction dispatch. The fixed regression must include unequal direction kernels
and distinct headings as well as the actual production preprocessing call.
