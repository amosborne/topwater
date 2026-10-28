# Per-carrier packet construction: corrected definition

**Superseded study warning.** The shared-envelope construction below is a
historical diagnostic under its original policy. It conflated the carrier sea
grid with a permitted packet-sideband dictionary and coupled different carriers
through one envelope. Its reported55s/12s ramps are not general requirements for
individual packets. Preserve that evidence; do not publish its figures as the
solution to the user's corrected per-carrier problem.

## Current definition

The carrier sea grid chooses the original frequencies and their amplitudes,
phases and directions. Each original carrier then gets an independently
constructed **unit packet** determined by its base frequency f0 and requested
dwell D, under a fixed documented design policy. All packets share the same
physical dwell start and duration; their ramps and spectral offsets may differ.

Use dimensionless time `s=f0*(t-t_dwell)` and dwell length `N=f0*D` cycles.
A useful representation is

`q_N(s)=g_N(s)*exp(-i*2*pi*s)`, with `g_N(0)=1` and `g_N approximately1` on `[0,N]`.

The original amplitude and specified phase multiply this unit packet afterward.
Thus, at the reference position and dwell start, every carrier has its specified
phase at full envelope amplitude. Direction enters the physical travelling-wave
components and dispersion; it need not change the scalar unit-envelope design.

For dimensionless offsets xi, actual component frequencies are `f0*(1+xi)`.
A symmetric fractional band `|xi|<=beta<1` keeps all physical frequencies
positive. The offset spacing is a packet-construction resolution, **not** the
spacing between original sea carriers. A separate hard bound on all final
physical components, if required, must be an explicit policy constraint rather
than silently inferred from the carrier grid's endpoints.

The pair `(f0,D)` does not mathematically select one spectrum: many coefficient
sets meet the same dwell. A fixed policy must also specify relative bandwidth,
dwell/quiet tolerances, overshoot/cancellation bounds, quiet horizon, spectral
resolution and the optimization objective. A deterministic algorithm then makes
the construction unique and reproducible. For example: fix those constraints,
find the shortest accepted ramp to a declared precision, then use a strictly
convex coefficient-norm tie-break on that feasible dictionary. This is an
algorithmic definition of one packet, not a proof that no other policy is better.

## Which parameter controls what?

- f0 sets the carrier period and, with dispersion, its wavelength. D does not
  set the minimum physical frequency.
- D controls the plateau length and a central envelope spectral scale of order
  `1/D`; this does not replace the carrier frequency.
- The transition duration controls the additional offset bandwidth needed for
  a given accuracy. Relative quantities are `r=f0*R` and `beta=B/f0`.
- Uniform packet spacing `delta_f` implies a recurrence period `1/delta_f`.
  Refining it chiefly changes that horizon/resolution; it is not itself the
  ripple control. Ripple and overshoot depend on the finite bandwidth and chosen
  coefficients, and should be constrained explicitly.

Finite discrete spectra cannot provide permanent calm. The design must state
how long the quiet interval is required to last, and the full simulation must
finish before packet recurrence or use a policy with a sufficiently long quiet
horizon. Unknown settling time is not supplied by f0 and D alone.

## Cheapest next demonstration to propose

Keep the100s dwell and use independent examples f0=0.03,0.10 and0.49Hz, giving
dimensionless dwells of3,10 and49cycles. Begin with a bounded scalar design,
not an optimization of every sea component. An illustrative fixed fractional
band beta=0.5 on a128cycle Fourier period would use K64,65 real cosine
coefficients, for each case. This is a **proposed comparison policy**, not an
adopted bandwidth or physical quiet horizon. It makes the three small problems
comparable without forcing them to share one envelope. Normalized solutions can
be cached by `N` and policy; amplitude/phase/direction are applied afterward.
That cache reduction applies only to a scale-invariant normalized policy; fixed
Hz bounds or fixed-second quiet intervals add their normalized values to the key.

First compare the resulting target envelopes, discrete spectra, dwell errors
and peak quiet tails, with specified phases referenced to the common dwell
start. Then verify positivity of all physical frequencies and convert normalized
ramps to seconds with `R=r/f0`. Only a later propagated-field check can establish
whole-hull startup/ending calm. No new implementation or numerical job is
authorized by this conceptual correction.

---

# Historical shared-envelope proposal — original policy preserved below

This is an executable diagnostic proposal for [issue36](https://github.com/amosborne/topwater/issues/36),
not a production input change or a solved hydrodynamic case. The requested
physical frequency grid is **0.01–1 Hz in 0.001 Hz increments** and the common
full-amplitude dwell lasts **100 seconds**. Specified carrier phases apply at
the common dwell start and fixed spatial reference point.

The helper is `artifacts/workflow/wave-grid100-design1/solve.py`, with a fixed
embedded design plan and only `--output FRESH_DIRECTORY` required. The coordinator
owns execution; this note does not claim any results before that run.

## Explicit provisional acceptance criteria

For the demonstration, use at most 1% complex-envelope error during dwell and
at most 1% envelope magnitude during the declared quiet interval. Set a 1.05
upper envelope bound and coefficient l1 bound of 10. These are **coordinator
study assumptions, not user-adopted physical thresholds**.

For a unit 1 m carrier, 1% corresponds to a 1 cm envelope-modulation bound. It
does not establish 1 cm total error for multiple waves: a conservative combined
bound is `0.01 * sum(carrier amplitudes)`. Actual hull-region calm, fluid speed,
forces, radiation reconstruction and settling require separate checks.

The Fourier period is 1000 s, so offset frequencies are integer multiples of
0.001 Hz. Quiet is checked throughout the rest of that period outside the
transition gaps; this does not create a permanently isolated pulse.

## Two common-envelope problems, not arbitrary single-carrier bandwidth caps

Use a real, symmetric envelope around the dwell centre:

`g(u) = c0 + sum(k=1..K) ck*cos(2*pi*k*u/1000)`.

The dwell is `|u| <= 50 s`. A ramp/transition gap R occupies the next R seconds
on each side. Quiet is `50+R <= |u| <= 500 s`. Impose `g(50)=1`; symmetry also
anchors the other dwell endpoint. A real positive envelope during dwell preserves
the specified original carrier phases.

Symmetric sidebands must fit on both sides of **every** original carrier.
Therefore the available offset bandwidth for a common envelope is the smallest
carrier distance to either physical band edge:

| Original carriers sharing the envelope | Full available symmetric offset band | K | Real cosine coefficients |
|---|---:|---:|---:|
| 0.03, 0.10, 0.25, 0.49 Hz | 0.02 Hz | 20 | 21 |
| 0.10, 0.25, 0.49 Hz | 0.09 Hz | 90 | 91 |

These choices use the full symmetric bandwidth available to each **set**. They
do not find the shortest individual packet for 0.25 or 0.49 Hz, which have more
room when considered separately.

At either endpoint carrier, 0.01 or 1 Hz, this symmetric construction has K=0.
Its envelope is constant and cannot satisfy both the dwell anchor of one and
quiet magnitude at most 0.01. Consequently no common envelope of this class
covers all endpoint-inclusive carriers. This is a restricted-class result, not
a proof that asymmetric complex envelopes or other physical inputs are impossible.

## Search and verification

Each LP minimizes a shared peak error delta on the dwell and quiet constraints.
The transition is free apart from a lower bound `-delta` and a guarded sampled
upper cap of 1.049. Auxiliary variables enforce the coefficient l1 limit.

Search integer-second R between 1 and 400 by bisection. A trial is classified as
infeasible only when its numerical dual lower bound exceeds 0.01; solver failure
or uncertainty is never treated as infeasibility. A finite-box Lagrangian
correction accounts for residual stationarity in returned dual multipliers,
with an additional numerical allowance. It is an engineering certificate, not
floating-point interval arithmetic.

The candidate envelope is evaluated densely by FFT. The conservative missed-peak
allowance is `max|g''| * h^2 / 8`, plus numerical allowance, using
`max|g''| <= sum |ck|*(2*pi*k/1000)^2`. All dwell/quiet boundaries lie on the
dense grid. Acceptance includes these allowances, the actual 1.05 global maximum,
lower bound, phase anchor and coefficient norm. An uncertain trial may add a
bounded set of dense constraint sites, or stop as inconclusive.

The search uses the slightly stricter 1.049 sampled cap to reserve room for the
continuous 1.05 limit. A final solve at the preceding integer uses the relaxed
1.05 cap. Only if its lower bound also excludes 0.01 can the accepted result be
called the shortest integer ramp in the full stated symmetric class. Otherwise
report the found acceptable ramp and restricted search bracket, without a
certified shortest claim.

## Execution bound and deliverables

There are two searches, each limited to 12 LP solves, six seconds per solve:
24 solves maximum and at most 144 s of LP time. The proposed root guard is
240 s, 1 GiB, one thread. The largest LP has 183 variables and initially 2185
inequalities; dense validation uses FFTs rather than a large cosine matrix.

Outputs are the embedded plan JSON, actual search history/metrics, numerical
coefficients and traces, and time/discrete-frequency plots. Representative unit
carriers illustrate the common envelopes; the actual sea is not remapped and
no BEM, WEC-Sim, propagated-wave or production validation is performed. No image
inspection or preview campaign is part of this diagnostic.
