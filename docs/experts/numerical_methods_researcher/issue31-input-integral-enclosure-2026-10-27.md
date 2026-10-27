# Issue 31: tighter enclosure of the known input tail

Independent method proposal, 2026-10-27. The original full-grid attempt stopped
after 128 queries: its combined estimate was 12.3233 mm, dominated by an
11.7158 mm input-tail triangle estimate. The actual candidate/reference difference
on those queries was only 0.00922 mm. See the versioned
[failure review](../../../artifacts/workflow/issue31-full-review1/acceptance.json).
This supports tightening the integration of the **known input**, without changing
the user 10 mm criterion, response envelope C, hydrodynamics, motion or runtime.
The solver owner authors implementation; this reviewer performs no numerical run.

## Finite interval enclosure

Use the exact Hermite input transform V derived in
[the preceding note](issue31-hermite-tail-2026-10-27.md). For a finite interval
[A,B] with maximum node spacing h and M bounding the magnitude of V'' throughout
that interval,

```text
integral_A^B abs(V(w)) dw
    <= trapezoid(abs(V at nodes)) + (B-A) h^2 M / 12.
```

The complex linear-interpolation error on an interval is at most
M (w-left)(right-w)/2. Its integral is M h^3/12. Convexity bounds the integral
of the magnitude of the complex linear interpolant by the trapezoid of its
endpoint magnitudes. Summing yields the stated upper enclosure; this is not
an unjustified trapezoidal estimate of a rapidly oscillatory magnitude.

For one transform term c exp(-i w t)/(i w)^p,

```text
abs(term'') <= abs(c) *
    (t^2/w^p + 2 p abs(t)/w^(p+1) + p(p+1)/w^(p+2)).
```

The initial a0/(i w)^2 contributes 6 abs(a0)/w^4. Sum the p=3 and p=4
jump contributions separately for each mode. Evaluating these nonnegative
decreasing expressions at A gives a valid uniform M. Add the original analytic
jump-tail integral above B. Finally multiply each mode's enclosed input integral
by C(point,mode)/pi and sum modes, exactly once. The response-envelope assumption
remains conditional; the tighter input calculation does not prove that C holds.

## Cheap evaluation and coverage

The coordinator should freeze the finite upper target and frequency resolution
before new results. B near 512 rad/s and a Fourier period 4096 or 8192 s were
proposed; these are input-quadrature settings, not new BEM frequencies or a new
physical time history. Compute once for the six inputs and reuse across queries.

For h=2 pi/T, include floor(48/h) h through ceil(512/h) h, or evaluate the two
exact endpoints explicitly. The former slightly overcovers the omitted band;
use its actual lower endpoint in M and its upper endpoint for the analytic
remainder. Never omit the small edge intervals because 48 and 512 are not bins.

Jump sums on the native 0.02 s grid can use FFT numerators. These finite
exponential sums are periodic in frequency even beyond native Nyquist. If an
rFFT is used, indices above its half range require conjugate folding. The
denominators (i w)^p always use the actual, unwrapped w. A frequency shift of
one native Fourier period does not make the whole input transform periodic.

If actual knot times differ from their rounded grid values by delta_t, the
extra finite-interval input integral is bounded by

```text
sum abs(jump2 * delta_t) * (1/A - 1/B)
+ 0.5 sum abs(jump3 * delta_t) * (1/A^2 - 1/B^2).
```

This follows from abs(exp(-i w t)-exp(-i w rounded_t)) <= w abs(delta_t).
Retain explicit numerical-evaluation/roundoff checks; an analytic interpolation
bound alone does not certify floating-point transform evaluation.

## Prospective controls and stopping point

Retain the original exact cubic pulse, multi-interval C1/future-join, zero and
scaling controls. Add direct exponential sums at predeclared FFT bins below and
above native Nyquist and above one complete native frequency period. Check
conjugation, actual denominators, endpoint coverage, nonnegative remainder,
and any time-rounding allowance. A finite-interval known-polynomial transform
can independently check the enclosure against adaptive magnitude integration;
the interpolation remainder, rather than apparent agreement, justifies its bound.

Preserve the failed triangle-bound attempt. Source and actual-control review
precede a separately allocated retry. The first-chunk accuracy/resource stop,
common-field identities, whole-domain coverage and physical startup/contact
gates remain unchanged. No threshold relaxation, fitted correction, response
frequency campaign or readiness promotion is part of this method proposal.
