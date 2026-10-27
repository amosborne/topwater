# Issue 31: validation-reference tail accounting

This is an independent method review, dated 2026-10-27. The coordinator selected
this bounded proposal after the seven-probe input-weighted screen. It is not a
physical correction, a numerical run, or acceptance of the whole surface.
Capytaine expert owns the implementation and its tests.

The user requires 10 mm maximum numerical elevation error, including reference
uncertainty, against the chosen linearized model. Runtime remains the intact
complex response applied to the saved velocity with its declared interpolation.
The following smoother input and endpoint separation are **validation only**.

## Declared reference

Interpolate every saved velocity and acceleration with a piecewise cubic Hermite
velocity. Keep zero velocity before time zero and the declared cubic continuation
from the final saved velocity/acceleration to final rest. Thus velocity and its
first derivative are continuous at internal knots and the future join. Initial
velocity must be exactly zero. Its right derivative at zero, a0, may differ from
the zero derivative in the preceding rest interval.

Use the previously reviewed finite-low discrete endpoint A, with native
s = -i omega, and define R = H - s A = H + i omega A. Evaluate the full complex
band-limited R on the Hermite input, then add A times the saved acceleration at
the original output knots. Do not add A to the production full-H runtime.
Measure the actual runtime-versus-reference difference at every output point
and time; neither the earlier 0.287 mm observation nor seven-probe results can
replace that full-domain measurement.

At t0, saved acceleration specifies the right-limit convention. If a0 is nonzero,
a symmetric Fourier derivative instead takes half that derivative jump exactly
at zero. State this convention, retain the runtime's raw initial value, measure
the discrepancy and preserve the existing physical startup gate. At internal
Hermite knots and the final saved-time join, the derivative is continuous.

## Independent tail derivation

Let V(omega) = integral v(t) exp(-i omega t) dt for the compactly supported,
piecewise cubic Hermite velocity. Distributional integration by parts gives

```text
V(omega) = a0/(i omega)^2
         + sum_j jump(v'')(tj) exp(-i omega tj)/(i omega)^3
         + sum_j jump(v''')(tj) exp(-i omega tj)/(i omega)^4.
```

The jump sums include the initial transition from rest, every internal native
knot, the saved-time/future-polynomial join, and the final transition to rest.
There is no missing velocity-jump term because v is continuous, and no other
first-derivative jump because the selected interpolant and continuation are C1.
Coefficients and signs follow the input transform exp(-i omega t); absolute
values remove sign dependence from the bound.

Define J2 = sum absolute jumps of v'' and J3 = sum absolute jumps of v''',
separately for each of the six input modes. If the declared engineering envelope
|R(point, mode, omega)| <= C(point, mode) holds for **all** omega above Omega,
then the omitted regular elevation is bounded uniformly in time by

```text
U_tail(point) <= sum_mode C(point,mode)/pi *
    (abs(a0_mode)/Omega + J2_mode/(2 Omega^2) + J3_mode/(3 Omega^3)).
```

The 1/pi factor accounts for both conjugate halves of the inverse transform.
Evaluate at Omega = 48 rad/s to cover all omitted regular frequencies, including
those above the FFT Nyquist frequency. No divergent absolute estimate for the
known linear-in-frequency A term is needed: that term is treated analytically
in this validation reference. Sampling, padding, retained-band frequency
interpolation, future and input-interpolation controls remain separate.

This inequality is conditional on C. The proposed C = 2 max_known |H - s A|
does not prove absence of an unseen resonance or faster growth. Calculate C
from each actual exterior query's authorized data; do not copy the seven-probe
maximum or an old 0.5 mm allowance over the grid. Sparse 64/96 evaluations may
enter only if separately allocated. Chosen operator/model accuracy remains
distinct from reconstruction accuracy and the deferred mesh review.

## Prospective implementation controls

Use the compact cubic input v(t) = t (1-t)^2 on [0,1], zero elsewhere. It has
a0 = 1; second-derivative jumps -4 at zero and -2 at one; third-derivative jumps
+6 at zero and -6 at one. Hence J2 = 6 and J3 = 12, and

```text
V(omega) = 1/(i omega)^2
         + (-4 - 2 exp(-i omega))/(i omega)^3
         + (6 - 6 exp(-i omega))/(i omega)^4.
```

Compare this expression to independent real/imaginary polynomial quadrature at
predeclared moderate, nonzero frequencies, for example 0.7, 1.3 and 2.7 rad/s.
Verify the computed jumps, endpoint inclusion and 1/pi integrated tail formula.
Also use a multi-interval C1 cubic control with an internal curvature jump and
a cubic future join, so passing does not rely only on a one-piece example.
Zero and scaled inputs check vanishing terms and linear scaling. These are
input-only controls; no numerical experiment was performed by this reviewer.

Stop at the independently reviewed producer and actual-result evidence. Whole
5929-by-751 acceptance still requires pointwise measured error plus the full
uncertainty accounting to stay within 10 mm, with unchanged incident/diffraction,
motion, exterior continuation and physical startup/contact gates.
