# Issue31: bounded radiation-surface diagnosis

2026-10-27. Capytaine expert source/data-inventory review only. No numerical job,
BEM solve, physical correction or public-artifact change was performed here.
Current user authorization is resolve issue31; it supersedes the completed
adoption goal in the older handoff. Root owns allocation and acceptance.

## Current evidence and conventions

Active branch checkpoint55e44a56281724a89bfb309eb7b28031ca5492b7 uses surface
source SHA256240198ee21b92334b2b195adc4e7d5a31ff347f9c0320a95d16569e832520251.
Installed distribution metadata identifies Capytaine3.0.0. Production constructs
modal-velocity transfer `Hv=i*eta_q/omega`, uses a real-part cosine inversion,
and recovers the response with `exp(+i*omega*lag)`. These agree with Capytaine's
`exp(-i*omega*t)` convention. They do not establish finite-band accuracy.

The current exact infinity acceleration term is zero at exterior free-surface
points. This does not imply every asymptotic contribution is zero: a conditional
`Hv~C/s`, `s=-i*omega`, corresponds to a displacement-to-elevation limit `C`.
The researcher is independently deriving that endpoint. Do not fit a derivative
term from finite-frequency growth or equate acceleration and displacement terms.

## Reusable stored inputs

| File | Contents / SHA256 |
|---|---|
| `artifacts/workflow/surface-response-attempt1/field/response.npz` |960 frequencies0.025..24rad/s, seven near/far/trace points, all six complex velocity transfers; six independent held-out frequencies; historical response arrays. SHA110380140cbab1dbcd2e6f72a6e5efc7ae0c0a6598af3a5fe62959065e7e8492. |
| `artifacts/workflow/refinement-checkpoint-attempt1/refinement/refined-transfer.npz` |486 frequencies, three far points, all six transfers and six held-outs. SHAec61104779972862575f9ceb7f926398436294e29eb47a1dd4148a4125e8a525. |
| `artifacts/workflow/surface-projection-attempt1/comparison.npz` |3001x7 historical signed-lag/positive/negative/cosine responses. SHA16e3e2670c3b7ce8255969c25856557c9dfdca2e747730cff0517d6af4157ac8. |

`scripts/pilot_surface_response.py` records the first cache construction and
held-outs `[.137,.687,1.337,2.637,5.137,9.937]rad/s`. The cache is historical
transfer evidence, not a newly run canonical trajectory. Reading ZIP/NPY headers
confirmed its array shapes without importing numerical libraries or computing.

## Ranked hypotheses and next discriminator

1. Finite bandwidth versus finite lag recovery: compare the same saved Heave
   transfer at near `[0,1.25]m` and far `[6,0]m`, cutoff12/24rad/s, lag60/100s.
   Longer Fourier recovery is distinct from the causal history needed for a60s
   transient. Keep windows below `pi/deltaomega`; no motion rerun is needed.
2. Missing finite-band tail/endpoint behavior: compare real and imaginary errors,
   not only magnitude, and derive the displacement endpoint independently. A
   correct endpoint is useful only with held-out complex-response validation.
3. Discrete near-surface Green cancellation: the default `low_freq` extraction
   combines exact Rankine integrals with numerical wave-part quadrature. An
   equivalent `high_freq` extraction is an existing Capytaine option and may
   discriminate cancellation/tail artifacts without changing physical geometry.
   This remains a hypothesis, not an established defect or selected correction.
4. Sign/factor/delay or cutoff effects: apply the identical production inversion
   to known causal `1/(1-i*omega)` and `s/(s+1)^2`, `s=-i*omega`, each undelayed
   and positive2s-delayed. The second has the relevant zeroDC/1s-tail structure
   and exact impulse `exp(-t)*(1-t)`. Independently integrate its two cascaded
   first-order ODE states and output `z1-z2`, driven by a smooth compact pulse.
   Finite-band negative-lag ringing is not evidence of physically acausal BEM data.

The equivalent Green expressions, asymptotic limits and default preference are
documented in [Capytaine3.0 theory](https://capytaine.org/v3.0.0/theory_manual/inf_depth_green_function.html).
Its integration note distinguishes numerical versus exact treatment of extracted
Rankine terms. That source motivates the discriminator, not a claim that switching
the option fixes this case. Broad mesh/irregular-frequency qualification remains32.

## Prepared first allocation and stopping point

Private helper: `artifacts/workflow/issue31-stored-diagnostic1/diagnose.py`.
Run after root claims and freezes the serial allocation:

```powershell
$env:PYTHONPATH = 'C:\Users\alex\topwater\src'
.venv/Scripts/python.exe artifacts/workflow/issue31-stored-diagnostic1/diagnose.py --input artifacts/workflow/surface-response-attempt1/field/response.npz --output artifacts/workflow/issue31-stored-diagnostic1/results
```

Freeze helper, cache, production surface module and imported package/dependency
identity. Proposed hard limits120s/2GiB/one numerical thread, internal115s budget.
Output is compact JSON; no images, model fitting, BEM or body integration. Each
point has its own transfer scale; phase is omitted when either response is below
1e-6 of its point's reference peak, and zero-peak relative metrics are undefined.
Reports preserve separate real/imaginary/amplitude/phase errors, analytic burst
rest/startup and finite-band raw negative-lag magnitude. No acceptance is inferred
from merely producing those diagnostics. Researcher/coordinator declare physical
replacement gates separately before candidate validation.

Stop after reporting this first discriminator and selecting ONE next cause test.
No production correction, WEC fork change, speed work, motion/forcing change,
full32-34 investigation or field regeneration is authorized by this note.

## Executed first diagnostic and next discriminator

Root executed the frozen helper under its serial allocation. The report
`artifacts/workflow/issue31-stored-diagnostic1/results/report.json`, SHA256
10c35890589b0e9ff43317e2b83da7aca2f019dc50614d228509afcb2ed705fb,
completed in4.594s (guard5.797s). Root reports terminal0, Jobempty, unchanged
frozen inputs, host process absence and release-last. No BEM/body integration,
source correction or published artifact change occurred.

Near Heave maximum complex error normalized by that point's held-out peak stayed
31.06--31.30% across cutoff12/24 and lag60/100; far stayed3.45--3.61%.
The near residual is almost entirely imaginary. Increasing lag60 to100s or the
tested band12 to24rad/s does not resolve it. At cutoff24 the undelayed analytic
endpoint control has maximum complex error2.69e-5, phase0.00315degrees and
smooth-burst peak error1.09e-6; delayed controls remain about2.2e-4 complex error
and1.5e-5 burst error. Those controls support the implemented sign/factor/delay
convention; they do not validate the BEM transfer or cover all modes.

The ONE next test uses existing native `low_freq` and `high_freq` extraction,
fixed Heave mode/facets/default quadrature and the same near/far points at
0.687,2.637,9.937,32,64rad/s. Cross-evaluate each identical source vector with
both field matrices, then separately resolve sources under each formulation.
Record bodyS/K, fieldS, source vectors, potential/elevation/Hv and their separate
real/imaginary differences. Also record the matched free-surface residual
`d_z psi - omega^2 psi/g`, with absolute values and scale equal to the larger
of the two individual terms; zero scales leave relative values undefined.
This distinguishes field evaluation from source-equation sensitivity without
a quadrature or mesh campaign. Equivalent splits may disagree discretely;
choosing whichever fits the cosine projection better is prohibited.

An infinity Heave solve records unit-velocity potential and its vertical
derivative. Normalize by `-i*result.omega` symbolically BEFORE converting to
numeric arrays. This supplies the independently derived conditional displacement
endpoint while retaining the zero acceleration feedthrough distinction.
Both extraction choices are natively overridden to high-frequency at infinity;
the report records that override explicitly.

Prepared helper: `artifacts/workflow/issue31-green-diagnostic1/diagnose.py`.
Inputs are the unchanged physical geometry JSON
`reference/capytaine/canonical-boat.json`, original17panel `native.mat` under
`artifacts/reference/minimum-20261024/native-gl2/`, the saved `response.npz`
above and its existing `field/green-cache` directory. The helper copies that
cache into its own output before use, so the historical cache remains immutable.
Prospective cap120s/2GiB/one numerical thread; no launch before root allocation
and independent helper review. No new production method is selected.

The first two Green diagnostic launches failed in diagnostic instrumentation
before publishing results: a post-solve engine cache returned an LU wrapper
instead of numericK, then an incorrect pre-v3 solver-level Green attribute was
used. Root verified terminal/frozen/process cleanup and lease release before
each revision; both failed attempts and exact helper versions are preserved.
Current extraction uses `solver.engine.green_function.evaluate` for raw body
matrices, with numeric type/shape checks. These were helper API errors, not
physical failures or evidence about the radiation method.

The installed symbolic array stack handler also drops the requested axis in
the public infinity-velocity route. The diagnostic avoids that route by
normalizing the symbolic source vector FIRST and evaluating numeric fieldS/fullK
matrix products with explicit component stacking. At finite frequencies it
cross-checks the same numeric-source path against the public potential/velocity
evaluators. No library or production correction is included; the endpoint is
ancillary to the unchanged finite-frequency discriminator.

## Successful Green discriminator and bounded follow-up

The third Green attempt completed in0.5s (guard3.329s); root verified terminal0,
unchanged inputs/process cleanup/release-last. Report
`issue31-green-diagnostic1/results3/report.json` has SHA256
9da95c7d1678841d8396f835cb613c42b9aa8160323fc53f603e8221ba1bc577;
independent diagnostic acceptance is `issue31-green-review1/acceptance.json`
(a28c7cf7). The two earlier instrumentation failures remain preserved.

At32rad/s, using the SAME low-split sources and changing only the field split
changes near Heave imaginaryHv by-0.287242, versus low matched+0.289287.
At64rad/s the corresponding change is-0.624741 versus low+0.620008. This exposes
large field-evaluation cancellation sensitivity before any source re-solve.
Source-equation differences are also substantial, especially at lower frequencies.
The low split's scaled near free-surface-condition residual is about99--101% at
32/64rad/s. The high split is also unqualified (about162%/35% there) and much
more ill-conditioned at0.687rad/s. Do not select high split as a blanket fix.
The exact infinity unit-velocity potential is zero while its vertical derivative
is[-0.8376657,-0.01604976] at near/far; this distinguishes displacement from
acceleration endpoint behavior without validating a finite-band tail correction.

The one next discriminator is native Gauss-Legendre2 versus centroid integration
on the SAME17facets, both splits/five frequencies/two points. It records the
Green-matrix free-surface residual before source contraction, source-solution
residuals and cross-evaluation. A potential-only one-sided derivative witness at
2.637rad/s near point uses frozenh=1e-4 and5e-5m to independently check the
public gradient path. This is a localized integration diagnosis, not a general
mesh campaign or permission to alter accepted dynamics. Prepared private helper:
`artifacts/workflow/issue31-quadrature-diagnostic1/diagnose.py`; root allocation
and exact-source independent review are required before execution.

The quadrature diagnostic completed in0.672s (guard2.547s), with root's
terminal/frozen/process/release checks. Report SHA256
2b9e3ce4aa74b68bcbe244478c8dfb506ebed11a1cdb7378105fe3907bf82075 is in
`issue31-quadrature-diagnostic1/results/report.json`; independent receipt
`issue31-quadrature-review1/acceptance.json` begins85f85a18.
GL2 reduces the low-split raw residual by roughly2.7x but does not qualify it;
the near64rad/s contracted residual remains95%. High-split raw residual is
frequency-independent,0.0372878 for centroid versus0.0152015 forGL2. The
potential-only derivative witness does not justify selecting the high split:
its GL2 difference from the public gradient is about100% at the tested near
point, while centroid-low is0.0184%. These are diagnostics, not surface-error
percentages.

Source algebra predicts the low-split free-surface residual
`2*k*(Q_R-I_R)` and high-split residual `2*(d_z I_R-d_z Q_R)` atz0, whereI is
the native analytical Rankine integral andQ is its panel-quadrature sum, both
including-1/(4pi). The observed frequency scaling matches those identities.
The next kernel-only helper fixes original panel0 and nearpoint, with three
frequencies; it verifies these identities without a BEM solve. At2.637rad/s
only, native no-tabulation theta1001/2001 evaluation and the declared derivative
h-pair separate integration mismatch from table interpolation. No-tabulation
uses safe one-entry radial/vertical arrays; zero-size grids are avoided because
release Fortran writes radial entry1 unconditionally. Prepared helper:
`artifacts/workflow/issue31-kernel-diagnostic1/diagnose.py`; independent review
and root allocation precede any run.

The release [matrix assembly](https://github.com/capytaine/capytaine/blob/v3.0.0/libs/libDelhommeau/src/matrices.f90)
and [wave evaluator](https://github.com/capytaine/capytaine/blob/v3.0.0/libs/libDelhommeau/src/Green_wave.f90)
show why the third option is incomplete as a repair: it moves one vertical
derivative Rankine contribution from numerical quadrature to analytical
integration, while leaving the potential unchanged. A correction must treat
potential and derivative together. For regular exterior field points, summing
the full elementary Green function and gradient with the same quadrature is a
concrete consistent representation to evaluate after diagnosis. No native
all-terms-quadrature switch was found; native mesh quadrature changes the wave
part only. Such a candidate would still need the independently frozen all-six
complex-response and causal-transient gates. It is not selected by this note,
and its spatial approximation remains distinct from full issue32 qualification.
