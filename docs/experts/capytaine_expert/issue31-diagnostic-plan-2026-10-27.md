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

## Verified cause and first candidate, frozen before validation

Kernel report `issue31-kernel-diagnostic1/results/report.json`, SHA256
185a1cbd2f864ced446d439d938f154617be84ebc87c2ab4608d73702f03e287,
completed in0.265s (guard1.891s). Root reports terminal0/frozen/Jobempty and
process cleanup. Predicted low/high extracted-Rankine residual identities match
within1.3e-7/1.9e-9 absolute on the fixed panel. At2.637rad/s, removing table
interpolation reduces the split-potential mismatch from0.00350 (centroid) or
0.000231 (GL2) to roughly1e-7. Potential-only derivatives also improve, with
finite angular-integration and finite-difference floors retained. These findings
identify mixed integration and separate interpolation errors; they do not
qualify the original17panel hydrodynamics.

Root authorized the first PRIVATE field candidate. Keep original low-split,
centroid, tabulated bodyK/source equations unchanged. At strictly exterior
regular field points, form the SAME full Green integral and full gradient with
the existing centroid points/weights:

```text
S_common = S_native - I_R - sign*I_ref + Q_R + sign*Q_ref
grad_common = grad_native - grad_I_R - sign*grad_I_ref
              + grad_Q_R + sign*grad_Q_ref
```

The reflected analytical integral is obtained at the reflected FIELD point;
its vertical derivative changes sign. The independent point-quadrature kernel
uses-1/(4pi*r), with gradient with respect to the field coordinate. Both potential
and gradient are changed together. This is a consistent regular-field quadrature
representation, not a fitted correction, motion change or blanket accuracy claim.
Original native transfers must replay the immutable saved data before any
candidate conclusion. Reject interior/boundary field points explicitly.

Private helper `artifacts/workflow/issue31-field-candidate1/diagnose.py` uses
Heave, the existing near/far points, training0.025..24rad/s and the six original
held-outs. Training uses native directtheta1001 evaluation; held-out reference
uses the equivalent high split/directtheta2001 with the SAME original sources.
Low-split directtheta2001 controls separate angular-integration sensitivity;
tabulated common-integral values and original native values remain separate.
The direct reference is independent of table interpolation, not continuum truth.
Production cosine inversion and lag60 remain unchanged; lag100 is diagnostic
only and cannot waive a failed production gate.

Root adopted these prospective gates before candidate results:

- EACH point hasS=max(abs(held-out reference transfer)); every complex error
  <=0.01S. Where abs(reference)>=0.05S, amplitude error<=2% and phase<=1degree.
  Exactly zeroS uses1e-12s absolute Heave-transfer allowance and undefined phase.
- Native directtheta1001 versus2001 held-out transfer difference<=0.001S.
- Actual LOW1001 candidate versus actual HIGH2001 reference AND LOW2001 versus
  HIGH2001 split consistency each<=0.001S. Forward candidate field transfer must
  also pass complex/amplitude/phase gates independently of kernel recovery.
- Regular nonzero modal Robin residual<=1e-4 of the larger of its two terms;
  absolute errors are retained, exact zero scales remain undefined relatively
  with1e-12 dimensional roundoff allowance.
- At2.637rad/s near point, BOTH fixed one-sided derivative errors at
  h=[5e-5,2.5e-5]m and their pair difference must be<=1e-4 times the declared
  same-point nonzero derivative scale. No best-h selection or retrospective
  tolerance adjustment. Zero scale uses1e-12 absolute allowance.

Fail closed unless the report contains exactly12 required Robin blocks and2
fixed-pair derivative blocks. A failed field/reference/angular control is not
evidence of body-source acausality. The original pre-review c535 helper is
preserved beside the revised helper; these gate-chain additions were adopted
before any candidate execution or results.

Proposed allocation300s/2GiB/one numerical thread, internal285s budget, with
independent exact-source/API review before root freezes and runs it. If the
original source operator prevents causal recovery, stop and report the block;
do not tune kernel/window settings or fit away the discrepancy. All-six modes,
causal burst tests and full saved-motion/publication acceptance remain required
after an informative first Heave candidate. No production implementation has
changed at this authoring checkpoint.

## First candidate outcome and staged elementary-point follow-up

The first candidate completed in14.796s (guard16.484s) with root's terminal0,
frozen/process/Jobempty checks. Report
`issue31-field-candidate1/results/report.json` SHA256
d4d1e1ea6111991c72a82dce646f93c0d87e68c24923c8607dc09f5d0f10ea90
has a FALSE verdict. Robin, broad forward-field and theta gates passed; strict
actual-reference/split and derivative controls failed. Far reference disagreement
reaches0.00021543s versus0.00015330s budget. Low1001 derivative pair fails while
low2001 passes; HIGH2001 reference derivative disagreement is about13.4%.
The remaining near cosine-recovery residual is0.070217s (22.7% of its transfer
scale), much larger than the near reference-value difference, but failed
reference controls prohibit attributing it solely to source causality.

The next PRIVATE helper `issue31-point-candidate1/diagnose.py` assembles the
native elementary wave kernel with direct/reflected POINT Rankine values and
their paired gradients, identical centroid weights and-1/(4pi) normalization.
It avoids analytical PANEL subtraction. Installed scalar/vectorized f2py API
signatures were verified from binary metadata and release `interface.f90`.
The original bodyK/source equations, geometry and cosine kernel remain unchanged.
Native `exp_e1` still has an internal single-precision special-function step;
point assembly does not promise to remove that inherited precision floor.

Stage1 evaluates ONLY six held-outs and the unchanged field/reference/theta/
Robin/fixed-pair derivative gates. Any failure publishes its negative evidence
and stops BEFORE dense960frequency training. No field failure is attributed to
source causality. A passing stage permits the same dense Heave near/far recovery
check. No new fitting, band/window selection, SciPy formula, production field
correction or physics output change is included in this helper.

Any eventually accepted shared field-operator correction must be checked for
both radiation AND diffraction evaluation, preserving incident components and
accepted boat motion/forcing. This does not authorize a separate issue33 load-
closure campaign. Current execution remains the single staged Heave candidate.

The native elementary-point attempt completed in0.609s (guard2.407s), with
root's successful cleanup/release. Report86ff5da6 under
`issue31-point-candidate1/results/report.json` has FALSE field controls and
correctly skipped dense training. Reference/derivative floors are essentially
unchanged after removing analytical panel subtraction. No production correction
or source-causality conclusion follows.

## Focused complex128 special-function witness

Private helper `issue31-exp1-witness1/diagnose.py` evaluates the SAME published
elementary Green integral at original panel0 centroid/nearpoint/2.637rad/s only.
It uses existing `scipy.special.exp1`, not a new E1 implementation or library fork.
Forzeta=z+i*r*cos(theta), theta in[-pi/2,pi/2], the imaginary coordinate is
explicitly nonnegative, including+0 at the negative-real branch endpoints.
The integrand is exp(zeta)*E1(zeta)+i*pi*exp(zeta). Published gradient identities
are differentiated with respect to the field coordinate. The exact Lipschitz
identity supplies I2=I1+2/r1; a separately quadratured I2 remains diagnostic.

Both low/high full-point potential and gradient are retained with-1/(4pi)
normalization. Independent outgoing-wave checks compare numerically integrated
imaginary components with BesselJ0/J1 identities. Native1001/2001 and SciPy
complex1281001/2001 cases retain the fixed derivativeh-pair, Robin and evaluator
budgets. The one-point witness does not qualify all source panels, frequencies,
Heave transfer, all-six reconstruction or production physics.

References: [Capytaine3.0 Green theory](https://capytaine.org/v3.0.0/theory_manual/inf_depth_green_function.html)
and [SciPy E1 branch convention](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp1.html).
Exact-source independent review and coordinator allocation precede execution.

The witness completed in0.094s (guard1.782s); root verified terminal0/frozen and
cleanup. Report SHA2560d6b3793ccdce207a095564a28eeda293c429a8a3ee4bd3d261cde42398d3093
is in `issue31-exp1-witness1/results/report.json`. Independent acceptance
`issue31-exp1-review1/acceptance.json` SHA52ddd89f1eaa06e62ad57c2bebdd79d1a076c43285ddc81d170086ab2537260f
accepts ONLY this fixed-point precision witness. SciPy FD errors are about
2e-11--9e-11, native about5e-6--2.6e-5; SciPy Robin error is about2.8e-17 and
independent Bessel errors are at floating-point roundoff. Attribution is the
native precision/approximation path, not uniquely E1: native angular-grid and
normalization arithmetic also contribute.

## Prospective adaptive full-Heave field candidate

Root confirmed adaptive angular integration BEFORE this candidate's results,
because near-surface source centroids create angular endpoint widths around
1e-4rad, below a fixed1001point grid step. The physical integral, original
sourceK/modal basis and cosine kernel remain unchanged. Existing complex,
amplitude/phase, Robin and fixedh-pair gates are unchanged.

Private helper `issue31-scipy-candidate1/diagnose.py` uses witnessed complex128
SciPy exp1 expressions with quad_vec, normmax, workers1 and subdivisionlimit2000:

| Role | Rule | epsabs | epsrel |
|---|---|---|---|
| Candidate | gk21 |1e-11|1e-10|
| Independent held-out reference | gk15 |1e-12|1e-11|
| Reference convergence control | gk21 |5e-13|5e-12|

Geometry-derived breakpoints resolve the narrowest angular endpoint layer;
they are not tuned to response phase. Every integration records error estimate,
status, evaluation count and interval count, and fails on unsuccessful convergence
or subdivision-limit exhaustion. Propagated angular transfer estimates and
candidate/reference/reference-convergence differences must satisfy the same
0.001S per-point budget. Independent Bessel imaginary checks remain active.
The angular-convergence name replaces theta-convergence because integration is
now adaptive; no transfer or derivative tolerance is relaxed.

Six held-outs and all evaluator/reference controls precede dense960frequency
training. On any failure, stop before training. If field controls pass, retain
the same60s production recovery and100s diagnostic window, with no best-window
waiver. This remains private Heave/two-point work; all-six, startup, complete
saved-motion and common radiation/diffraction qualification remain ahead.

Before execution, root required every evaluator call (including derivative,
reference and dense-training calls) to pass its existing independent Bessel
control. Every returned angular estimate must be finite/nonnegative and no
larger than max(epsabs,epsrel*maxabs(integral)), with finite integral values;
solver success alone is insufficient because the returned estimate includes
roundoff. No multiplier or tolerance waiver is applied. Original pre-review
a90 helper bytes remain preserved. Reported estimates remain numerical estimates,
not rigorous physical error bounds.

The first adaptive launch stopped after3.672s at the self-imposed tighter
angular request9.937rad/s: status2 roundoff, estimated error5.39675e-13 versus
requested3.21566e-13. This negative attempt and exact a3 helper remain preserved;
root verified terminal99/frozen/process cleanup/release-last before revision.
The reference requests above were then prospectively revised, with independent
review agreement: candidate remains1e-11/1e-10, reference1e-12/1e-11 and tighter
control5e-13/5e-12. Only these primitive quadrature requests changed. All physical
gates, status/error-budget checks, Bessel controls and derivative tolerances remain.

Saved original Heave source vectors give9.937rad/s sum(area*abs(sigma_v))=4.12857.
The corresponding angular-error-to-transfer factor is6.69959s, so the observed
roundoff floor propagates to about3.62e-12s, versus the unchanged far-point
0.001S budget about1.53e-4s. The revised primitive requests remain far below
that budget while permitting SciPy's internal tol/8 convergence test above the
observed floor. No failed integration is accepted or silently waived. Root owns
the new results2 allocation; no rerun was performed during source revision.

The second adaptive run completed in158.562s (guard160.422s). Its results2 report
SHA98036466d170b33dd2418f5cde019928bafaabdc685fe76f8276a9711c6301e1 passes ALL
field/reference/angular/Robin/derivative/Bessel controls, while cosine recovery
still differs by0.0702219s near (22.66% of reference peak) and about1.66% far.
The100s lag diagnostic is nearly unchanged. This isolates a reconstruction/data
question, but does not establish a bodyK defect: the shallowest source centroid
depth0.000922m implies a Gaussian frequency scale around103rad/s, above band24.

Root's next allocation is ONLY a cheap known-causal point-Green bandwidth
control. `issue31-green-bandwidth1/diagnose.py` imports the extracted qualified
primitive `issue31-point-primitive-v1/green.py`; its normalized mathematical
diff is only a shape-derived centroid count. Original404 and all prior outputs
remain untouched. Actual source centroids0/1, fixed nearpoint, input volume flux
and transfer units s/m^2 are explicit. Analytic Bessel real transfer drives the
same24rad/s cosine and60/100s recovery. Independent finite-band KK24 uses Cauchy
PV quadrature. Compare cosine-minus-KK24 separately from KK24-minus-full-Green;
no high-band scan or body/source correction is allocated.

## User-directed end quantity: saved-trajectory surface height

Latest steering asks how to establish the surface answer, analogous to motion
accuracy, for visualization/diagnostics and a future native-solver benchmark.
Postprocessing cannot change the accepted WEC-Sim motion, but it can change
surface/contact/wetted-set diagnostics. Internal transfer and derivative tests
are diagnostics, not a substitute for a height requirement.

Root's1cm maximum surface-height target is a PROPOSAL, not user acceptance.
The proposed comparison freezes the saved trajectory and declared displayed
domain/times, comparing TOTAL elevation at identical valid exterior coordinates
against an independently qualified causal fixed-mean linear reference. Report
maximum/RMS difference and reference numerical uncertainty separately; an adopted
budget should cover both. Endpoints and near-waterline points remain included.
A representative-point screen can reject cheaply but cannot prove a full-domain
maximum. This would qualify reconstruction of the chosen linear baseline, not
true-fluid accuracy or the full deferred mesh/applicability investigations.
No new bodyK correction, expensive campaign or production change is currently
allocated. Finish the authorized cheap bandwidth-control review, then stop at
the height-reference proposal for coordinator/user prioritization.


## Accepted height requirement and actual-motion evidence ? 2026-10-27

The user subsequently CONFIRMED <=10mm maximum numerical surface-height error
INCLUDING reference uncertainty against the chosen LINEARIZED model. Earlier
wording above describing1cm as a proposal is historical. No nonlinear truth
comparison is requested. Keep accepted motion/forcing and physical inputs fixed.

The completed ALL6 seven-probe screen recomputes current cosine radiation from
NEW accepted native velocity/acceleration; it never reuses old radiation history.
Cosine24 minus full-complex reference-candidate reaches42.6603mm at44.68s on a
continuation trace; the direct near exterior reaches31.9271mm. Its uncertainty
sensitivity sum reaches2.91712mm, mostly interior future-cut controls. These are
method differences with an unqualified reference, not certified physical errors.
FFT sign, exact periods/time origins, startup, smooth future continuations,
interpolation and original motion-knot preservation were checked. Response above
24rad/s was explicitly unknown. Root verified terminal/frozen/process cleanup.

The completed same17-facet native ALL6/7-probe24-to48 extension took26.422s.
Report `artifacts/workflow/issue31-tail-extension1/results/report.json` SHA256
7de711731e238e3bb65a9dd15854ae02a39364ec4d4759cfacccc355053a75a2.
Original0..24 samples remain unchanged; native20/22/24 overlap error is exactly0.
Maximum saved-motion48-minus24 change is0.0591576mm. The36-to48 absolute
input-weighted spectral envelope is at most0.0386197mm and bounds ONLY that
represented finite band, not response above48. Updated observed sensitivity sum
is at most2.86777mm. The selected-case tail evidence is encouraging but does not
certify an infinite-tail bound or whole-domain accuracy. Reference remains
unqualified; whole5929-by751 TOTAL surface acceptance remains false.

This job used a300s/4GiB/one-thread guard, internal240s cap; actual child26.422s,
terminal0/guard0/Jobempty/frozen, root independently verified host process absence,
idle readback and release-last (finish36942). No body integration, motion, sourceK,
mesh or forcing change occurred. Existing scientific failures remain preserved.

Next bounded source-only work is a private causal signed-lag reconstruction
comparison using the SAME saved ALL6 response and accepted trajectory. Derive
and independently review the finite-low native operator's instantaneous
acceleration endpoint before inversion; do not fit it or substitute the symbolic
infinity branch. Restore saved acceleration exactly once. Compare reference and
candidate using the SAME regular memory/direct term and report differences from
the previous spectral-derivative FFT reference explicitly. Smooth-future
sensitivity is not a bound on the actual negative-lag contribution. No production
method is selected until measured actual-height evidence and review support it.


## Causal-cut rejection and bounded implementation ? 2026-10-27

Causal attempt3 completed in1.078s with numerical sign/delay/direct controls
passing. Report `artifacts/workflow/issue31-causal-attempt3/results/report.json`
SHA0e773e274430264d4a0a868a8b88271f905592bd95723d9ab8a9c27b991d5489.
The analytically derived finite-low discrete acceleration endpoint has condition
4920.88 and source residual4.44e-16. Removing that term before signed inversion,
cutting negative lags and restoring savedA*a changes actual height by194.933mm,
versus42.593mm for the original cosine policy. Its36-to48 sensitivity is11.125mm.
Reject this cut; do not tune band, endpoint or model to rescue it. The full-H FFT
reference and regular-memory-plus-savedA*a reference differ by only0.286675mm.
Sparse regular response grows at96rad/s; no universal no-growth/tail proof exists.
The two instrumentation failures (missing tabulation directory and symbolic DOF
array use) and their exact helper/plan bytes remain in attempt1/2-source folders.

Root authorized implementation of intact full-complex offline FFT(v), with no
runtime endpoint extraction, causal cutting, physical motion/model change or
separate direct acceleration addition. Scope is `_capytaine_surface.py` and
focused field tests. The generic API preserves raw startup, saved velocity knots
and full complex phase; fixes zero velocity prehistory, an explicit smooth final
v/a-to-rest continuation and a declared period/grid/frequency band. Spatial chunks
reuse the native source bank and bound FFT/transfer storage. Native incident,
diffraction, exterior traces and geometry/startup gates remain unchanged.
One active readiness flag is FALSE pending whole-domain qualification. A private
qualification process may explicitly override it with recorded scope/source;
there is no public bypass. No tests or numerical jobs were run while authoring.

Whole-domain qualification plan: freeze accepted NEW motion, native operator,
sea inputs, point/trace mapping and original5929-point/751-frame clocks. Replay
all6 radiation inputs in spatial chunks with full H, comparing production results
against independently assembled full-complex Fourier output at IDENTICAL points.
Retain exact same incident/diffraction histories, and verify the TOTAL difference
and unchanged motion/forcing explicitly. Run declared period512/1024, dt.02/.01,
future10/20, interior futurecuts20/40 and band36/48 controls over ALL valid points,
plus startup/pre-response and savedA*a first-sample diagnostic. Reuse transfer
chunks across controls, not repeated BEM. Sum declared grouped sensitivities,
implementation differences, direct/interpolation semantics change where relevant
and a separately justified remaining-tail allowance. Require max error plus
reference uncertainty <=10mm; seven-point evidence alone cannot pass this gate.
If reference uncertainty or current physical startup gates fail, report the
specific conflict and stop. Do not silently zero initial radiation, alter the
trajectory, fit phases or expand into body/sourceK/mesh corrections. Tail policy
is an explicit engineering assumption awaiting its separate bounded review.


## Reviewed gated implementation and measured-pilot proposal

Independent source review accepted the new production full-H API, chunking,
metadata and focused regressions. Root's focused test run passed42 tests in17.81s;
prior sandbox temporary-directory setup failures are preserved. Root verified
terminal/frozen/process absence/idle/release and committed d3630fe. The active
production readiness gate remains FALSE; tests do not certify physical accuracy.

The separate input-tail screen completed in0.188s. Root reports max ASSUMED
full-H tail plus known uncertainty4.457671mm at the7 existing probes. It remains
conditional on its growth envelope and beyond-Nyquist treatment; it does not
qualify5929-point coverage or justify copying a seven-probe budget everywhere.

Fresh private source: `artifacts/workflow/issue31-whole-surface1/pilot.py` with
plan.md alongside. Exact current native geometry/context identity is asserted,
128 actual unique exterior queries combine extent and nearest-source coverage,
and all unchanged0.025..48rad/s6DOF sources are computed once and retained. Then
actual field evaluation reuses that bank with an unchanged RHS count. Production
FFT and four control costs are timed; independent inverse Fourier sums check
actual field spots. Root guard OS peak working-set/commit measurement MUST
accompany the pilot before a larger allocation. Source-only authoring/AST checks
have run; no pilot or new scientific calculation has been launched by this expert.
The plan specifies pointwise whole-grid uncertainty, verified common incident/
diffraction, startup/direct-term checks, existing geometry gates and <=10mm
combined criterion. Stop after pilot review; no speculative large allocation.


Before any pilot data, root required transfer-frequency interpolation uncertainty
separately from high-cutoff convergence. The pilot and full recipe now compare
native.025 spacing with coarsened.05/.1 from the SAME stored H, with no new solves.
Whole pointwise uncertainty must include this grouped sensitivity and observed
convergence; only a measured-height-error justification can propose targeted
midpoint checks. No automatic denser BEM campaign is authorized. Existing guard
records continuous kernelJobcommit peak, sampled aggregate working-set/private
and per-process lifetime memory peaks; pilot PID binds those external measurements.
No new native memory-wrapper implementation was needed. Pilot source/plan remain
pending independent final-byte acceptance, not numerical execution or qualification.


## Actual pilot and full-runner source checkpoint

Root's pilot completed62.641s, with source warmup50.484s, native cached field
128queries5.813s, productionFFT1.093s. Actual domain5899unique exteriorqueries,
5929targets,384continued nodes.11520RHS sourcebank is retained in pilot-results.
Guard peakJobcommit543100928B and sampled aggregate working set491053056B with
no gaps. Direct inverse-sum error8.67e-16m. Source/context identity passed. Pilot
spacing.05/.1 maxima0.140293/0.707945mm, finalfuture10max0.956670mm; no pilot
accuracy rejection. Root verified terminal/frozen/CIM/idle/release. Full-domain
cost/uncertainty still needs the actual full runner; no sourcebank warmup repeats.

Root explicitly allocated only64/96 endpoint-envelope calibration (12RHS plus
2field frequencies/chunk), no band sweep. Full source-ready runner is qualify.py
and qualification-plan.md alongside the preserved pilot. Prospective root/NMR
cap1200s/3GiB/chunk128/internal1140, justified by measured pilot costs. It measures
and checkpoints actual firstchunk controls/oldcosine work and rejects a failed
combined budget or necessary resource projection early. Conservative budget
failure does not establish actual physical error>10mm.

The revised VALIDATION reference is fullcomplexR=H+i*w*A on C1-Hermite velocity
through exact savedv/a, plus exact saved/right-limitA*a at original output knots.
Runtime is unchanged fullH linearv FFT with rawstartup; no A is added there.
At0 symmetric spectral derivative and saved right limit differ ifa0!=0; measured
candidate-reference difference includes this convention. No startup overwrite.
The independently derived compact-cubic jump formula bounds all omitted regular
continuum above48 (includingaboveNyquist) by C/pi times
abs(a0)/48+J2/(2*48^2)+J3/(3*48^3), summed over modes. C is pointwise2maxknown|R|
fromdense36..48 plus64/96; unseen-growth/no-resonance remains an explicit ROOT
engineering assumption, not measured proof. The new hermite_tail.py and tiny
check_hermite_tail.py controls use pulse and multi-interval polynomial quadrature,
zero/scaling and exact endpoint-jump tests; source review and root execution are
separate from the full wrapper. No scientific imports/run occurred while authoring.


## Full attempt1 conservative rejection and variant2 source

`issue31-full-attempt1/results/report.json` records an EARLY REJECTION after128
queries: maxcombined12.323328mm, tail11.715798mm, future0.298435mm, spacing0.600117mm,
interpolation0.006531mm and dt0.0000472mm. NMR measured actual candidate/reference
max0.009218mm in the retained checkpoint. The failure is the loose absolute-jump
uncertainty bound, not demonstrated12mm physical error. Root proved terminal0,
frozen/CIM/idle/release; all original helper/runner/plans/rootpolicy1 bytes were
preserved under issue31-whole-surface2/preserved-attempt1 with hashreceipt, and
original paths are untouched. Firstchunk took18.485s with140.436s preloop setup.

Root/NMR prospectively selected a tighter mathematical INPUT-only enclosure,
not a solver/kernel/physics or criterion adjustment. Fixed coefficient-transform
period8192s and targetupper512rad/s; runtime responseFFT remains1024s. Cover
floor(48/h)h..ceil(512/h)h. Exact jumpDFT numerators are periodic/conjugate folded,
but denominator powers use ACTUALfrequency. Magnitude trapezoid bounds complex
linear interpolation by convexity; explicit perinterval M*h^3/12 covers its
second-derivative remainder, an explicit phase allowance covers rounded input
coefficient times, then the original jump bound covers continuum beyondactualB.
Multiply by SAMEpointwiseC/pi once. No term is simply dropped or fitted.

New a3ec Hermite helper/control sources and tail-enclosure-plan are independently
reviewed before root's tiny allocation. Optional cached128 diagnostic freezes
acceptedmotion and retainedfirstchunk, verifies oldcombined from unchanged
candidate/reference/allU, replaces ONLY oldtail, retains10mm/C/otherU and hashes/
arrays; fullscope qualification stays false. It avoids another field/setup job
if the genuine tighter bound still cannot qualify the known chunk. No author
scientific imports or numerical execution occurred while source was drafted.
New fullvariant01f5 also instruments exact preloop stages (policy/input/cache/
geometry/Green/bank/trace/reference/endpoint/prior7/train/incident spatial+temporal).
No broad performance changes, BEM sweep, WEC rerun or readiness flip is included.


## Final selected-case acceptance and canonical delivery, 2026-10-27

Authoritative final method review accepts the full5929-target/751-time surface
(5899 actual exterior queries,384 continued targets):0.341183mm maximum runtime/
reference difference and5.784425mm combined difference plus reference uncertainty,
within10mm. All751 unchanged physical geometry/contact/opening/startup gates pass.
Accepted motion/forcing, chosen fixed-mean operator, incident/diffraction/radiation,
clocks/basis and trace continuation are preserved. No WEC/body rerun or physical
source correction was needed. Evidence: issue31-final-method-review1/acceptance.json
and issue31-independent-check-attempt1/results/report.json under artifacts/workflow.

The conditional pointwise factor2 regular-response growth/no-unseen-resonance
assumption remains explicit; numerical reconstruction acceptance is NOT nonlinear,
experimental or mesh/model-resolution qualification. Issues32-34 remain open.
The standard canonicalNPZ/PNG and truthful cached-surface result.log are published
with prior artifacts backed up; result SHA prefixd3438810 and PNG275ac.... Six
unchanged-motion calm-restoring warnings and the new18.27% peak wet-set warning
retain their provenance; retired42.7mm projection warning is excluded. Video was
regenerated:751 frames,12.5fps,800x800, source0..60s,142.641s render, no previews
inspected. Publication/recovery reports are preserved under their attempt folders.

Root promoted ONLY the readiness Boolean, source SHA prefix6d8edb3; numerical
function AST is unchanged. Final integration regressions pass50 tests in24.20s,
with no skips and independent receipt review. Repository/issue bookkeeping
remains coordinator-owned. Recovery geometry
cost540.829s includes wet-set integration355.109s and disturbance metrics157.192s;
these measured costs guide the next speed phase without authorizing optimization
here. This final documentation update ran no workflows and changed no source.
