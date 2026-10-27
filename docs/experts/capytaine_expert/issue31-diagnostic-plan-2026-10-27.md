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
