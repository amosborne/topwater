# Issue 31: independent method proposal and prospective gates

2026-10-27. Bounded source/literature review for the newly authorized goal
`resolve #31`. This note proposes diagnostics and acceptance criteria; the
coordinator selects execution and adopts gates before validating a candidate.
No numerical imports, calculations, BEM/body runs or production changes were
performed. The reviewer must not author the physical correction being reviewed.

## Scope and verified facts

Only radiation-surface reconstruction in `topwater solve --solver=wecsim` is in
scope. Saved six-DOF motion/forcing, incident train and diffraction remain fixed.
Speed work and the full investigations in issues32-34 are later work; a small
discretization discriminator can be necessary to identify this issue's cause.

The installed distribution metadata identifies Capytaine3.0.0. Current
`src/topwater/_capytaine_surface.py` (SHA256
`240198ee21b92334b2b195adc4e7d5a31ff347f9c0320a95d16569e832520251`)
implements displacement-to-velocity
conversion `Hv=i*eta_q/omega`, positive-sign transfer recovery, real-part cosine
inversion and trapezoidal zero-prehistory convolution. It is the current source,
not the historical source hash recorded by the original investigation.
Prior evidence/definitions are in
[the original phase investigation](causal-phase-issue31-investigation-2026-10-25.md)
and [provisional method review](surface-reconstruction-review-2026-10-25.md).

No new cause is established here. The known17-30.4% recovery residuals have a
per-mode maximum denominator over sampled points/frequencies; they are not
pointwise percentages, total-surface errors or trajectory errors.

## First-principles endpoint audit

Capytaine's [theory manual](https://capytaine.org/stable/theory_manual/theory.html#free-surface-elevation)
uses `exp(-i*omega*t)`, displacement-normalized body radiation and
`eta=i*omega*Phi/g`. Its linear surface condition is
`g*dPhi/dz=omega^2*Phi`. The following algebra applies separately to each mode
and regular exterior observation point; rotational and translational units
must remain separate.

Let `psi` be potential per unit modal velocity. Then

```text
Phi_q = -i*omega*psi
Hq = eta/q = omega^2*psi/g
Hv = eta/v = i*Hq/omega = i*omega*psi/g
Ha = eta/a = -Hq/omega^2 = -psi/g
```

Equivalently, with causal Laplace variable `s` and `s=-i*omega` on the frequency
axis, `Hv=-s*psi/g`, `Ha=-psi/g`, and `s^2*psi+g*dpsi/dz=0` at the surface.
These are algebraic transformations of the same operator with matching initial
conditions, not different physics or a reason to change the body dynamics.

**Low frequency:** if the unit-velocity potential has a finite static exterior
limit, `Hv(0)=0` and `Ha(0)=-psi0/g` is finite. A generic `1/omega` pole in `Hv`
does not follow from the division in its conversion. Verify this limit rather
than fitting a pole or merely declaring the supplied zero endpoint correct.

**High frequency:** the documented
[Green-function limit](https://capytaine.org/stable/theory_manual/inf_depth_green_function.html#singularities-extraction)
is Rankine minus reflected Rankine. At a regular exterior `z=0` point these
cancel, consistent with zero acceleration feedthrough. A zero value of the
potential does not imply a zero normal derivative. IF the exterior derivative
has finite limit `C=dpsi_inf/dz`, the surface condition gives

```text
Hv ~ C/s,  Ha ~ C/s^2,  Hq -> C.
```

Thus a nonzero instantaneous displacement term and finite initial velocity-kernel
value can coexist with zero instantaneous velocity/acceleration feedthrough.
The bounded-derivative and convergence assumptions require checking on the
actual discrete representation. The infinity source normalization is symbolic;
do not numerically divide a literal infinity. Capytaine documents its
[endpoint evaluation conventions](https://capytaine.org/stable/user_manual/post_pro.html#pressure-velocity-free-surface-elevation).
A derived endpoint can be compared against fixed-point/mode trends, never fitted
from whichever finite-frequency value happens to be largest.

For a regular real causal kernel, `H(omega)=integral h(t)*exp(i*omega*t) dt`.
Real-part cosine inversion is exact only with appropriate whole-band data,
endpoint terms and sufficient lag coverage. Finite-band inversion instead obeys
`h_cos(t)=h_raw(t)+h_raw(-t)` at positive lags: it folds finite-band negative-lag
ringing forward. Negative lags alone do not establish noncausal physics.

## Cheapest discriminator before a replacement

Use saved Heave transfers at near `[0,1.25]` and far `[6,0]`, unchanged input
identity, existing independently solved held-out frequencies, cutoffs12/24 and
lag windows60/100s. Keep100 below `pi/domega` for the0.025rad/s grid; do not
confuse finer lag sampling with additional frequency information. Split real,
imaginary, amplitude and phase residuals by point and frequency, including the
low-frequency/DC behavior. No new BEM or body run is required for this stage.

Run the SAME diagnostic on known causal controls. In addition to a lowpass and
its delayed version, include `G(s)=s/(s+a)^2`, `a>0`, whose exact causal impulse
is `exp(-a*t)*(1-a*t)`. It has zero DC, a `1/s` high-frequency tail, initial
impulse value1 and zero impulse integral. This matches important surface-transfer
endpoint structure; it can expose a failure hidden by a generic lowpass control.
A positive delay multiplies the transfer by `exp(-s*tau)` and shifts the causal
impulse later. Use an analytic displacement burst with analytic velocity and
acceleration; do not numerically differentiate saved motion to create a reference.

Compare the production pipeline against an independent evaluation of the SAME
finite-band integrals first. Then compare both with the known full analytic
answer. This separates implementation/sign errors from band/window effects.
Require the discrete implementation comparison within `1e-10` of a declared
nonzero control scale; exact-zero cases use a separately declared dimensional
roundoff allowance. Do not demand the truncated answer equal the whole-band
analytic answer to roundoff.

Stop with an attributed failure or ONE next experiment. If analytic controls
explain the residual, fix the demonstrated representation/window defect. If
they pass while a localized BEM residual persists, one fixed-mode/point/frequency
check of existing Green extraction/quadrature choices is appropriate. It is
not permission to begin general mesh qualification. Pole fitting is later:
stable poles alone do not establish accurate physics.

## Proposed candidate qualification, to freeze before validation

These are reviewer-proposed engineering reconstruction tolerances, not existing
user-specified universal accuracy limits. Root must adopt or replace them with
a stated rationale before looking at candidate validation results.

1. Freeze the candidate's frequency/point scope, independent references, fit
   frequencies and held-out frequencies. Start with existing held-out
   `[0.137,0.687,1.337,2.637,5.137,9.937]rad/s`; retain reported failures at high
   frequencies rather than quietly shrinking the scope. Expand validation to
   all six modes at the existing near/far and admissible exterior waterline
   points. Declare interpolation and endpoint extrapolation separately.
2. For EACH point/mode define `S=max|H_reference|` over the frozen validation
   band. Require complex error at EVERY held-out frequency `<=0.01*S`.
   Where `|H_reference|>=0.05*S`, also require amplitude error `<=2%` and
   phase error `<=1degree`. The complex criterion still applies near zeros;
   phase there is undefined and must not be reported as either failure or zero
   error. If `S=0`, use a predeclared dimensional absolute roundoff criterion.
   Never pool mode scales or hide local errors under a global maximum.
3. Manufactured analytic causal tests must use smooth displacement bursts from
   rest, analytic velocity/acceleration, a delayed start, and an independent
   exact/converged reference. Require peak elevation error `<=1mm` AND
   `<=1%` of the reference burst peak when that peak is at least1mm; smaller
   signals retain the absolute criterion. Amplitude, duration and units must
   be frozen before candidate validation. Use a meaningful amplitude rather
   than scaling the control down to pass. Require no pre-response beyond the
   separately frozen absolute numerical tolerance.
4. For BEM-data bursts, require causality, complete startup history and
   documented time-step/band/history sensitivities. Unless an independently
   derived causal time-domain reference exists, do not claim an exact transient
   error bound. The raw signed-lag inverse is NOT such a reference. A fit's own
   state-space evaluation is not independent validation of that fit.
5. Validate saved-motion reconstruction on the full accepted history without
   changing motion or forcing. Preserve incident/diffraction/radiation terms,
   rest/geometry/contact/failure guards and log qualifications. Report the new
   surface change as a change, not as a measured physical error. Independent
   review must accept the actual correction and its scope before publication.

A rational model is an alternative only after the cause is understood. The
[original vector-fitting method and official implementation description](https://www.sintef.no/en/software/vector-fitting/)
support real stable state-space approximation of complex responses. For this
project convention the data belong at `s=-i*omega`; a routine expecting
`+i*omega` needs conjugate data for a real system. Surface elevation versus modal
velocity is not a power-conjugate impedance, so force-impedance passivity
constraints cannot be transferred without a derivation. No dependency is
selected, and no fit is approved by this note.

## Outcome and stopping point

The most economical next action is the stored-data/analytic-control discriminator.
Cause and method choice remain unproven; current provisional warnings remain.
This review stops at the proposal and prospective gates. Execution, production
implementation and issue publication belong to the coordinator/producer owner.

Source-only review of private diagnostic helper
`artifacts/workflow/issue31-stored-diagnostic1/diagnose.py`, SHA256
`8ee71e2b7f9f60ff475ed22224396b00d75714d5cb0974f43e25f870e94ba928`, found no
blocking issue for that bounded diagnostic. Its lowpass/delay ODE controls,
point-specific scales and lag-window comparisons are useful discriminators.
Its descriptive phase threshold is not the prospective candidate gate above;
Heave/two points do not qualify all modes or the full surface. It does not yet
contain the zero-DC analytic control or an independent same-finite-band integral
comparison. Those omissions limit interpretation but need not block the first
stored-data diagnosis. No execution or numerical result is claimed by this review.
