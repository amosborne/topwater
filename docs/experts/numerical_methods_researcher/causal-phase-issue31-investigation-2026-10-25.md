# Radiation-surface phase investigation: issue 31

Scope: `topwater solve --solver=wecsim`, specifically Topwater's surface adapter.
This is separate from the future native solver and is not an established
upstream WEC-Sim defect. Related independent issues: [reference resolution](https://github.com/amosborne/topwater/issues/32),
[field/load consistency](https://github.com/amosborne/topwater/issues/33), and
[linear applicability](https://github.com/amosborne/topwater/issues/34).

2026-10-25. Independent source/literature/report review only. Read
[issue 31](https://github.com/amosborne/topwater/issues/31) including its complete
body; it had no comments when read. No numerical calculations, imports, tests,
BEM/body jobs, production edits or GitHub writes were performed. The active
canonical attempt 4 continues under the coordinator's separate lease. This note
proposes a later bounded discriminator, not another execution allocation.

## Finding

The discrepancy remains unresolved. No new sign, normalization or omitted-field
error was established. The provisional warning remains appropriate. Finite-band
projection and finite lag-window effects should be separated before fitting a
replacement or attributing the residual to the BEM representation. The existing
small spacing/cutoff sensitivities do not distinguish those possibilities.

Source reviewed: `_capytaine_surface.py` SHA256
`2858b847efea11a97c50bd9d12c744bf3069723e32e5fa401309c4e2ee0f352a`.
Capytaine **3.0.0** was verified from installed distribution metadata and
`uv.lock`; no package import was used. Historical response producer identity is
separately frozen in the
[response manifest](../../../artifacts/workflow/surface-response-attempt1/manifest.json).
Its source hash is different; current geometry/logging changes are not evidence
that the historical radiation operator was corrected.

## What the evidence establishes

The [response report](../../../artifacts/workflow/surface-response-attempt1/field/report.json)
has near-point per-mode full-complex recovery discrepancies about 17--30.4%.
The denominator is each mode's maximum transfer magnitude over the sampled
points/frequencies. This is neither a pointwise percentage nor a total-elevation
error. Spacing sensitivity is at most 3.26 micrometres and the tested cutoff
sensitivity at most 0.753 mm in the actual-motion radiation response. The
discrete infinity acceleration candidate is zero. Shared-Green/public evaluator
agreement to 1.72e-16 verifies equivalent evaluation paths using the same model;
it does not independently validate the transfer data.

The [projection report](../../../artifacts/workflow/surface-projection-attempt1/report.json)
shows full signed-lag versus cosine differences up to 42.665 mm maximum and
12.511 mm RMS, persisting before 55 s, with future-velocity contribution up to
77.265 mm. Its raw finite-band inverse is not a causal reference or an error
bound. The older
[far-point dense report](../../../artifacts/workflow/refinement-checkpoint-attempt1/refinement/report.json)
separates real and imaginary residuals: cutoff 24 / spacing .05 has maxima
0.000299981 and 0.00532523 respectively. That supports investigating phase
reconstruction, but does not supply the missing componentwise near-point audit.

## Conventions and a useful finite-band identity

Capytaine uses displacement-normalized radiation, time dependence
`exp(-i omega t)`, and `eta = i omega Phi/g`. Its installed radiation boundary
condition is `-i omega * displacement_on_face`. Therefore
`Hv = i eta_q/omega` is the correct modal-velocity transfer. These conventions
are explicit in the official [theory manual](https://capytaine.org/stable/theory_manual/theory.html#forces-on-body-surfaces)
and [surface-elevation equations](https://capytaine.org/stable/theory_manual/theory.html#free-surface-elevation).

For a real causal regular kernel, the project's convention gives

```text
H(omega) = integral_0^infinity h(s) exp(+i omega s) ds
h_raw(s) = integral_0^Omega [Re(H) cos(omega s) + Im(H) sin(omega s)] d omega / pi
h_cos(s) = 2 integral_0^Omega Re(H) cos(omega s) d omega / pi
```

Direct algebra, with the SAME frequency quadrature, gives
`h_cos(s) = h_raw(s) + h_raw(-s)` for positive lags. Thus the cosine construction
folds the raw negative-lag part into positive lags; it does not merely discard
negative lags. This explains why the signed-lag positive contribution and causal
cosine response can differ substantially without a convolution sign error.

Also by Fourier algebra, a sharp band truncation convolves a causal kernel with
a two-sided sinc kernel. Negative-lag ringing can therefore occur for perfectly
causal input data. It is not by itself evidence of acausal BEM physics. Causal
real/imaginary relations require the relevant frequency coverage and endpoint
terms; the original ship-motion
[Kotik--Mangulis paper](https://doi.org/10.3233/ISP-1962-99701)
states its added-mass/damping relation using knowledge over all frequencies.
That force result motivates the consistency check; it does not independently
prove a surface-height transfer or its discrete approximation is qualified.

## Ranked questions and cheapest discriminators

| Priority / possible cause | Established limit | Bounded next discriminator, after allocation |
|---|---|---|
| 1. Finite-band and lag-window reconstruction | The kernel uses only real transfer data through 24 rad/s; `recover_transfer` integrates lags only through 60 s. A kernel covering all requested transient history need not cover its entire Fourier-transform tail. | Use saved transfers for Heave at the existing near `[0,1.25]` and far `[6,0]` points. Separate held-out real/imaginary residuals by frequency and point. Compare the current lag window with ONE longer diagnostic window safely below `pi/DeltaOmega`, independently of the already saved cutoff 12/24 choices. Do not change the production duration or motion. |
| 2. Band-end treatment versus underlying transfer inconsistency | The same sharp cutoff can ring; a taper changes the represented response. Neither extra zero padding nor a small quadrature-step delta supplies absent high-frequency information. | Run the identical inversion/recovery diagnostic on a known causal analytic transfer, such as `1/(a-i omega)` with positive `a`, plus its delayed version `exp(i omega tau)/(a-i omega)` with positive `tau`. Record band/window-induced real/imaginary and pre-response residuals. A taper may be a comparator, not a hidden repair or accuracy proof. |
| 3. Near-hull discrete transfer/quadrature and irregular-frequency limitations | The unchanged 17-panel representation is unqualified. Capytaine's current producer uses default panel quadrature; mesh/irregular-frequency warnings occur. Near points are more sensitive than the old far-point sample, but that does not establish a cause. | Only if the stored-data check leaves a localized discrepancy, select ONE mode/frequency/near-far pair for an independently allocated numerical discrimination. Change one integration/discretization choice, preserving physical geometry and the saved motion. Do not launch a mesh or high-frequency campaign by default. |
| 4. Endpoint, direct or distributional terms | The exact discrete acceleration endpoint is zero. This does not establish the rate/uniformity of finite-frequency approach, all remaining tails, or every possible direct term. Increasing finite-frequency magnitudes are insufficient grounds to fit a derivative feedthrough. | Inspect fixed point/mode component trends and analytically derived endpoint quantities, not maxima that can switch mode/point. Reject an endpoint prescription if it contradicts its derivation or held-out complex data. Low-frequency pole terms likewise require a derivation. |
| 5. Conversion, convolution or input-frame error | Source inspection finds the expected `Hv`, positive-sign forward recovery, endpoint half weights, zero prehistory and modal velocities. It does not make numerical implementation infallible. | The analytic control above cheaply checks signs, factors, delay direction and rest startup. Keep modal units separate; do not compare rotational and translational coefficients under an unlabelled common absolute scale. |

The official [Capytaine quadrature documentation](https://capytaine.org/stable/user_manual/mesh.html#defining-an-integration-quadrature)
distinguishes analytical Rankine integration from numerical wave-part integration
and provides the optional Gauss--Legendre 2 scheme. This supplies a possible
targeted numerical discriminator, not evidence that enabling it will fix this
case. The documented [infinite-frequency Green-function limit](https://capytaine.org/stable/theory_manual/inf_depth_green_function.html#singularities-extraction)
is consistent with the already evaluated zero acceleration endpoint at regular
exterior free-surface points. It does not justify inventing a finite-band tail.

## Candidate methods, only after diagnosis

An acceleration formulation is algebraically possible:
`Ha = Hv/(-i omega) = -eta_q/omega^2`. With matching zero-prehistory conditions,
it must represent the same output. It weights different transfer components
when inverted and requires its own low-frequency limit, high-frequency direct
term and integration-by-parts boundary treatment. It is a useful complementary
diagnostic, not an automatically correct replacement. Do not differentiate a
noisy saved velocity to manufacture acceleration or change the accepted motion.

A real stable rational/state-space approximation can fit BOTH components. With
standard Laplace notation it must satisfy `G(-i omega) = Hv(omega)`; a fitter
expecting data on `+i omega` instead needs the corresponding conjugate data for
a real system. Use conjugate pole/residue pairs and left-half-plane poles.
Stability is not physical accuracy, and velocity-to-point-elevation is not a
power-conjugate impedance: do not impose positive-real/passivity constraints
borrowed from radiation-force models without a separate derivation.

The original [Gustavsen--Semlyen method](https://www.sintef.no/globalassets/project/vectfit/vf_paper.pdf)
and [official SINTEF description](https://www.sintef.no/en/software/vector-fitting/)
support rational approximation and stable state-space realization. They do not
validate extrapolation, a particular model order or the input BEM data. The
[VFIT3 implementation](https://www.sintef.no/en/software/vector-fitting/downloads/vfit3/)
also has an explicit commercial-embedding restriction; no implementation reuse
or dependency addition is selected here.

## Proposed stopping point and acceptance

First allocate only the stored-data/analytic-control discriminator above after
the current process is terminal and its lease released. It needs no BEM or body
rerun. Stop with per-point/per-mode held-out real/imaginary residuals, window/band
attribution limits and a recommendation for ONE next action. A failed causal
fit alone would not prove noncausality of incompletely sampled data.

Before validating any replacement, the coordinator must declare amplitude/phase
and transient-elevation tolerances, the frequency/point scope, and absolute
handling of near-zero responses. Hold out the existing direct frequencies from
fitting. Then require both transfer components, causal smooth-burst startup and
actual saved-motion response checks; expand across six modes only after the
initial discriminator is informative. Relative phase at a response zero is
undefined, not a failed angle comparison or a zero-error pass.

Keep the six-DOF motion, forcing, incident train, diffraction and physical
radiation normalization immutable. Field/load discrepancy, moving-contact
continuation and export geometry are separate limitations; none explains away
a radiation-only held-out transfer residual. Retain the current provisional
warning until independent evidence supports a replacement. No pause/restart or
new numerical campaign follows from this note.
