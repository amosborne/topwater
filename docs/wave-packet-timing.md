# Wave-packet timing: a 1D construction example

This illustration addresses [issue #36](https://github.com/amosborne/topwater/issues/36):
use ramps measured in wave cycles while giving different packets the **same
full-amplitude dwell**. The requested CLI duration is intended to mean dwell;
the production CLI has not yet been changed.

The example compares two original wave periods, 4 s and 6 s. Each ramp lasts
two periods, giving ramp durations of 8 s and 12 s. Both envelopes reach full
amplitude at 14.5 s, stay there for exactly 10 s, and start falling at 24.5 s.

| Packet | Wave period | Each ramp | Ramp-up interval | Full-amplitude dwell | Ramp-down interval |
|---|---:|---:|---|---|---|
| A | 4 s | 8 s | 6.5–14.5 s | **14.5–24.5 s** | 24.5–32.5 s |
| B | 6 s | 12 s | 2.5–14.5 s | **14.5–24.5 s** | 24.5–36.5 s |

The different ramp lengths extend outward from the common dwell. Thus total
packet lengths are 26 s and 34 s even though both dwells are 10 s. This is one
candidate alignment, rather than an implemented or qualified test policy.

![Two packets with common dwell, shown as time traces and discrete frequency amplitudes and phases](images/wave-packet-dwell-time-frequency.png)

## Construct the envelopes and waves

For each packet, use the same smooth transition

\[
S(u)=10u^3-15u^4+6u^5,\qquad 0\leq u\leq1.
\]

It joins zero to one with zero first and second derivatives at both ends.
Stretch that transition over the packet's ramp duration, hold the envelope at
one during the shared dwell, then reverse the transition during ramp-down.

Multiplying the envelope by its carrier wave gives the target time trace:

\[
\eta_j(t)=A_j E_j(t)\cos\!\left(2\pi(t-a)/T_j+\phi_j\right),
\]

where a=14.5s is the common dwell start. Both A_j=1m and both specified
phases phi_j=0 in this illustration, so both waves are crests at dwell start.
The envelope, rather than the oscillating height, defines dwell: a zero
crossing during dwell still represents a full-amplitude wave. Phases need not
be adjusted to make the wave end at a zero crossing.

## Phase reference

The specified phases apply at the reference point x=0 and **common dwell start**,
using local time tau=t-a. That is when all components begin the full-amplitude
part of the case with their prescribed phases. The ramps lead into those
phases without resets. A component specified at90degrees crosses zero at dwell
start while its envelope is still at full amplitude.

The example was regenerated to meet this convention; this changes the carrier
waveforms relative to the earlier global-t0 version, while preserving envelopes,
periods and dwell. The older figures remain in Git5aa0170/84fbc8f.

The FFT operates on a global-time sample record. To display each bin's phase
relative to dwell start, use `X_dwell = X_global * exp(+i*2*pi*f*a)`.
Restore the negative phase factor before a global-grid inverse transform.
Bin phases still include the envelope and need not all be zero even when the
specified original carrier phases are zero. The phase panels explicitly use
the dwell reference, not globalt0. Production input conventions have not yet
been changed; this is the agreed convention for the redesigned test.

## Read the frequency plots

A finite envelope has a frequency distribution rather than a single spectral
line. Multiplying it by a carrier shifts copies of that distribution to the
carrier frequencies. The resulting fixed complex spectral coefficients encode
both magnitude and phase; summing them reproduces the packet. The plotted
magnitudes alone are insufficient for reconstruction.

The figure separates the envelope spectrum from the wave-packet spectrum so
that ramp effects can be distinguished from the change of carrier frequency.
Spectra are plotted as bars at the discrete bin frequencies. Their heights are
the amplitudes of the actual cosine components, rather than a connected curve
or a continuous spectral density. Phase panels give the corresponding bin
phases. The inverse transform uses the full complex coefficients and is checked
against the constructed time traces.

Padding samples the spectrum more finely; it does not add physical information
or make a finite packet monochromatic. The discrete amplitudes and bin spacing
belong to the declared192s Fourier-series record and must be used together.

## What this does not establish

These are target **1D time traces**, not a propagated water solution or a WEC-Sim
run. Actual production packets are constructed upstream in space, decomposed
spectrally and propagated with dispersion. Their observed envelope can change,
and finite spectral support leaves tails. The nominal first ramp beginning at
2.5 s is therefore not a certification of the desired calm-water interval.

Issue #36 must still establish actual threshold crossings over the hull,
multidirectional timing, spectral-support sensitivity, measured boat settling,
and the effect of the physical ending on surface reconstruction. The accepted
quintic joins do not by themselves settle those timing questions.

## Executed numerical checks

Both global- and dwell-referenced inverse forms reproduce the new target traces
with maximum error4.44e-16m. An independent selected-time direct component sum
differs by at most2.00e-15m. Both envelopes are exactly one over the same sampled
10s dwell; both specified zero-phase carriers equal1m at dwell start.

The discrete cosine amplitudes are `2*abs(X_k)/N` for positive interior bins;
DC and Nyquist use `abs(X_k)/N` once. The negative-frequency partner is therefore
already accounted for in each positive-frequency cosine amplitude. Envelope
amplitudes are dimensionless and wave amplitudes are in metres. The record is
192s at128samples/s, giving1/192Hz bin spacing. The phase display omits bins
below a stated small-amplitude threshold; that is a plotting mask, not a
modification of the reconstructed packet. All saved complex coefficients remain
in the inverse transform, including bins outside the displayed frequency range.

See the [executed evidence](experts/numerical_methods_researcher/wave-packet-timing-demo-2026-10-27.json) for parameters, normalization and checks. These checks demonstrate the spectral construction, not physical propagation or boat settling.

## Actual sea periods and cycle-based ramp lengths

The current unchanged sea.npz has384 carriers spanning0.03024–0.49481Hz,
or2.02–33.07s periods. The variance-weighted median period is9.80s; the central
approximately90% of source variance spans4.31–15.97s. Weights are A², not an
unweighted count of components or a measurement of energy in the finite packet.

| Cycles per ramp | Median ramp | Central approximately90% variance range | Longest ramp |
|---|---:|---:|---:|
| 1 | 9.8 s | 4.3–16.0 s | 33.1 s |
| 2 | 19.6 s | 8.6–31.9 s | 66.1 s |
| 3 | 29.4 s | 12.9–47.9 s | 99.2 s |

These are seconds **per ramp**, outside the requested dwell. They explain why
two cycles should not become a default without a design reason. No carrier
was removed, rephased or rescaled in this read-only statistics calculation.
[Period/variance evidence](experts/numerical_methods_researcher/wave-ramp-period-stats-2026-10-27.json).

## Preliminary common-ramp bandwidth screen

Equal ramp durations in seconds, with equal dwell and shape, share an envelope
spectrum. Each carrier shifts that spectrum to its own frequency. Equal cycle
counts at different periods do not share the same absolute envelope spectrum.

A cheap reference-point temporal screen used all actual source carriers and
their phases interpreted at dwell start, a60s dwell, and common quintic ramps.
It projected the target onto the existing positive wave-frequency band through
0.533459Hz, with DC explicitly removed. This is not the spatial construction,
not a propagated field over the hull, and not a qualified minimum ramp.

| Each ramp | Maximum target-to-bandlimited elevation difference |
|---|---:|
| 1 s | 98.1 mm |
| 2 s | 13.6 mm |
| 4 s | 7.0 mm |
| 8 s | 2.4 mm |
| 16 s | 0.39 mm |

The existing requested positive-frequency grid has671nodes at omega-step0.005,
independent of ramp duration. A sharper target does not automatically increase
that count: it increases the difficulty of representing the desired target
inside the same band. Expanding the band or refining its sampling changes cost.
This grid count is not a complete BEM/WEC runtime measurement.

The1024-to2048s padding checks at1/2s ramps changed dwell elevations by up to
0.321/0.199mm respectively. These are screen sensitivities, not reference-error
bounds. Calm thresholds and allocation of packet-construction error within the
overall error budget remain to be selected. No source frequencies were pruned.
[Bandwidth-screen evidence](experts/numerical_methods_researcher/wave-ramp-bandwidth-screen-2026-10-27.json).

The next investigation may discard the polynomial constraints. See the
[three researched alternatives](research/minimal-pulse-spectra-2026-10-27.md):
constrained minimax/equiripple, Slepian concentration, and sparse frequency
selection. Their objectives differ; no optimal ramp or numerical savings have
yet been established for this case.
