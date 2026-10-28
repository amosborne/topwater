# Minimal spectra for a finite wave event

**Latest user choice:** use an optimizer-free fixed-quintic construction and
retain the largest N Fourier coefficients. The requested
[N10/20/40 comparisons](pulse40-three-examples.md) are generated with phases
defined directly at dwell midpoint. The optimization proposals below remain
research alternatives, not the selected implementation or a queued assignment.

Research proposal, 2026-10-27. A subsequently requested bounded
[minimax chart comparison](minimax-packet-comparison-2026-10-27.md) is now executed;
it is an illustration, not a qualified production wave generator.
Scope: [issue36](https://github.com/amosborne/topwater/issues/36), allowing the
polynomial shape and derivative constraints to be discarded. Current production
inputs and accepted results are unchanged. The coordinator owns experiments.

## What should be minimized?

The latest clarification separates the desired sea-carrier grid (proposed
0.01–1Hz/.001Hz) from packet-construction offsets and spacing. Each carrier
gets its own deterministic packet based on base frequency and requested dwell,
under a documented policy. Do not impose the carrier range as a hard final-
component bound. All generated physical components still need consistent
frequency/direction handling and adequate hydrodynamic coverage. No production
policy is implemented or qualified by the earlier illustrations.

The useful target is a nearly unit **complex envelope** during the requested
dwell and small elevation tails over declared quiet intervals. It is not
constant real water height: the carrier continues oscillating, with its specified
phase at dwell start. A rounded short pulse does not automatically meet this
flat-dwell requirement.

Three different costs must remain separate: occupied bandwidth, number of
retained spectral components, and actual producer work. A narrow band can contain
many finely spaced lines; a few separated lines can span a wide band. Frequency
components may combine before hydrodynamic work, so line count alone does not
prove runtime savings. Producer cost mapping needs separate review.

There is no exact isolated, flat-topped pulse made from finitely many frequencies:
a finite exponential sum is analytic, and equality to a constant on an open
interval makes it that constant everywhere. Lines on a common Fourier grid repeat
after its period; arbitrary discrete lines do not provide permanent
calm. Define a finite observation/quiet horizon, allowed dwell ripple, allowable
tail height, transition gap and coefficient/overshoot bounds before asking for
a minimum. Removing derivative constraints does not remove this tradeoff.

## Top three options

| Option | What it optimizes | Tradeoffs for this test |
|---|---|---|
| **1. Constrained minimax / equiripple: first recommendation** | The worst dwell/quiet error at a chosen spectral budget. Related FIR design and bounded complex-sum synthesis are established: [Remez](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.remez.html), [Lebret–Boyd](https://stanford.edu/~boyd/papers/ant-pat.html). | Directly specifies a nearly flat complex dwell and maximum tails; transitions are free. Search bandwidth or transition length separately. Moderate offline optimization, ordinary spectral summation afterward. Usually dense: it does not minimize line count. Reject if useful bounds require excessive bandwidth, overshoot or coefficient magnitudes. |
| **2. Slepian / DPSS concentration** | Integrated energy concentration in prescribed time/frequency intervals, with an established eigensystem implementation: [Slepian–Pollak](https://doi.org/10.1002/j.1538-7305.1961.tb03976.x), [SciPy DPSS](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.dpss.html). | Attractive bandwidth/energy-tail benchmark, but a leading pulse is rounded rather than flat. Small integrated leakage does not bound peak tails. A constrained combination can enforce dwell/quiet limits, adding optimization and possibly losing its compression advantage. Fewer basis vectors need not mean fewer physical frequencies. Reject on peak/ripple failure regardless of concentration ratio. |
| **3. Sparse frequency selection** | Number of retained frequencies through a constrained l1/reweighted-l1 surrogate: [Candès–Wakin–Boyd](https://web.stanford.edu/~boyd/papers/rwl1.html). | Keep the same peak dwell/tail constraints while promoting zero coefficients, then prune and refit. Most relevant when active lines reduce actual work. Additional offline solves, no generic minimum-count guarantee, and potentially fragile cancellations. It need not minimize bandwidth. Reject if savings vanish after refitting/off-grid checks or do not remove producer work. |

These are three distinct objectives, not three guarantees of the same optimum.
Options1 and3 adapt established optimization ideas to this project; no cited
paper establishes their wave-packet performance. Dolph–Chebyshev is a useful
special-case minimax comparator, not a fourth recommendation: its classical
main-lobe/sidelobe optimum does not independently prescribe our dwell width and
ripple. [Official definition and endpoint caveat](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.chebwin.html).

SciPy documentation consulted is v1.18.0; these are reference implementations,
not a claim about the installed environment or a request to change dependencies.

## Recommended first formulation

For fixed envelope offsets nu_k, write

`g(tau) = sum c_k exp(i*2*pi*nu_k*tau)`, with `tau=0` at dwell start.

Choose coefficients to satisfy `|g-1| <= eps_dwell` during dwell and
`|g| <= eps_quiet` in the stated early/late quiet intervals. Fix `g(0)=1` to
anchor the specified carrier phases. Allow a transition interval rather than
imposing a quintic. Bound overshoot and coefficient norm to prevent an apparently
excellent finite-grid result bought by huge cancellations outside that grid.
For a selected bandwidth/count, feasibility versus transition length gives a
meaningful shortest-transition question; bandwidth minimization is a separate
outer search.

Before trusting grid constraints, check between nodes. The elementary bound
`|g'| <= 2*pi*sum |nu_k*c_k|` bounds the missed excursion by that slope times
half the grid spacing. Dense sampling alone is evidence, not a continuum proof.

Carrier j uses `A_j*g(tau)*exp(i*(2*pi*f_j*tau+phi_j))`. Its absolute frequencies
are `f_j+nu_k`; impose the actual physical frequency limits on those, not merely
on envelope bandwidth. Band-edge carriers can restrict a common symmetric
envelope. Preserve directions and the real-field conjugate representation, then
evaluate the dispersed spatial field over the hull: a successful scalar trace
does not certify spatial calm. An antenna-array analogy supplies mathematics,
not water-wave physics or a reason to change forces/damping.

## Bounded first experiment to propose

Use the closed-form truncated Fourier series of a rectangular dwell as the
**simpler baseline**, with the same allowed transition gap. It exposes ringing
without another arbitrary polynomial. Place the rectangle's jumps inside those
transition gaps, not at the dwell
endpoints, where the Fourier series would approach half-height and create an
unfair endpoint error.
Compare it with constrained minimax at identical dwell, quiet horizon, phase
anchor, frequency dictionary and tolerances.
First use one carrier; then crossed carriers and the actual directional input.
Report peak dwell error, peak quiet height, overshoot, coefficient norm, occupied
band and retained-line count separately. Translate normalized ripple into metres
against an explicit amplitude bound, not an accidental cancellation at one point.
Reject minimax before BEM/WEC if it fails the allocated height/calm limits or
offers no useful improvement over the baseline. A DPSS candidate must meet the
same peak limits; only then compare its energy-tail benefit. Attempt sparse
selection only after a feasible design exists, and reject it unless retained
frequencies can actually be omitted from production work. No optimum ramp seconds
or count is claimed until tolerances and these checks are fixed.
