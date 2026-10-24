# Bounded spectral packet proposal, 2026-10-24

Coordinator question: give the consumer one concrete exact-dispersion incident
construction worth a cheap diagnostic under the already frozen numerical-calm
and first-arrival gates. No numerical execution or consumer implementation by
this researcher. Feasibility is unproven until full-hull eta/velocity/excitation
diagnostics pass; analytical origin estimates below are not acceptance evidence.

For carrier `omega_j=2*pi*f_j`, choose positive half-bandwidth `B_j<omega_j`, and
group-envelope delay `tc_j`. Let

`W_j(nu) = cos(pi*nu/(2*B_j))^4` for `abs(nu)<B_j`, otherwise zero.

Its integral is `3*B_j/4`. At propagation coordinate
`s=x*cos(theta_j)+y*sin(theta_j)`, prescribe

`eta_j(s,t) = a_j*Re[exp(i*phi_j) * integral W_j(nu) *`
`exp(i*(omega_j+nu)^2*s/g - i*(omega_j+nu)*t + i*nu*tc_j) dnu] / integral W_j`.

Each spectral component satisfies the one-way deep-water dispersion relation
exactly. Use that identical field at the reference origin for import and across
space for the displayed incident surface. The spectral group phase delays the
envelope, preserving the canonical carrier phase at the origin. It is not a
global time-ramp multiplication of an existing wavefield.

At time zero, stationary phase puts spectral energy upstream at
`s=-g*tc_j/(2*(omega_j+nu))`. This construction uses a dispersed/chirped initial
packet rather than assuming a narrow envelope rigidly travels at one speed.

Candidate controls communicated before diagnostics:

| Canonical carrier | B (rad/s) | tc (s) | Role |
|---|---:|---:|---|
| 0.1 Hz | 0.30 | 100 | Retained long, late packet |
| 0.25 Hz | 1.47 | 10 | Candidate first arrival within ten seconds |
| 0.5 Hz | 0.18 | 200 | Retained narrow, late packet |

The first-arrival packet's support is approximately 0.0160 to 0.4840 Hz; the
highest packet's support is approximately 0.4714 to 0.5286 Hz. Both respect the
approximately 0.533 Hz deep-water limit from wavelength at least twice beam.
Canonical amplitudes `[.02,.01,.005] m`, headings `[0,45,90] deg` and carrier
phases `[0,pi/2,pi]` remain explicit. These envelope choices are physical input
choices, not a claim of a developed three-component sea within ten seconds.

At the reference origin, the integral has analytic envelope

`E(t-tc) = 4*pi^4*sin(x) / [x*(x^2-pi^2)*(x^2-4*pi^2)]`,

`x=B*(t-tc)`, with removable limits handled analytically. It is one at the peak
and decays as the fifth inverse power away from it. The 0.25 Hz origin initial
envelope estimate is about `6e-4`, with envelope close to one near nine seconds.
Spatial chirp over the hull and velocity/excitation tails may spoil the
initial-calm gate; this proposal must be rejected or adjusted if they do.

Coverage correction before diagnostics: the first unrun suggestion was `B=1.55`,
`tc=9.5`, but its minimum angular frequency was below the producer's `0.1 rad/s`
minimum. The revised `B=1.47`, `tc=10` has minimum angular frequency about
`0.100796 rad/s` and remains inside native coverage. The initial and full-hull
estimates are prospective; neither suggestion has been numerically validated by
this researcher.

Cheapest decisive diagnostic by the consumer: evaluate this spectral sum at the
same frozen frequencies and times for the full hull footprint and submerged
velocity region, including an explicit spatial sampling remainder, and report
the unmodified initial eta/fluid speed. Then evaluate all six initial excitation
tails with native hydrodynamic data and compare the retained-history native IRF
path against an independently evaluated spectral transfer diagnostic. Long late
packets must not create a false pass merely because the finite retained kernel
truncates physically relevant noncausal support. No initial samples or force
tails may be clipped, and body motion must start at exact rest without cropping
a moving run. Report all three centers/onsets even when their peaks occur later.

Stopping point: this proposal and independent review of owner diagnostics. No
claim of a feasible packet or motion result yet.

## Independent first executed diagnostic review

Inspected owner evidence
`artifacts/reference/minimum-20261024/packet-attempt2/diagnostic.json`.
The candidate is correctly REJECTED under unchanged prospective gates:

- Full hull bounding volume: 398,898 samples, maximum spacing 0.025 m,
  elevation bound `12.346 um < 35 um`, velocity bound
  `42.824 um/s < 43.982 um/s`. Bounds include local sampled gradients and a
  global spectral Hessian remainder, not vertices alone.
- First packet arrival: 92.2 percent canonical 0.25 Hz amplitude within ten
  seconds. All three carriers retained. Entire spectral support is within the
  selected native hydrodynamic band and wavelength restriction.
- Initial native/direct heave force `0.37135/0.37061 N` and pitch torque
  `0.13014/0.12827 Nm` exceed the frozen absolute 0.1 limits. The direct-H check
  rules out treating omitted noncausal kernel support as the explanation for a
  false pass. All six carrier force scales are nondegenerate.

No body was simulated and this result does not permit motion acceptance. A
prospective fixed-order cheap screen at focus times `[11,11.5,12,12.5] s`, same
bandwidth/carriers/gates, is an appropriate next physical packet-input test.
Screen initial native/direct forces and arrival, then evaluate full hull-region
eta/velocity on the selected passing candidate. No gate relaxation, tail
clipping or fitting to body motion is implied. A feasible candidate remains
unproven.

## Selected packet executed review: wet region passes, full planform pending

Inspected `packet-attempt3/diagnostic.json`, source snapshots and frozen manifest
under `artifacts/reference/minimum-20261024`. The selected first passing
fixed-order candidate has focus times `[100,11.5,200] s`, unchanged bandwidth,
canonical carriers and calm gates. Execution ledger records terminal exit zero.

Recomputed hashes of current source, controls, canonical physical input, corrected
native netCDF and consumer BEMIO matrices all match the frozen manifest. Source
review confirms no wave or force samples are clipped/overwritten; the sole
endpoint operation makes the analytic zero spectral weights exactly zero.

The declared wet-hull bounding volume passes: elevation bound `3.3055 um`,
velocity bound `24.4842 um/s`, all six retained-IRF and direct-H initial forces
below their prospective relative/absolute limits, arrival `59.4093%` for the
0.25 Hz carrier, and allowed physical frequency band. The maximum native/direct
initial force is approximately `0.07430 N` (direct heave). No body is simulated.

The continuous-region bounds are mathematically adequate for their declared
domain: sample spacing at most 0.025 m, sampled local gradients, and global
spectral Hessian Taylor remainders. The 17 triangular faces and all their
interiors lie in that box; coverage is more than vertex samples.

Coverage limitation: the box `x=[-2.8233,3.3659]`, `y=+/-1.09854`,
`z=[-0.402142,0] m` derives from the wet mesh and excludes above-water hull
overhang in the full 24 by 9 ft hull planform. A supplemental calm-only bound on
the full mean-world hull planform with the same wet-depth range is recommended
before unqualified full-hull startup acceptance. Keep the identical packet and
gates; existing force/arrival evidence need not be recomputed or refitted.

The retained-IRF versus direct-spectral heave force maximum difference is
`0.0286456 N` over the first ten seconds. Report this finite-band/kernel
diagnostic limitation; it is not independently a motion-accuracy acceptance
gate. The existing evidence does not yet establish body initial state or video.

## Final independent incident-input startup signoff

PASS for numerical calm, forcing tails, physical band and arrival of this frozen
incident input. Inspected supplemental
`packet-full-planform/diagnostic.json`, frozen source manifest, archived/current
sources and terminal ledger. Recomputed all six manifest input hashes and the
three archived source/control hashes: all match.

The full hull is transformed by subtracting material CG, applying the canonical
active body-to-world quaternion, and adding mean world CG. Its bounding planform
`x=[-2.8249821,4.4901874]`, `y=+/-1.3716 m`, combined with full wet-depth range
`z=[-0.402142,0] m`, covers the complete planar triangular hull projection and
submerged region. The 587,412-point grid and reviewed Taylor remainder bound the
continuous region at elevation `3.323565 um` and fluid speed `25.495529 um/s`,
below unchanged `35 um` / `43.982297 um/s` gates. Spectral `c` and `w` arrays are
asserted exactly equal to accepted packet attempt three; there is no forcing
change or refit. Terminal exit zero is recorded for exact session 80979.

Together with attempt three's all-six native/direct initial-force gates,
physical frequency band and 59.4 percent first arrival, this closes incident
input startup review. It does not prove exact body q/v/radiation-memory initial
conditions, actual Simulink force application, motion or video. The native
excitation-mask multiplication issue discovered subsequently must be resolved
and independently reviewed before body reliance.
