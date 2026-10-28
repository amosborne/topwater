# Shared physical and producer/consumer contract

## Per-component packet interface, corrected user decision2026-10-27

Inputs:basefrequencyf0,phasephi DIRECTLYatdwellMIDPOINT,rest,ramp,dwell,andN.
Timing:rest→ramp→dwell→ramp→rest; total2rest+2ramp+dwell. Fixedquintic
ramps are reversed copies; no optimizer. Each sea component retains amplitude,
phase and direction; the current examples are1Dtime traces only.

N counts POSITIVE-frequency components. DC and negative-frequency conjugates
are free. Use the REAL target a(t)*cos(2pi*f0*(t-mid)+phi), preserving its mean
and Hermitian Fourier pairs; no amplitude/phase rescaling. Compute coefficients
directly, retain N largestpositive magnitudes and their conjugates. Bandwidth,
placement and achievedquality are outputs. No imposedfinalband/sharedenvelope
acrosscarriers. Current example:f0.01/.1/1Hz,phi0/45/90at65s,rest10/ramp5/
dwell100/ramp5/rest10,N10/20/40. A phase reference shift is not an inputfeature.

The previous complex positiveprojection dropped necessary signed content and
was wrong for this reconstruction task. Preserve that study under itsoriginal
policy; do not describe its projection loss as physical necessity. Correctreal
coefficient:Cpos=e^(iφ)G(f-f0)+e^(-iφ)G(f+f0);DC=G(f0)cosφ. Store fullFFT
coefficients Cpos/2 with conjugates. Validate against an independent realtarget
FFT and checkIFFTimaginaryresidue/inverseagreement. FullrealFourier converges;
finiteN introduces measuredtruncationerror. No predeterminedphysicalqualitypass.

Envelope symmetry belongs to the prescribed a(t); fulloscillatingwave neednot
beeven for arbitraryphi. Its spectrum is Hermitian aboutzero; no extra paired
frequency requirement aboutf0. Do not mistake analytic-signal magnitude for
the prescribed amplitude envelope in broadband/lowcarrier cases. Plot actual
realreconstruction against actualrealtarget. Generatedphysicalcomponents still
need producercoverage/spatialreview before productionuse. Canonicaldata,
production solver and maintainedWECfork are unchanged by theseexamples.

## Simulation timing integration remains separate

CLI requested time means dwell. Actual body settling must be measured; it is
not implied by the requested incoming-wave rest interval and need not return
the body to its original position/heading. Future integration needs an explicit
did-not-settle outcome/cap and unchanged physical startup/contact guards.
Do not add damping merely to force rest or infer calm from one zero crossing.
Current runtime timing policy remains unchanged until separately implemented.

## Surface reconstruction accuracy, user decision 2026-10-27

For issue31, the user requires maximum numerical water-elevation error <=0.01m,
INCLUDING reference numerical uncertainty, against the chosen linearized model.
Preserve the completed WEC-Sim trajectory/forcing and the chosen mean-hull data;
this requirement does not ask for comparison with nonlinear physics or experiments.
Compare actual total elevations at the same coordinates/clocks over the full60s
and declared valid exterior domain. No phase alignment or motion modification.

The reference must preserve complex response amplitude/phase and qualify startup,
frequency/time integration, padding/future independence and numerical uncertainty.
A finite-band signed-lag output or agreement of two unqualified approximations
alone is insufficient. Report unresolved tail uncertainty rather than call it zero.
A few-probe screen is not a whole-region maximum. Keep the continuation policy
explicit; never sample invalid mean-hull interior Green functions silently.
Internal transfer-function percentages and derivative tolerances are diagnostic
checks, not replacements for this accepted height requirement. They do not change
this user's1cm criterion or qualify mesh/model accuracy tracked separately.

## Canonical NPZ workflow extension, 2026-10-25

### Authorized WEC-Sim force-ownership transition

Latest user decision, 2026-10-26: pursue native spectral excitation instead of
imported-elevation convolution. WEC-Sim receives the complete physical wave
components and computes their harmonic force sum. Topwater may construct and
encode physical components and coefficient inputs, but must not supply applied
force histories. Preserve the original weighted complex component coefficients,
frequencies, exact propagation headings, origin, absolute time and phases.
The coefficients already include packet quadrature and dispersion Jacobian;
do not add those weights again or replace coherent components with energy averages.
Equivalent spectrum-density/phase encoding is acceptable only if it reconstructs
those same components, without random phases, renormalization or a temporal ramp.

The first bounded native test uses three nonuniform frequencies and two headings,
explicit phases and nonunit directional-bin widths. Independent expected elevation
and force sums use the prescribed integrated amplitudes/phases directly, not the
production density/spread calculation. Include heading-dependent and purely
imaginary transfer values. This test excludes low-frequency and radiation changes
so failures isolate phase/sign/weight handling. Verify the actual native force path.

For the later unchanged3001x6 canonical spectral force replay, maximum error must
be <=1e-10 times each reference channel's absolute peak and time-weighted RMS
error <=1e-10 times its reference RMS. Use trapezoidal integration over0..60 s,
all endpoints and the original basis, with no fitting, mean removal or denominator
floor. Exactly zero reference channels require exact zero by default; any separate
dimensional roundoff allowance must be independently frozen before execution.
Report component-sum/peak and component-sum/RMS conditioning diagnostics without
using them to waive a failed gate. These prospective numerical replay limits were
independently reviewed before the new spectral route results.

Preserve the old complex bilinear coefficient interpolation and quadratic
low-positive-frequency bridge. Native spectral excitation removes the need for
an excitation IRF; the existing radiation CIC, Ainf, dry mass/inertia, restoring
operator, initial rest/zero memory and physical startup limits remain unchanged.
Actual force application must still be checked exactly once before motion.
The subsequent motion compatibility limits below remain applicable. The older
convolution handoff and coarse/refined kernel gates below describe retained
studies, not requirements to carry excitation convolution into the spectral route.
Dependency promotion and merges remain subject to the user's existing holds.

The user now requires WEC-Sim to calculate wave excitation as well as boat
motion. Topwater prepares boat/sea inputs and coefficient data, runs the tools,
and reconstructs the surface after motion. It must not synthesize the applied
wave-force histories in the production WEC path. The previous spectral force
calculation may remain only as an offline comparison oracle. This concerns the
WEC-Sim integration; Topwater-native remains a nonworking prototype.

Preserve the complete sea realization, amplitudes, headings, phases, spatial
reference and time zero. Directional elevation records, including opposite
propagation branches, must reconstruct the unchanged incident elevation at the
native origin over output and padded times. Group only identical headings. Do
not replace direction-specific records with one scalar sum. Preserve true local
calm/rest/zero radiation memory and zero artificial ramp at startup.

The agreed standard-H5 handoff first computes radiationIRF from the unchanged
0.1..12 rad/s data, retaining its kernel/time/frequency arrays and direct Ainf.
The original transition baseline used genuine common raw tables on0.1..3.4 rad/s
with all25 periodic headings. The subsequently accepted full-band candidate uses
genuine0.1..12 rad/s raw tables at all25 headings, including FK/diffraction;
its common low-band data and full radiation objects are exactly preserved.
Official H5 readback verifies that expanded candidate. The excitation kernel
uses the genuine operator and explicit low-frequency bridge on its own grid;
no unsolved coefficients or new reader schema are invented. Readback must prove
full-band radiation memory and added mass are retained. Candidate acceptance
does not promote the held runtime or qualify whole-case forces and motion.

First compare WEC-computed forces with the offline old spectral reference on
identical coefficients and wave records, including startup and endpoints, and
verify actual downstream force application exactly once. Quantify each force and
moment channel's initial/peak/RMS/maximum differences and numerical-window effects.
Then compare all six motion histories with the unchanged public baseline. Do not
assume finite-kernel convolution is identical to direct spectral synthesis, silently
relax the requirement, overwrite canonical results, or promote the held pin before
the differences are reviewed. The eight-file public example contract is unchanged.

### Prospective canonical WEC-Sim compatibility gates

These are coordinator-selected engineering agreement criteria, independently
reviewed before the new actual-case results. They qualify the fixed canonical
transition only; they are neither general CLI requirements for other boats nor
physical-accuracy guarantees. Finite-probe kernel admission does not guarantee
them for the packet's complete spectrum and directional cancellations.

Use the unchanged 0..60 s reference, dt .02, all original samples/endpoints and
declared frame/origin mappings. Do not fit a time shift, phase, gain or rotation.
For each force/moment channel, define P as the maximum absolute reference value
and R as its time-weighted RMS, using trapezoidal integration over the full 60 s
without mean removal. Both coarse and refined candidate histories must satisfy
maximum absolute error <= .001 P and RMS error <= .001 R. Refined-minus-coarse
sensitivity must satisfy .0005 P and .0005 R using those same reference scales.
An exactly zero reference channel requires exact zero unless an independently
derived dimensional roundoff allowance was frozen beforehand; no denominator
floor or after-the-fact relaxation is permitted.

Report every channel's initial/final value and error, signed extrema and times,
absolute peak, RMS, maximum error and its time, and RMS error. Tiny nonzero
reference channels retain their original denominators and are identified.
Existing local calm, force/moment startup, initial rest, zero radiation memory
and zero artificial-ramp gates remain separately mandatory. Complete padded
history must cover both simulation endpoints and the declared kernel support.

Actual downstream wave-excitation samples must match the prepared native total
exactly once: maximum difference <= 1e-10 times that prepared channel's peak,
with exact zero for a zero prepared channel. This implementation replay gate
cannot substitute for reference-force agreement. On a force-gate failure, stop
before motion and inspect retained worst samples and a few independent scalar
directional contributions; do not tune physical inputs or budgets.

After force acceptance, maximum world-position difference vector norm must be
<= 1 mm and maximum relative-orientation rotation angle <= .01 degrees against
the saved canonical trajectory. Use rotation matrices or normalized quaternions
without fitting an alignment. Report all six raw motion channels, velocities,
initial/final differences and peak-error times. Hold inertia, radiation operator
and integration settings identical; isolate any separate consequential dynamics
correction before attributing trajectory changes to excitation. Dependency
promotion remains subject to the user's existing hold and merge instructions.

Latest user authorization is workflow implementation, not native dynamics.
Current physical inputs are examples/simple-boat/boat.npz and sea.npz, preserving
the accepted boat and stronger measured spectrum. Public data are numeric NPZ;
the current array/frame contract is docs/npz-format.md. Historical references
below retain their original physical inputs and policies.

All application boats are convex physical triangular hulls with exactly one
inferred opening. Caps are internal only. Result exports a complete hull snapshot
(possibly remeshed), CG/absolute pose history and mandatory total water-air
geometry. Native nonlinear output is not constrained to a linear decomposition.
The WEC adapter reconstructs its linear incident/diffracted/radiated field from
the existing model and computed motion without feedback or another body solve.
Latest user steering: preserve that total field and its producer mesh on export.
Do not cut or remesh water triangles to fit the hull. Portions inside the posed
hull are invalid fluid and are implicitly excluded. Compute the physical-hull
waterline intersection in a separate shared geometry routine. The format still
accepts a native solver's own evolving surface mesh; it imposes no linear field
decomposition or fixed-grid requirement on that solver.
Geometric moving-contact continuation is an explicit approximation, not nonlinear
waterline physics. Quantify total-minus-incident surface delta and force-field
consistency; do not tune to a small delta or alter accepted body forcing.

Physical applicability departures warn at the agreed10% restoring/wetted-set
engineering reporting thresholds where the backend makes those assumptions.
Fatal errors prevent a defined complete result and preserve previous outputs.
No output metadata/solver marker or physical-accuracy certification is implied.

## Scope

The incoming-packet/startup contract below belongs to the completed REF-001
minimum reference and remains immutable for that result. Latest user steering
authorizes REF-002, one fresh buoy-generated strong-motion reference, described
at the end of this file. Do not require REF-002 to preserve REF-001's sea.

Reference requirements (authorization is in the [handoff](../checkpoint-handoff.md)):
standalone Capytaine -> official BEMIO -> WEC-Sim external reference. Ten-second six-DOF motion, isometric video, then variation of those
motions with solver mesh. Topwater solver implementation is held.

Authoritative neutral physical input:
[canonical-boat.json](../../reference/capytaine/canonical-boat.json).
It carries the minimum full24ft hull and frozen17-triangle wet surface, mean
equilibrium pose/world CG, full dry mass/inertia including products, canonical
three-component mild sea and provenance hashes. Density1025kg/m3, gravity9.81m/s2,
infinite depth, zero forward speed. Keep minimum full hull geometry separate
from derived solver meshes; remeshing preserves the physical surface, equilibrium
and mass properties. The retired canoe and archived P1 do not qualify this boat.

World axes: x forward, y port, z up; right-handed. Modes Surge/Sway/Heave/Roll/
Pitch/Yaw, rotations about mean WORLD CG. Capytaine complex convention
Re(amplitude*exp(-i*omega*t)); headings are propagation radians CCW from+x.
Capytaine excitation is diffraction PLUS Froude-Krylov per unit elevation.
WEC-Sim directions are degrees, with its own phase convention and WAMIT
normalization. Its official readCAPYTAINE negates excitation imaginary parts;
verify the interface once, never manually conjugate a second time.

Producer, accountable capytaine_expert: native netCDF, native hydrostatic files,
physical body/mesh/CG/inertia/DOF identities, exact frequencies/headings, versions,
code/input/source hashes and reproducible export commands. Current producer
readiness and unfinished work are in the producer current-state and handoff.

Consumer, accountable wecsim_expert: official BEMIO parsing/normalization and H5,
hydrostatics and full inertia/frame mapping, transient excitation/radiation
handling, freely floating six-DOF WEC-Sim case, native output grids/rotation
representation, incident visualization and video. Both experts sign off the
small interface with a reproducible example. The main coordinator does not
reverse-engineer their boundary from raw code.

Startup: boat at rest in locally calm water. Waves begin UPSTREAM and roll in;
no global amplitude/force step beneath the hull, no imposed temporal ramp for
the first result. Freeze the incident packet/train, leading edge and actual
arrival timing consistently with transient excitation. Retain startup; settling
inside10s is not required. A moving visual wavefront over globally forced or
steady-state motion is invalid. Label prescribed incident surface distinctly
from a computed scattered/radiated fluid field.

Minimum-case startup interpretation, frozen prospectively 2026-10-24: finite-band
deep-water packets have noncompact tails. Report those tails without clipping.
Numerical local calm at t=0 requires hull-region maximum |eta| <= 1e-3 sum(a)
(3.5e-5 m), maximum fluid vector speed <= 1e-3 sum(a*omega)
(4.39823e-5 m/s), and every excitation component |F_i(0)| <= 1e-3 S_i AND
<= 0.1 N or 0.1 Nm, where S_i = sum_j a_j |H_i(omega_j,heading_j)| using
dimensional native carrier excitation. Report degenerate S_i explicitly.
Body displacement/velocity and radiation-memory initial conditions remain zero.
At least one incoming packet must reach 25% of canonical amplitude at the hull
within ten seconds with nontrivial motion; retain initial/peak ratios. These
are startup gates, not motion-accuracy or convergence claims. Independent
reasoning and padding checks: numerical_methods_researcher/startup-review-2026-10-24.md.

First obtain minimum-mesh end-to-end motion, then characterize changes in all
six histories with refinement. Do not replace that milestone with whole-library
audits, coefficient gates or optimization. Independently reproduce reported
engine defects on the selected version; resolve confirmed result-invalidating
issues or withhold the affected result. Do not blindly carry old workarounds.

Only resulting six time histories determine numerical motion variation and
reference uncertainty. Intermediate A/B/Ainf/coupling differences are diagnostics,
not standalone0.1% acceptance gates. Missing/corrupt data and invalid solves do
not establish a result. Preserve old failures under original policies.

Historical environment setup (not numerical validation): Capytaine3.0.0;
MATLABR2025b with installed/licensed
Simulink/Simscape/Multibody; upstream WEC-Simv7.1.0 commit
a23d2993807410fe78cc63546c1ba9d81912024f. Library/importer/class loading checked,
NOT an actual boat simulation. See [setup evidence](../capytaine-wecsim-environment.json).
The maintained fork pin and validation status are in [wecsim-fork.json](wecsim-fork.json).
Broader Topwater comparison/runtime goals and major speed campaigns remain held
until the required physical comparison; local residuals, energy checks, fitted
loads and ODE tolerances do not establish motion accuracy.

## REF-002: buoy-derived strong-motion canonical case

Latest regroup: stop unproven bug investigations and deliver ONE less energetic
collected-buoy-derived60s case and video with true locally calm incoming-wave
startup. No official PTO-zero or further mesh job is allocated. Use a genuinely
milder observation, not a scaled old realization. Quantify resulting all6motions
and applicability; do not tune amplitudes/phases to force a peak criterion.
The close camera view must not limit the physical incident-wave construction.
WEC-Sim here integrates body motion using open-water hydrodynamics, not a finite
water-domain simulation. Displayed incident field is prescribed, not a computed
scattered/radiated fluid field.

Keep immutable generator parent amplitudes/directions/phases and source metadata.
A derived physical finite incoming train necessarily broadens the stationary
spectrum; report its actual support, energy/directional changes and any removed
unsupported components. Do not claim exact stationary Hm0/bin closure for that
finite train unless actually established. No global temporal force/elevation
ramp, history seed or post-hoc amplitude/phase fitting. Native ramp remains zero.

Latest user steering,2026-10-24: final validation duration is60s, with all startup
retained. Prune unsupported incident frequencies using the actual physical model
and defensible accuracy evidence, not an arbitrary WEC-Sim cutoff. Capytaine's
full radiation-diffraction excitation is FK+diffraction, not an FK-only model.
The existing wavelength>=2beam/.6Hz project support envelope remains; distinguish
it from numerical mesh/time/IRF coverage. Pruning incident components does not
repair an unreliable radiation-memory operator.

Final startup: body initially at rest in TRUE locally calm water, incoming waves
initially upstream, no prior body motion, zero velocity and zero radiation memory.
Body rest in an established stationary sea is insufficient. Declare and verify
full-hull initial surface/fluid-velocity tails and all six initial excitation
loads. REF001's original tiny-wave0.1N/Nm cap remains immutable for that result.
For the new REF002 case, the coordinator declares an engineering calm-load
tolerance epsilon=.001: initial force-vector norm<=epsilon*m*g and moment-vector
norm<=epsilon*m*g*Lchar (Lchar=max full-hull length/beam), with each component also
<=epsilon*S_i in the SAME fixed native basis as its carrier-load scale. Never
rotate a positive scale vector as though it were a physical wrench. Use the
original full-hull eta/u numerical-calm bounds with the new carrier scales.
This is a scale-based engineering assumption, not a universal precision or
linearity proof; report actual residuals. Preserve all old rejected attempts
under their original0.1N/Nm policy. No wave/force/history changes follow from
assessing an immutable train under the new explicitly recorded acceptance policy.
No hidden temporal force ramp. Declare incident field and initial pose consistently.
No stationary-history q/v/CIC seed in that final startup. Viscous damping remains
off unless separately justified and authorized for matching physics in both models.
No damping fit, force cancellation, restraints or prescribed post0 trajectory.

Use the same physical boat/mass/CG/inertia, density/gravity/depth and frames.
Existing collected-buoy generator/catalog is the source. Numerical pruning may
remove supported whole measured frequency bands; retain their original amplitudes,
MEM directional spreading and phases without renormalization. Report original and
retained variance/Hm0/period/direction, mask, rationale and source identities.
Final case is ONE frozen state/seed/heading, not a suite. Target roll about5-10deg
plus appreciable other5DOF. A peak target is not sustained-response statistics.

Default bow faces dominant arriving energy; any chosen offset must be explicit.
Keep global source bank, retained subset, nauticalFROM/globalpropagation/bodyframe
conventions and R/t mappings distinct. Never rotate only the rendered box.

Assume a linear operating range a priori from exact full-hull buoyancy/righting,
restoring/waterplane/wet-set curves against fixedmean geometry. Report combined
heave/trim and contact/flooding. These justify an engineering assumption, not
nonlinear dynamic proof; no nonlinear trajectory model/adapter is requested.
Only six complete motion histories establish motion variation/accuracy; coefficient
and energy diagnostics guide rejection/investigation, not arbitrary coefficient
thresholds. No actual radiation implementation bug has been established.

Old full-bank and nonzero stationary-history prototypes retain their original
policies and evidence; they are superseded for the final requested input/startup.
Topwater solver implementation/tests/simulation/optimization remain HELD. Existing
generator reference-input work is allowed; externalreference loads/trajectories
never replace independently computed Topwater physics.
