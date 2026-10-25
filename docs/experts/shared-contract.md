# Shared physical and producer/consumer contract

## Canonical NPZ workflow extension, 2026-10-25

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
