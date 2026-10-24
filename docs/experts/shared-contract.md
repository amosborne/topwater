# Shared physical and producer/consumer contract

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

Use the same physical boat/mass/CG/inertia, density/gravity/depth and axis
conventions. Generate a credible new sea from topwater.sea_states and the
collected buoy catalog, with full observed spectrum, MEM directional spread,
seed and component phases retained. Select a genuinely higher-energy record;
do not multiply the old packet or fit phases to motion. Internal candidate
selection is permitted; final deliverable is ONE frozen case, not a suite.
Target approximately5-10deg roll plus appreciable motion in the other5DOF.

Default forward orientation faces the dominant arriving direction. Export both
buoy FROM and world propagation conventions, measured/realized directional
first moments/resultants, chosen body heading and explicit heading offsets.
Do not conflate rotating the rendered box with rotating hydrodynamic modes.
Preserve the generator bank in its global frame and document body-frame mapping.

Assume a linear operating range a priori from exact full-hull buoyancy/righting,
restoring, waterplane and material-wet-set curves relative to WEC-Sim's fixed
mean geometry. Report heave/trim coupling, contact/flooding and sign failures.
These curves justify an engineering assumption; they do not prove nonlinear
dynamic accuracy. No new nonlinear trajectory model/adapter is requested.
Do not change the boat or add damping/restraints to force the target response.

Before dynamic launch, declare the new sea, heading, duration/initial state and
producer frequency/direction coverage. A stationary generator realization is
not calm water or an initialized nonlinear fluid field. Choose and document
canonical initialization consistently for the future comparison; do not silently
reuse the old packet, impose hidden force ramps, drop directional energy or
claim body/sea initialization not represented by the chosen model.

The supported physical envelope remains incident wavelength at least twice beam
and frequency at most0.6Hz; enforce the existing generator's whole-band policy.
For9ft beam andg9.81 the wavelength condition is the tighter upper bound near
0.53346Hz. This limit governs physical input, not the hydrodynamic integration
band needed by the radiation/excitation model. Breaking/slamming/flooding/loss
of contact remain outside supported physical scope.
