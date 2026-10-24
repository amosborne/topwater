# Shared physical and producer/consumer contract

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
