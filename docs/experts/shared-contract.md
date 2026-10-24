# Shared physical and producer/consumer contract

Current scope: standalone Capytaine -> official BEMIO -> WEC-Sim external
reference. Ten-second six-DOF motion, isometric video, then variation of those
motions with solver mesh. Topwater solver implementation is held.

Authoritative neutral physical input:
[canonical-boat.json](../../reference/capytaine/canonical-boat.json).
It carries the minimum full24ft hull and frozen17-triangle wet surface, mean
equilibrium pose/world CG, full dry mass/inertia including products, canonical
three-component mild sea and provenance hashes. Density1025kg/m3, gravity9.81m/s2,
infinite depth, zero forward speed. Solver remeshing must preserve this physical
surface and physical properties. The retired canoe is not this boat.

World axes: x forward, y port, z up; right-handed. Modes Surge/Sway/Heave/Roll/
Pitch/Yaw, rotations about mean WORLD CG. Capytaine complex convention
Re(amplitude*exp(-i*omega*t)); headings are propagation radians CCW from+x.
Capytaine excitation is diffraction PLUS Froude-Krylov per unit elevation.
WEC-Sim directions are degrees, with its own phase convention and WAMIT
normalization. Its official readCAPYTAINE negates excitation imaginary parts;
verify the interface once, never manually conjugate a second time.

Producer, accountable capytaine_expert: native netCDF, native hydrostatic files,
physical body/mesh/CG/inertia/DOF identities, exact frequencies/headings, versions,
code/input/source hashes and reproducible export commands. Keep input geometry
distinct from derived meshes. Existing export_case.py is an unfinished/unrun
prototype with undeclared packet/export controls; it is NOT accepted data.

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

Known integration questions, not yet resolved physical results: per-direction
imported components; excitation IRF noncausality relative to its wave reference
point; padding must exist in actual convolution timeline (extra CSV rows can be
discarded); explicit simu.rampTime=0 instead of default100; preserve the first
ten seconds from actual rest, never a later already-moving cropped window.
WEC-Sim7.1 disallows elevationImport with nonlinearHydro=1; label the initial
case linear hydrodynamics. Imported surface display is qualitative. Check missing
hydrostatic-file zero defaults, full inertia products and frame signs.

Only resulting six time histories determine numerical motion variation and
reference uncertainty. Intermediate A/B/Ainf/coupling differences are diagnostics,
not standalone0.1% acceptance gates. Missing/corrupt data and invalid solves do
not establish a result. Preserve old failures under original policies.

Pinned known-working setup: Capytaine3.0.0; MATLABR2025b with installed/licensed
Simulink/Simscape/Multibody; upstream WEC-Simv7.1.0 commit
a23d2993807410fe78cc63546c1ba9d81912024f. Library/importer/class loading checked,
NOT an actual boat simulation. See [setup evidence](../capytaine-wecsim-environment.json).
The maintained fork and its tested commit are owned/recorded by wecsim_expert.
