# Paired body-only refinement: prospective protocol

Scope: assess body discretization on the unchanged canonical calm-water surface
using the existing Topwater flat-halfspace inertia path and Capytaine infinite-
frequency radiation. This is a prerequisite investigation, not the complete
finite-wave reference, nonlinear RHS, time-domain validation or speed benchmark.

The previous 48-panel seed has large numerical odd-parity couplings. Declare one
symmetry correction before solves: clip world y>=0, canonicalize only roundoff
coordinates within1e-12 m, then reflect that half with reversed winding. Verify
the original physical facets, area moments and conforming edges; reject any lost
facet or physical asymmetry. Do not average matrices, force coefficients to zero,
fit loads or select a better half after viewing results. Symmetry is a verified
property of this reference fixture, not an assumption in the runtime solver.

Construct four nested levels by0/1/2/3 midpoint refinements of the same frozen
half seed, reflecting each level. Farthest-first vertex order preserves orientation
and supplies Capytaine's intended maximum radius. Record actual seed/level counts,
quality and proofs in a geometry artifact BEFORE hydrodynamic solves. Preflight
all levels against unchanged12 GiB Capytaine workspace estimation and512 MiB
Topwater S/K storage. No silent finest-level omission. Archive actual producers,
the pinned complete parent bundle, calm body/pose/full dry inertia, density1025,
gravity9.81, infinite depth and zero speed. Serial jobs, one native thread,
<=180 s for the complete study, checked between solves and after each stage.

At each level retain all36 raw entries for both solvers and all six auxiliary
potential fields. Normalize Capytaine's symbolic infinity-frequency traces by
(-i*omega), explicitly canceling the symbolic factor; preserve normalized source,
potential and prescribed normal data. Reapply original S/K before their cache
becomes LU. Require field/operator and independent load reconstruction at1e-10
scaled tolerance; this is internal consistency, not physical accuracy.

Use the original entrywise normalization sqrt(Mdry_ii*Mdry_jj). Capytaine reference
acceptance is unchanged: both finest increments <=0.001*abs(Afinest)+1e-6, and
latest <=previous+1e-6. U=max(abs(A4-A3),abs(A3-A2)), empirical not certified.
Topwater's finest mesh increment must fit0.01*abs(Cfinest)+1e-5, and
abs(Cfinest-Afinest)+U must fit0.01*abs(Afinest)+1e-5 for every entry. Candidate
agreement cannot qualify while the reference gate fails. Same-mesh correlation
is not convergence. Record every failed entry and all raw cross-couplings;
no Richardson extrapolation without evidence for its assumptions.

Compare with the previous finest reference as a diagnostic only; it remains
unqualified. A successful Ainf-only gate does not qualify finite A/B, excitation,
phase, periodic fluid fields, initial compatibility or time-domain motion.

Deliverables: reviewed geometry freeze, complete paired study or concrete rejected
family, focused symmetry/nesting/gate tests, independent adversarial audit,
durable producer/input/result evidence, issue #30 update and fresh handoff.
Stop after interpreting the full paired study. No free-surface/domain remeshing,
new formulation, optimization or one-second integration in this checkpoint.
