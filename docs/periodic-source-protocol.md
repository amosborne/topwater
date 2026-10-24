# Periodic source recovery: bounded checkpoint

Scope: recover the indirect Capytaine source densities and hull potential traces
for the already frozen mild periodic case on all four immutable mean-hull meshes.
This is a necessary reference-data dependency, independent of the eventual fluid
grid. It does not select a new Topwater formulation or qualify a restart.

Before launch, archive the complete producers, this protocol's JSON controls,
the mild case, and the pinned raw-reference bundle. Run serially, single native
thread, at most 1800 seconds checked between solves, with the existing 12 GiB
estimated workspace ceiling. There are exactly 84 solves: four meshes, three
frequencies, six unit-displacement radiation modes and one selected diffraction
heading at each frequency. Preserve every basis source, potential, pressure and
Neumann boundary trace, not just the combined state.
Reapply the exact cached indirect operators to every source basis: S*sigma must
reproduce its saved potential and K*sigma its prescribed normal velocity at the
trace tolerances below. Preserve these produced traces separately from input
boundary data. Apply the elapsed budget after Green construction, after every
solve, after operator verification and before accepting completion as well.

For each original complex load component and all 36 A/B entries, require
abs(recovered-original) <= 1e-12*physical_scale + 1e-11*abs(original).
The A scale is sqrt(M_ii M_jj), B scale is omega times that, force scale is
rho*g*[L^2,L^2,L^2,L^3,L^3,L^3]. This is a reproducibility check, not a new
continuum acceptance tolerance. Require p=i*rho*omega*phi for every saved trace
at relative/absolute-normalized tolerance 1e-12. Reintegrate pressure independently
through the mesh normals, areas and rigid-mode vectors, including moment arms.

Use exp(-i*omega*t), zeta=a*exp(i*phase), xi=RAO*zeta. Compose scattered source
sigma=zeta*sigma_diff+sum(xi_j*sigma_rad,j). Incident potential is analytic and
is not a hull source. Total hull potential includes it. Require the resulting
total normal-velocity trace to equal -i*omega*sum(xi_j*mode_j dot n), and the
integrated total pressure to match the retained radiation plus excitation loads.
All source units refer to Capytaine's pinned indirect Green-function convention;
they cannot be passed to Topwater's Rankine operator as if they were compatible.

Exclusions: evaluation on a proposed nonlinear grid, continuation inside the
mean hull/above the mean free surface, a nonlinear geometry map, field spatial
uncertainty, compatibility projection, startup transients, runtime force driving
and speed optimization. No sampled eta0/phi0 restart is produced by this stage.
The four missing mild-case readiness flags remain false. The failed original
coefficient-convergence decision remains failed.

Deliverables: immutable basis/composed source data with all mesh geometry and
Green settings, raw-load reproduction and phase/normal checks, focused tests,
independent adversarial review, durable bundle, issue #30 update and fresh
handoff. Stop after these checks or a concrete failed recovery. Next work remains
an admissible nonlinear geometry map/grid and field uncertainty; an actual
coupled setup/RHS screen must reject routes projected above 200x before a long run.
