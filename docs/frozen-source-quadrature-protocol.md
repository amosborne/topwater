# Frozen-source pressure integration: prospective diagnostic

Scope: separate pressure-output quadrature effects from the remaining measured
body-refinement increment. Hold the paired study's864-panel Capytaine normalized
indirect sources fixed. No new BEM solves, geometry/model change, runtime forcing,
remeshing or speed optimization. The previous failed gates remain failed.

Use the pinned complete paired-body evidence bundle. Verify its study JSON/NPZ,
actual producer graph, physical-input bytes, mode order and source normalization.
Use Capytaine3.0, the same infinite-frequency Delhommeau settings, world CG, full
dry inertia, density1025, gravity9.81 and infinite depth. Topwater direct potential
traces are not source densities and must not enter this calculation.

First evaluate the frozen sources at native Capytaine face centroids and reproduce
all six potential traces and all36 native load entries at1e-10 scaled tolerance.
Then subdivide each pressure-integration triangle into4/16/64 equal-area children,
leaving the source mesh/densities unchanged. Evaluate only strictly interior
descendant centroids on their owned physical hull triangle, below z0. Single-layer
potential has a continuous body trace; no interior-volume continuation or gradient
jump is used. Verify ownership, area sums and finite values before native work.

Pressure matrix A_Q=rho*sum(w*N_fluid*chi), with N_fluid=[-n_body,
(x-CG) cross(-n_body)]. Use actual Capytaine coarse-face areas divided by the
declared child count, and its unchanged facet normals. Record complex raw matrices;
their imaginary part must fit1e-10 scaled before reporting real Ainf entries.
Retain every target, owner, weight and all six evaluated potentials.

Apply the original0.001 relative +1e-6 normalized policy to the two finest pressure-
quadrature increments and their trend as a stabilization DIAGNOSTIC only. Preserve
all36 signed differences, absolute-entry magnitudes, ratios and failure masks.
No quadrature stabilization can qualify source discretization or the full reference.

For each quadrature level preserve exact telescoping:
A_f-A_c=(A_Q-A_c)+(A_f-A_Q), where A_c is native864-panel and A_f native3456-panel.
Call the last term an unresolved remainder: it includes changed source/collocation
solution and finest-output quadrature error, not uniquely pure source error.
Do not infer a continuum bound from cancellation or measured increment accounting.

Archive actual producers, controls, pinned parent and physical inputs before
native evaluation. One native thread, blocks<=100points by864sources, at most
180 wall seconds including preparation; check before/after each native block and
at final completion. Required targets:864 native +3456/13824/55296 descendants.
Missing cases/dependencies, invalid points or failed centroid reproduction cannot
pass. Retain a concrete failure; do not silently skip the expensive level.

Deliverables: complete diagnostic or rejection, focused ownership/weight/phase/
pressure-gate tests, independent adversarial review, immutable durable evidence,
issue #30 update and fresh handoff. Stop after interpreting this diagnosis. Full
reference, finite-wave/initialization, nonlinear RHS, motion and runtime remain
unqualified regardless of its numerical outcome.
