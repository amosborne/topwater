# Pressure discrepancy attribution: bounded diagnostic

The first canonical closed-boundary pressure probe has large all-six added-mass
discrepancies against the retained finest Capytaine reference. Before changing
resolution or optimizing, separate error owners using exactly the same wet
physical triangles and world CG. No comparison in this diagnostic is convergence
or promotion; the reference's original failed gates remain failed.

Archive the completed pressure probe, complete parent evidence, actual producers,
calm physical inputs and the farthest-first body triangles before native solves.
Prove those triangles preserve the canonical wet facets, area moments and edges.
Run serially, one native thread, <=60 s checked between solves:

1. Existing Topwater instantaneous_added_mass on those triangles with the exact
   flat Dirichlet halfspace and infinite depth (six auxiliary columns together).
2. Capytaine infinite-frequency inertia on the identical farthest-first triangles
   (six solves). Preserve all 36 raw entries and the actual Green settings.
3. Compare each with the closed-boundary probe and retained fine-reference raw
   matrix. Retain the fine-reference mesh envelope separately. Normalize each
   entry by sqrt(Mii*Mjj); report all entries, not just diagonal agreement.

Use the already-declared candidate 1% relative +1e-5 normalized entrywise policy
only as a diagnostic comparison, explicitly never a physical qualification pass.
Same-mesh correlation is not convergence. A large closed/halfspace difference
implicates finite-surface/domain representation; a large same-body Capytaine/fine
difference implicates body discretization; disagreement between the two same-body
solvers can additionally expose formulation/discretization error. Do not assign
all error to one owner from a single difference or silently relax a gate.

Exclusions: new formulation, body remeshing, domain refinement, speed optimization,
nonlinear RHS or motion qualification. Stop after the attribution or a concrete
failure, then independent review and the current checkpoint/handoff. Preserve all
failed gates and incomplete baseline status.
