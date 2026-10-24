# Common-output source-resolution diagnostic: prospective protocol

Scope: compare archived 54/216/864-panel infinite-frequency Capytaine source
solutions for the SAME canonical calm physical surface. This is a three-source
diagnostic, not a replacement or qualification of the original failed four-level
study or finite-frequency reference. No BEM solve, remeshing, formulation change,
optimization campaign or motion run. Full/source qualification stay false.

Use ONE common pressure-integration grid: retained 864-source Gauss12/16/20
targets, actual facet fluid normals and weights. Regenerate exact targets/owners/
weights and compare against the SHA-pinned order-20 continuation. Verify its
complete archived producer graph, manifest/NPZ bytes in that bundle, physical/
source/Green identities and required integration pass. Pin its historical
Gauss001 bundle as well. Archive the actual input bundles/code/controls before
native evaluation. Historical producer code stays separately identified.

Verify each source mesh's per-facet containment/moments, conformity and complete
physical surface against the same frozen parent. Repeat its own native centroid
potential and all36-load reproduction for54/216 at1e-10 scaled tolerance. Reuse
864 fields only from the bound continuation and reconstruct its loads. Targets
and normals never change when the source resolution changes. The shared bounded
full S*sigma evaluator independently evaluates both reflected target halves;
no mirrored fields, fitted loads, zero couplings or matrix symmetrization.

For EACH source, retain all six complex fields at all common target orders and
all36 pressure matrices. The existing empirical integration policy applies:
both12→16 and16→20 increments <=0.001*abs(A20_normalized)+1e-6, latest increment
<=previous+1e-6. Entrywise normalization is sqrt(Mdry_ii*Mdry_jj). Retain physical
quadrature envelope Uq=max(abs(A20-A16),abs(A16-A12)), not a continuum bound.
Nonfinite values, nonzero imaginary matrices beyond1e-10 scaled, missing orders,
changed targets/inputs or failed reproduction cannot pass.

At source20-order outputs Q0,Q1,Q2 define Dprev=abs(Q1-Q0), Dlast=abs(Q2-Q1),
Eprev=Uq0+Uq1, Elast=Uq1+Uq2, all normalized entrywise as above. Preserve signed
source differences and separate absolute/envelope/charge arrays. Source-policy
diagnostic requires Dprev+Eprev AND Dlast+Elast <=0.001*abs(Q2)+1e-6 for ALL36
entries. For a quadrature-robust trend require
Dlast+Elast <=max(Dprev-Eprev,0)+1e-6. Report this explicitly stronger uncertainty
check separately from the original uncharged native-study trend. Any failed or
missing integration control prevents a source-policy pass. Preserve source
increment envelope Us=max(Dprev,Dlast) separately from integration uncertainty.
No Richardson, assumed spectral rate or certified continuum-bound claim.
Even a diagnostic pass leaves the fourth source level and finite waves required.

Run serially, one native thread, <=100 targets per block. Newly evaluated fields
are2*(124416+221184+345600)=1382400 points plus54+216 native controls. Projected
native cost195 s from measured pair counts, a screen not a timing claim. Declare
360 wall seconds INCLUDING preparation, check budget before/after each native
block and final completion. Exact process must be terminal before source/input
changes or another numerical job. Reference costs are offline, never Topwater
simulation speed. Separately preflight3456 sources before any future launch;
order20 alone on this fixed grid projects1247 s plus required controls.

Deliver complete diagnostic or concrete rejection, raw fields and all36 signed/
absolute source increments with separate charged Uq/Us and failure masks/ratios,
actual immutable producer/input graphs, focused gate/evaluator tests, independent
prelaunch/postrun review, reproducible durable evidence and tracker/handoff.
Commit/push the checkpoint and stop after selecting the next source-control step.
