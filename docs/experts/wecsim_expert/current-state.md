# Current state: WEC-Sim / BEMIO expert

Setup, fork creation and source-map/reproduction planning are complete.
Reference work is queued, not authorized; follow the
[handoff](../../checkpoint-handoff.md). WEC-001 scope and acceptance are in
[assignments](../assignments.md); numerical scheduling follows coordination.md.
No engine corrections or numerical tests have occurred. No motion/video exists.

Maintained fork: https://github.com/amosborne/WEC-Sim, local
artifacts/dependencies/topwater-wec-sim, branch topwater/reference-v7.1.0,
published at official v7.1.0 commit a23d2993807410fe78cc63546c1ba9d81912024f.
See [fork pin](../wecsim-fork.json); recheck checkout identity before relying on it.
The historical setup note records its earlier local-only state.

Marine-gyro's recorded engine commit was independently matched to clean checkout
0753b2e47f2457c078751dcfe5d251d1767b80ab. Yaw order, unequal postprocessing grids,
angular-bin normalization and row/column directions remain unreproduced leads
on selected v7.1.0. No workaround is accepted without independent reproduction.

BEMIO supports a separate hydrostatics directory as its third argument (suffix
with leading separator); joint executed producer/consumer signoff is pending.
Open first-case questions: full inertia/frame/Ainf mapping, per-direction imports,
excitation-IRF noncausality and actual convolution prehistory/padding, upstream
arrival shared by forcing and display, and explicit simu.rampTime=0. Keep the
first ten seconds from actual rest. v7.1.0 disallows elevationImport with
nonlinearHydro=1; initial case is linear, with qualitative imported surface display.
Preserve the unfinished exporter; its packet/export controls are not approved.

[Setup/source map and reproduction plans](native-setup-2026-10-24.md),
[read manifest](native-setup-read-manifest.json),
[historical leads](marine-gyro-leads.md) and
docs/capytaine-wecsim-environment.json contain evidence and limits. Recorded
MATLAB/engine setup is historical, not live readiness or an execution lease.
Next: await user-authorized reference work and a bounded coordinator assignment;
independent review is required before relying on physical-result corrections.
