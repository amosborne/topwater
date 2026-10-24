# Current state: WEC-Sim / BEMIO expert

Minimum-mesh reference work is authorized by the 2026-10-24 user goal; follow
[handoff](../../checkpoint-handoff.md). WEC-001 scope and acceptance are in
[assignments](../assignments.md); numerical scheduling follows coordination.md.
Stock imported-excitation cropping is independently reproduced on selected
v7.1.0: actual wave setup discards past/future input before noncausal convolution.
Evidence: artifacts/reference-imported-excitation-stock/result.json and archived
source/manifest in artifacts/reference-imported-excitation-stock-source.
Scalar independent expected/actual maximum errors are [2,16,0] for future,
past and interior impulses. The allocated job terminated exit0; lease released.
The retained-history and component-sum corrections passed independent executed
review and are committed/pushed on the maintained project branch at
2f43f547301b223cd9a07e84dc88c9375249af0e. Exact tests, failures and scope are in
[the imported-excitation checkpoint](imported-excitation-2026-10-24.md).
No motion/video exists. Accepted input and inertia facts are in the
[accepted-input checkpoint](accepted-input-2026-10-24.md).

Maintained fork: https://github.com/amosborne/WEC-Sim, local
artifacts/dependencies/topwater-wec-sim, branch topwater/reference-v7.1.0,
based on official v7.1.0 commit a23d2993807410fe78cc63546c1ba9d81912024f.
See [fork pin](../wecsim-fork.json); recheck checkout identity before relying on it.
The historical setup note records its earlier local-only state.

Marine-gyro's recorded engine commit was independently matched to clean checkout
0753b2e47f2457c078751dcfe5d251d1767b80ab. Yaw order, unequal postprocessing grids,
angular-bin normalization and row/column directions remain unreproduced leads
on selected v7.1.0. No workaround is accepted without independent reproduction.

BEMIO supports a separate hydrostatics directory as its third argument (suffix
with leading separator). Joint producer interface is approved for provisional
0.1:0.1:12 rad/s, 0/45/90-degree headings; executed consumer signoff is pending.
Native corrected files under artifacts/reference/minimum-20261024/native-gl2
passed independent physical review. reference/wecsim/prepare_bemio.m executed
official import/write with exact native MAT metadata and phase checks; H5 and
interface-check.json are in artifacts/reference/minimum-20261024/consumer.
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
Incoming packet startup passed independent full-physical-hull review; the exact
packet-attempt3 spectral/input files remain frozen. Inertia boundary correction
is independently reviewed and published at fcd7c365a1201bbb56ce7f571a5baf9b2a94d89b;
see fork pin for R2025b serialization limits. Native excitation parent mask may
dispatch the precomputed imported total once per wave object; actual-mask
reproduction/correction and stale-clone/mixed-group policy now block body run.
Cold-start model authoring is in progress, no simulated warmup or trajectory.
Next: complete actual-mask gate, then exact native10s case/force/frame checks and
video. Four historical marine leads
remain unconfirmed on selected release; passive-yaw/full-direction spectrum
paths are unexercised by the proposed linear imported-packet case. Equal-grid
fixed ode4 output is proposed; output-grid behavior still requires verification.
No historical workaround has been carried over. Independent review is required
before relying on physical-result corrections.
