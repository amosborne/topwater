# Imported excitation and native interface checkpoint

Selected upstream v7.1.0 base `a23d2993807410fe78cc63546c1ba9d81912024f`.
Maintained fork correction `2f43f547301b223cd9a07e84dc88c9375249af0e`, published
only on user-owned `topwater/reference-v7.1.0`. No body simulation/video yet.

## Independently reproduced defects

1. Stock waveClass setup resamples imported elevation only on simulation times;
   bodyClass then convolves this cropped data with a noncausal kernel. Future
   impulse at11s and past impulse at-1s are omitted from requested0:10s force.
   Interior impulse at5s establishes indexing. Independent direct scalar sums
   versus stock maximum errors `[2,16,0]` are in
   `artifacts/reference-imported-excitation-stock/result.json`; exact source and
   manifest are in `artifacts/reference-imported-excitation-stock-source`.
2. Stock hydroForcePre uses only waveObj(1) for elevationImport. Two known
   imported components should sum3K but stock givesK, maxerror32. Executed actual
   preprocessing with radiation stubbed solely to isolate excitation, using
   preserved upstream source. Evidence:
   `artifacts/reference-multiple-imports-stock-attempt5/result.json` and
   `artifacts/reference-multiple-imports-stock-source/manifest-attempt5.json`.

Both corrupt first-case forces if left unresolved. The correction retains raw
file history, integrates at explicit physical `t-tau`, rejects insufficient
coverage, and sums each component using its own heading. It leaves output-grid
elevation intact; no body prerun, force switch or packet overwrite is used.

## Fixed execution and independent review

`reference/wecsim/check_imported_excitation.m` executes independent scalar sums
for dt1/.5, both edges and interior, ramp/output-grid invariance, fail-closed
history coverage, distinct0/90-degree unequal kernels/amplitudes through actual
hydroForcePre across all6DOFs, and asymmetric[-1,2] lag support. Final errors are
all0. Evidence: `artifacts/reference-imported-excitation-fixed-attempt2/result.json`,
exact source/hash manifest in `artifacts/reference-imported-excitation-fixed-source-attempt2`.
Independent numerical_methods_researcher source plus executed review passed for
this exercised linear imported-excitation path. This is not body-dynamics validation.

Preserved failures: engine discovery/recovery attempts are bootstrap-only;
stockmultiple attempt4 cleared its root path variable before test; fixed attempt1
exposed missing per-component matrix allocation in the actual production path.
Their terminal logs remain under artifacts and are not relabelled successes.

## Official native interface

`reference/wecsim/prepare_bemio.m` executed official readCAPYTAINE, radiationIRF,
excitationIRF and writeBEMIOH5 on independently reviewed `native-gl2` exports.
Exact native MAT CG/CB/volume/K supersede rounded legacy text metadata; K is
normalized by rho*g, direct infinite-frequency added mass by rho after IRF
generation. Independent named-DOF FK+diffraction comparison checks the importer's
single imaginary sign conversion. Both real and imaginary residual maxima are
1.4552e-11 dimensional. Evidence:
`artifacts/reference/minimum-20261024/consumer/interface-check.json` and
`artifacts/reference-bemio-source/manifest.json`. Excitation lag+-20s/radiation10s
are prospective for the ten-second packet case, pending startup diagnostics.

## Remaining exercised-path questions

No startup packet has passed the prospective eta/fluid-velocity/force calm gates.
No compiled physical inertia mapping has been accepted. A selected-version lead
is that WEC-Sim's public products are [Ixy,Ixz,Iyz] but the actual Simscape block
expects [Iyz,Izx,Ixy]. The
[R2025b Inertia documentation](https://www.mathworks.com/help/releases/R2025b/sm/ref/inertia.html)
defines these as off-diagonal tensor entries (already negative integrals).
This needs an independent actual-path reproduction before any correction/reliance;
the canonical dry tensor's sole middle product does not establish generic coverage.

Four historical marine-gyro leads remain unconfirmed on selected release:
yaw mapping, unequal postprocessing grids, full-direction angular-bin measure,
and direction row/column shape. Passive-yaw and full-direction spectrum paths
are not selected for the linear imported-packet case. Fixed ode4 equal grids are
proposed; actual sensor/output mapping still requires case evidence. No old
workaround has been copied and unexercised paths are not claimed unaffected.
