# Accepted minimum-case incident input and inertia boundary

This checkpoint validates reference inputs and two exercised engine boundaries;
it does **not** establish motion or video. Body execution remains held for the
newly identified imported-excitation mask dispatch integration issue.

## Frozen incoming field

Canonical carriers/amplitudes/phases/headings are unchanged:0.1/0.25/0.5Hz,
0.02/0.01/0.005m,0/pi/2/pi,0/45/90degrees. Each carrier has a normalized compact
cosine-fourth spectral envelope with half-width0.3/1.47/0.18rad/s and phase
`canonical_phase + (omega-carrier_omega)*focus_time`. Focus times100/11.5/200s
place the field upstream at t0. Every spectral component obeys exact one-way
deep-water dispersion; no time ramp, force switch, clipping or overwritten
sample is used. Actual frequency support0.0160422..0.5286479Hz lies within native
hydrodynamic support and the minimum-wavelength physical bound0.5334587Hz.
The later low/high-frequency packets remain in the frozen input with canonical
amplitudes; their peaks are not claimed to arrive within this first10seconds.

The first leading focus10s candidate failed both direct and retained initial
heave/pitch force gates, preserved as packet-attempt2. Fixed-order cheap screening
of11/11.5/12/12.5s chose the first passing11.5s candidate; all four rows and
per-component contributions are retained in
`artifacts/reference/minimum-20261024/packet-screen1/result.json`.

Frozen input arrays/files: `artifacts/reference/minimum-20261024/packet-attempt3`.
The packet.npz SHA256 is
`8925a52abf598eac55c12662cb3a7f53cc1d3bb175169201b47b5000cb5428ea`.
Its scalar/input gates pass; the supplemental whole mean-world hull planform
times wet-depth:0 coverage also passes,587412grid samples at spacing<=0.025m.
Conservative sample/gradient/global-Hessian bounds are3.323565um elevation and
25.495529um/s fluid vector speed, below35um and43.982297um/s prospective bounds.
The world region is x[-2.824982,4.490187],y[-1.3716,1.3716],z[-0.402142,0]m.

Absolute initial retained-IRF forces/moments are
`[.011597,.013599,.049249,.009122,.056153,.000928]`, and independent native-H
spectral values are `[.011856,.013017,.074297,.008029,.056100,.000988]`.
Both pass each axis's unchanged `min(0.1,1e-3*carrier_force_scale)` bound.
Leading packet elevation reaches59.4093% of its canonical amplitude within10s.
Finite+-20s excitation-kernel versus direct spectral maximum discrepancy over
the first10s is `[.000675,.000885,.028646,.001702,.004532,.0000909]`N/Nm;
this is reported as a diagnostic, not hidden or treated as motion accuracy.

Independent numerical_methods_researcher read-only source/manifest/array/region
review passed final incident-input startup. Evidence: packet-attempt3 and
packet-full-planform diagnostic.json, with exact archived source/manifests in
their corresponding `-source` directories. All accepted wave and spectral files
remain immutable for first-body execution. Exact body rest and actually applied
force still require native runtime evidence.

## Corrected native solid boundary

Upstream/internal WEC-Sim products are[Ixy,Ixz,Iyz]. R2025b
[File Solid](https://www.mathworks.com/help/releases/R2025b/sm/ref/filesolid.html)
expects[Iyz,Izx,Ixy], using tensor off-diagonals (already negative integrals).
Actual official H5/bodyClass hydroForcePre/adjustMassMatrix plus actual library
parameter evaluation reproduced errors0 for canonical dry,6 for distinct dry
products,1071.643398kgm2 for native full dry+added mass. Original nonsymmetric
added-mass residuals were retained, never blanket symmetrized.

Proper Simulink API correction reordered exactly three File Solid boundary
expressions with[3 2 1]. Public/internal properties and residual matrices remain
unchanged. Fixed actual matrix roundtrip errors are0/0/1.819e-12. Independent
capytaine_expert source/math/executed/semantic review passed selected rigid path.
Maintained fork project branch now published at
`fcd7c365a1201bbb56ce7f571a5baf9b2a94d89b`.

Saving via R2025b upgraded the original R2020b library serialization. Review
independently verified146system connectivities and39embedded MATLAB scripts
unchanged, all531migrated port counts equal; equivalent empty-choice migration
and an inactive Flex sensor default were serialized. The active rigid sensor
did not change. All variants are not compiled/physically validated. Evidence:
`artifacts/reference-inertia-stock-attempt2`, `reference-inertia-fixed` and exact
source snapshots/manifests; library-semantic-diff.json documents the save changes.

## Remaining integration gate

Native excitation parent mask creates one force branch per wave object, while
each elevationImport branch reads the same precomputed body total. Source
inspection predicts3times total for3imports. Actual mask reproduction is not yet
complete because native cold-start model authoring required initialization-order
fixes. No body simulation has been launched. Independently reproduce/correct
dispatch, handle reinitialization/stale clones and mixed-type groups explicitly,
then compare native applied excitation against all six accepted retained histories.
