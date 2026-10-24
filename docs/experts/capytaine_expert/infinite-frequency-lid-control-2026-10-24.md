# Infinite-frequency API audit and isolated lid proposal

Audit outcome: no producer API misuse identified. Control source independently
reviewed and allocated24 calls completed once, terminal0; executed result below.
Earlier sources and raw
outputs remain unchanged; this is not a coefficient correction.

## Installed3.0.0 semantics

Reviewed exact installed source, whose implementation/binary hashes match the
original frozen native manifest. The producer's RadiationProblem(omega=np.inf,
water_depth=np.inf,rho=1025,g=9.81) is the documented native infinite-frequency API.
The [official3.0.0 intermediate cookbook](https://capytaine.org/v3.0.0/examples/B_intermediate_cookbook.html)
uses that construction explicitly. The
[problem setup manual](https://capytaine.org/v3.0.0/user_manual/problem_setup.html)
distinguishes angular omega from freq, period and wavenumber and gives the free
surface default0. Infinite depth and infinite angular frequency do not set
free_surface=inf; that separate parameter removes the free surface entirely.

Installed `bem/problems_and_results.py` `_get_frequencies` lines130-145 maps
omega=inf to symbolic infinity for omega/freq/k, with zero period/wavelength.
Finite deep water uses k=omega²/g. Forward speed is zero in our problems. The
RadiationProblem boundary condition contains the symbolic omega, and the native
solver's symbolic-aware linear solve preserves its finite coefficient. Pressure
adds another symbolic omega. RadiationResult.added_mass divides force by omega²
through matching symbolic factors, not literal floating-point inf/inf; damping
is explicitly zero at infinity. See `tools/symbolic_multiplication.py` and
`bem/solver.py` lines165-207. The exporter reads native named added_mass values in
the same influenced-row/radiating-column convention used for finite frequencies.

`green_functions/delhommeau.py` defaults gf_singularities='low_freq'. At k=inf it
automatically switches to 'high_freq' before calling the native binary. Finite
problems retain the default decomposition. Both are native representations;
this numerical representation difference plus coarse spatial discretization is
a possible explanation to investigate, not proof of a bug or justification to
fit Ainf. The installed Fortran source is not distributed, only pinned binaries;
this audit does not independently verify the compiled limiting kernel. No finite
large-omega replacement or change of GF method was made.

Cheapest already-scoped next mechanism control is the SAME68 physical hull with
native internal lid, at1/2/12/inf, all six modes. It probes irregular-frequency
contributions without claiming to isolate every source of finite/infinite mismatch.

## Native lid boundary

The [official3.0.0 mesh manual](https://capytaine.org/v3.0.0/user_manual/mesh.html)
documents `mesh.generate_lid()` for irregular-frequency removal. Partial interior
coverage is permitted and affects efficacy; panels on the free surface remain
experimental. This control fixes z=0, resolution0.25m and maximum512 lid faces,
reports actual area/coverage, and never treats numerical cap area as physical hull.

Installed `meshes/meshes.py:655` generates downward horizontal interior panels;
`bodies/bodies.py:85` preserves the hull separately, orients the lid down and builds
a hull mask. The six physical DOFs remain defined only on hull faces. Radiation
boundary velocity on lid entries is zero. Solver pressure is explicitly indexed
by hull_mask before force integration (`bem/solver.py:206`), so lid pressure does
not enter physical force moments or mass/hydrostatic metadata.

At infinite k the native GF still selects its high-frequency branch, including
when a lid is attached; the caller must not override that or silently substitute
a different free-surface problem. A native failure/degeneracy at this configuration
would be preserved and reported, not repaired or automatically retried.

## Exact source and resources

`reference/capytaine/export_lid_control.py`, controls
`reference/capytaine/lid-control-controls.json`; archived draft
`artifacts/reference/lid-control-20261024/proposal-draft.json`.

```powershell
.venv/Scripts/python.exe -u reference/capytaine/export_lid_control.py --controls reference/capytaine/lid-control-controls.json --output artifacts/reference/lid-control-20261024/attempt1
```

Exactly24 radiation solves, no excitation/new bank/new frequency sweep/272 refinement.
One CPU, maximum300s checked between solves; estimated20-90s, not measured for the
generated lid. Estimated memory<1GiB, not measured. Same GF initialization/solver
as no-lid scout. Freeze actual source/controls/input/dependencies/Python identities;
independent source PASS and coordinator allocation precede READY and atomic claim,
separate successful claimed readback and Run. Every nested command must exit0.
Exact terminal/children, idle under lease, separate idle proof, final release.

Expected invariants: exact68 source triangles/vertices/faces, DOF order and meanCG;
unchanged physical hull area/full M/K/CB/volume; lid horizontal/downward, interior
vertices, area≤exact mean waterplane, zero lid boundary velocity. Waterplane convex
hull area must equal exact GL2 K33/(rho*g) before it is used as the inside test;
this rejects an invalid convex fill for any different boat. Actual coverage is
reported. Arrays remain all36 raw A/B at1/2/12 plus direct Ainf; exact readback,
finite values and source/output hashes; retain all warnings and full difference
matrices versus pinned no-lid68 results.

If negative diagonal B(12) disappears while low-band changes are small, that
supports an irregular-frequency contribution. It does not qualify a radiation
kernel, high-frequency tail, mesh uncertainty or actual motion. If it persists,
do not assume lid-only repair or silently change lid/mesh/method. No positivity
clipping, symmetrization, empirical Ainf fit or production replacement is included.

## Executed24-call result

Exact source9a0f7dc8bb38894a61c5ab275b083c9f8def1fc21bc0c367ccf1e44643a2e8e5
reviewed PASS; READY status was the only controls change. Frozen source commit
3d5cedbb16a70df2de0974c8b8ec4e75a9ff4e7a, READY controls026c207a4dd18ba3e2811796b090e3934bd08deefccb427f7e4fffc4d74fc720,
proposalda8d99386adfae0f688875d2f27238b50f483027a813e9e534da70acea943ad5.
Output `artifacts/reference/lid-control-20261024/attempt1`; manifest
32dfcc286cef13171e1402e28efe8372a542174a6293ab4fbc40b32632f73c46,
postcheckd9927be87d5869135dd6666bfec599f519f322e0fa839d9459218bd77bdcc88f.

Claim750cf3 exit0, SEPARATE Verify96751f exit0, Run session99662; producer30420,
guard36072, UTC launch18:35:25.146810. Exact24 calls55.1033316s, one thread;
guard waited producer terminal0, tool terminalc6759b exit0, exact identities absent
after terminal1342db exit0. Idle64c5d5 exit0 under lease, separate idle readback
fb89e6 exit0, final removal973b54 exit0. No further numerical allocation.

Actual83 downward horizontal native lid panels cover9.754764m2,77.2667% of the
12.624790m2 interior waterplane. Physical68 hull faces and13.979549m2 hull area,
M/K/CG/CB/volume stay unchanged. All raw arrays/readbacks/identities finite and
unchanged source/output hashes pass. Native warnings retained.

| Quantity |No lid |Native partial lid |
| --- | ---: | ---: |
| Aheave(1),kg |19242.534 |19522.662 |
| Aheave(2),kg |12927.137 |13337.231 |
| Aheave(12),kg |-506.527 |-1454.250 |
| Bheave(12),kg/s |-48117.677 |48954.918 |
| Bpitch(12),kgm2/s |-90211.661 |217757.119 |
| Bsway(12),kg/s |-1035.022 |-88.178 |
| Byaw(12),kgm2/s |-170.621 |-4290.817 |

Direct Ainf changes only roundoff, max2.39e-12 over all entries (independent
reviewer). This native-lid sensitivity supports a numerical mechanism contribution,
but negative sway/yaw diagonal damping remains, so a passive operator was NOT
obtained. Partial coverage means persistence cannot disprove irregular-frequency
effects. Neither sign changes nor invariant Ainf uniquely attribute the entire
finite/infinite discrepancy or qualify a kernel. No corrected coefficients adopted.

Independent executed integrity review PASS: reviewer audited every frozen/output
hash, raw MAT metadata, lid geometry/coverage, terminal/idle proof and lease absence.
No independent solver rerun or physical correction. Next decisions remain with
coordinator; no refinement, full-band export or additional control is automatic.
