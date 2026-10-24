# Native harmonic energy witness: exactly3 solves, source only

Status: independent final source/math/API review PASS at the source/controls/archive
identities below; no READY
status, numerical fixtures or native solve executed. Latest physical requirements
are rest/zero velocity/zero radiation memory, viscous off, full60s startup and
physics-based pruning. Old unpruned stationary-history prototypes and stopped
geometry audits retain their own policies; none becomes the final case.

## Model and sign basis

Capytaine3.0.0 solves full linear potential-flow radiation/diffraction, with
incident FK PLUS diffraction excitation; it is not an FK-only long-wave model.
The [official theory](https://capytaine.org/v3.0.0/theory_manual/theory.html)
assumes inviscid/incompressible/irrotational flow, small waves and body motions,
linear free surface and flat bottom. The project's wavelength>=2beam support
rule is separate from numerical spatial/time/IRF coverage and is not a universal
diffraction-theory cutoff. Physical incident pruning does not repair radiation
integration outside the incident band.

Installed3.0.0 source checked:
`bodies/bodies.py:238` integrates fluid-on-body force as minus pressure times
outward normal and DOF motion; `bodies/dofs.py:90` uses axis cross(r-CG) for rotations;
`bem/problems_and_results.py:448` applies normal velocity=-iomega*dof.normal;
`bem/solver.py:191` gives p=iomega*rho*Phi, retaining physical hull pressure only;
`bem/problems_and_results.py:582` extracts A=Re(F)/omega² and B=Im(F)/omega.
Thus Ffluid=(omega²A+iomegaB)q, v=-iomega*q. Translation/rotation generalized forces
are work-conjugate N/Nm with metre/radian displacements about the same meanCG.

For real coupled c and dimensionless amplitude zeta, q=c*zeta. Mean body-on-fluid
power is minus.5Re(v*F)=.5omega²zeta²*c.T*B*c. For real c the raw reactive A term
contributes zero mean work even if its numerical matrix is nonsymmetric. General
complex motion instead requires Herm(B-iomegaA), not an unqualified B-only test.
No raw arrays are symmetrized or overwritten. Symmetric-part quadratic reductions
are diagnostics only. Instantaneous memory return or negative component work is
not the same as negative total steady cycle-averaged power.

## Exact stored witness

Read-only reductions use exact raw17-panel `radiation.mat` rows, no interpolation:
omega2.2(row43) and6.2(row123), indexed from zero. At2.2, every diagonal B entry is
positive but the fixed real coupled witness has c.T*B*c=-185.0183778668kgm2/s,
outgoing mean power=-447.7444744376W per zeta². At6.2 the previously stored real
witness gives-1756.8214690876 and-33766.1086358632 in the same units.
Source values: `generated-bank-native-20261024/attempt1/radiation.mat` and
`generated-bank-consumer-20261024/case-jan90-attempt1/stored-passivity-diagnosis.json`
under artifacts/reference. The latter supplies the FIXED vector, not a new fit.

Normalization D=diag(2.7432m,2.7432m,2.7432m,1rad,1rad,1rad), c=D*y. D is a
work-conjugate coordinate scaling, not a change to boat dimensions. Fixed c:
[2.2379192986776766,-.2710067725448315,-.14130201465424522,
-.1395146161773907,-.49395942403515825,-.2420474994948622] metres/radians per zeta.
Its induced physical surface displacement is c_translation+c_rotation cross(r-CG).
The same field and c are used on17 and68 triangles without re-optimizing.

Omega2.2rad/s is inside the project's physical wavelength envelope, but17 panels
are numerically unqualified. This is a discrete-operator concern, not proof of
a Capytaine implementation bug or a justification for an arbitrary band cutoff.

## Exactly3 native solves

1. Same17 physical centroid hull, composite prescribed displacement at2.2rad/s.
2. Same17 hull, pure heave at1rad/s as positive sign/postprocessing control.
3. Same68 four-way-subdivided physical surface, SAME composite at2.2rad/s.

Composite WorkWitness is a prescribed linear combination of the existing six
rigid modes for this frequency-domain diagnostic, not an additional physical
degree of freedom in a boat dynamics model. Full six physical force components
remain saved. No lid,272 panels, extra frequency, diffraction solve, time-domain
motion or automatic cascade. Positive control does not qualify an entire band.

All three use native indirect solver, same Delhommeau defaults/centroid integration,
zero forward speed, free_surface0, depthinf, same rho/g/CG/full physical metadata.
Per-parent oriented area/first/second moments, plane, winding and containment prove
unchanged17→68 surface; refined triangles also match the frozen prior68 artifact.
keep_details=True retains potential, pressure, source strength and BC. Each native
call is counted immediately after solve, before postprocessing can fail.

Independent normal-area/moment integration reconstructs physical forces and panel
work .5Re(sum p*conj(vnormal)*area), compared with minus.5Re(v*F). Native extra-mode
force is checked against c-projected six physical forces for composite cases. The
two17 cases additionally compare all six forces and total work with exact stored
A/B at the identical omega, without interpolation. These nearfield equalities
verify the interpretation of one discrete solution, not physical convergence.

Native `post_pro/kochin.py:21` requires retained indirect sources and divides its
centroid source integral by4pi. Compute H at360/720 periodic angles and independently
repeat that pinned source reduction. Save nonnegative integrals of abs(H)² and
angular-grid differences. This norm alone is tautologically nonnegative; its
nonzero outgoing field alongside negative total input is an energy inconsistency
only with the outgoing-wave sign/model convention justified. The two reductions
use the same BEM source solution and are not an independent numerical solver.

**Absolute Kochin power is explicitly NOT claimed.** The printed official far-field
asymptote equations have an unresolved amplitude-factor inconsistency, as the
independent reviewer noted; positive constants agree in sign but differ in power
normalization. Installed GF normalization must be independently derived/verified
before any magnitude closure claim. No constant is fitted to nearfield damping.

## Source, controls and resources

Source `reference/capytaine/native_energy_witness.py`, SHA
60c618b779caaa089d8b6ec880f94d5ccb2f0227cbafb3d3900bed2182773bfc.
Controls `reference/capytaine/native-energy-witness-controls.json`, SHA
c7a3605d1141cb389888d7d36af409141cc191a5f5ed2f54eb3fa313bb52bb9f.
Archive `artifacts/reference/native-energy-witness-20261024/proposal-source-only.json`,
SHA8689cc83cf6874c275e2b0dbbe88caff6c6056bc3a699ea5e685d1e7dc8c1c96.
Original exporter/scout helpers and native dependency implementations/binaries are
pinned. AST only passed; prepared algebra fixture (positive B work and zero real-A
reactive mean work) has not executed. No numerical setup demonstration.

```powershell
.venv/Scripts/python.exe -u reference/capytaine/native_energy_witness.py --controls reference/capytaine/native-energy-witness-controls.json --output artifacts/reference/native-energy-witness-20261024/attempt1
```

Prospective one CPU/max180s, expected55-75s dominated native GF initialization;
memory designed below1GiB, not a measured resource result. Three small systems,
bounded720angle arrays. Checks between calls; exact process must still be monitored
to terminal on budget/failure. Freeze sources/controls/physical data/native evidence/
dependencies before any later allocated claim, separate verified readback, run.
Exact terminal/children, idle under lease, separate idle proof, final release.

Outputs: per-case details.mat (raw pressure/potential/sources/BC/geometry/modes/H/
physical metadata), work.json, manifest/source snapshots/raw warnings and hashes.
Any negative power, unchanged sign or mesh sensitivity stays raw. Three points
cannot select a frequency cutoff, qualify a band/IRF, uniquely attribute a cause
or justify a coefficient correction. Stop at reviewed source readiness pending
coordinator grant; no native solve or followup authorized by this document.

Independent reviewer verified all12 frozen input hashes and every listed installed
dependency file hash, the actual reactive-A fixture, full six-force projection,
immediate solve counter and exactly three definitions/same fixed68 field. No
fixture or solve was run by reviewer or author for readiness. READY status and
fresh proposal/HEAD freeze still require the coordinator's later allocation.
