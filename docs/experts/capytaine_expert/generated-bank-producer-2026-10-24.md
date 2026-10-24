# Generated-bank native producer proposal

Status: SOURCE/INTERFACE REVIEW PENDING; no numerical launch authorized by this
document. Scope is one bounded producer batch for all three unchanged generated
records, DEFAULT incoming-aligned offset0 only. Coordinator owns allocation and
selection of ONE later WEC-Sim case. Producer owns no response integrator.

## Immutable geometry and sea identity

Source: `reference/capytaine/export_generated_banks.py`; controls/interface:
`reference/capytaine/generated-bank-contract.json`. Same canonical physical input
ea18844904e03bbfd72abc3396240c97b1994b3584e36c8149f8ae5f60a99e4f,17 unchanged
triangles from accepted `native-gl2/native.npz`, native centroid hydrodynamic
quadrature, mean native modes/CG. Exact accepted GL2 K/M/CB/volume and hull STL
are reused as metadata; new radiation/diffraction is computed directly.
All Capytaine implementation/binary hashes must match the accepted source.
Capytaine3.0.0 API source was inspected for RadiationProblem, DiffractionResult
forces and native `froude_krylov_force`. There is no new geometry or force correction.

Banks are the three records in
`artifacts/reference/generated-sea-candidates-20261024/manifest.json`. Only their
`heading-offset-0deg.json` files are consumed. All384 components per record are
retained, INCLUDING zero-amplitude entries; no frequency/heading Cartesian product,
direction binning, interpolation of H, phase fitting, spectrum scaling or cropping.
Each exact omega/theta_body pair gets a direct native diffraction solve and native
FK integration. Bank component order, amplitudes and local phases are copied.

Native frame stays the existing mean boat frame. Global bank is unchanged;
theta_native=theta_global-yaw modulo2pi. Global position is Rz(yaw)*x_native+t,
with R/t copied from reviewed headingfile, and phase_local already includes k.t.
The producer must not apply another yaw or phase shift. Native x is roll axis;
vectors/tensors and eventual rendered world poses transform exactly once.
This is fixed initial mean geometry/modes, not adaptive heading hydrodynamics.

## Outputs and computation count

| Product | Direct new native calls | Meaning |
| --- | ---: | --- |
| Radiation at omega=.05:.05:12rad/s | 1440 |240 frequencies x6 radiating modes, all36 couplings |
| Direct infinite-frequency mass |6 |all6 radiating modes |
| Heading0 diffraction on radiation grid |240 |REAL compatibility data for official BEMIO mandatory variables |
| Three exact384-component pairs |1152 |one diffraction solve per original pair, FK from native API |
| Total |2838 |no old BEM coefficients reused, no hidden refinement |

`canonical.nc`, `Hydrostatics.dat`, `KH.dat`, `native.mat`, `boat.stl` support
official importer. Heading0 excitation is compatibility data; generated forces
MUST use paired tables. Do not fabricate missing Cartesian pair cells in NetCDF.
`body_name` compatibility alias preserves native `body`. Native MAT retains full
dimensional M/K/cg/cb/rho/g/volume with fresh direct Ainf. `radiation.mat` stores
omega(N,1), A/B(N,6,6), Ainf/M/K and ordered DOF names.

Each `<observation_id>-pairs.mat` contains component_index1:384, omega_rad_s,
heading_native_rad, amplitude_m, phase_local_rad (allN,1), fk_complex,
diffraction_complex, excitation_complex (N,6), ordered Surge/Sway/Heave/Roll/Pitch/Yaw,
CG/rho/g, globalbank/headingfile SHA256, yaw_global_rad, global_from_native_R/t.
H=FK+diffraction, dimensional per unit elevation (N/m or Nm/m), exp(-iomega*t).
No BEMIO imaginary-sign conversion is applied to paired complex H. Consumer
reconstruction is Re[sum a_j H_j exp(i phase_local_j-iomega_j*t)] at the native
spatial origin. Native H already includes mean hull spatial dependence relative
to that origin; CG must not be inserted as an extra phase correction.

Radiation grid is separate from physical incident band. Consumer may use declared
componentwise linear omega interpolation of full A/B for a preliminary harmonic
RAO ranking ONLY. Those are interpolated values, not direct native pair radiation
or a radiation-kernel convergence check. No interpolation is performed in this
producer. A/B are not symmetrized and no offdiagonal terms are dropped. Any later
direct pair radiation solves require a scoped amendment/allocation. The requested
40s CIC memory/dt.02 and60s record are consumer proposals, not certified by this
export. Recurrence2pi/.05=125.66s alone does not establish tail/cutoff accuracy.

## Command, limits and evidence

```powershell
.venv/Scripts/python.exe -u reference/capytaine/export_generated_banks.py --contract reference/capytaine/generated-bank-contract.json --output artifacts/reference/generated-bank-native-20261024/attempt1
```

One CPU thread; proposed elapsed budget600s, estimated180-300s from baseline
1086calls/69.98s; estimated peak memory<1GiB for17panels, no measured bound claim.
Time budget is checked between solves; a slow individual solve must be monitored
to actual terminal, never inferred killed at budget. Unique output directory,
partial manifest/checkpoints/source snapshots retained on failure, no auto-rerun.
No new BEM launches until joint source/interface approval and sole owner allocation.

Before claim, archive exact command, source/contract/canonical/bank/accepted-native
and all numerical dependency identities. Use the coordination CreateNew lease,
successful claimed ledger write, SEPARATE successful matching readback before
launch. Inspect every actual nested command exit code. Record exact process/time,
wait children terminal, write idle under lease, SEPARATE idle-under-lease proof,
release final. The geometry attempt1 protocol deviation remains separately recorded.

Acceptance:2838 direct calls terminal0; all finite raw coefficients/full metadata;
all bank arrays unchanged and MAT-readback exact; H=FK+diff; all17triangle identities;
official NetCDF variable/DOF shapes and exact GL2 metadata; Ainf finite/full6x6;
source/dependency/input identities unchanged; output/source hashes and terminal
evidence. Consumer performs official BEMIO roundtrip and independent complex-phase
reconstruction, researcher reviews frame/convention/counts. No coefficient tolerance,
mesh accuracy, temporal accuracy, target-roll success or nonlinear validity is
claimed by successful export. Geometry curves already show relative-flatwater5/10deg
is a substantial fixed-K departure; actual relative-wave geometry must be checked
before accepting any global-roll interpretation.

Declared minimal post-export checks (no additional BEM): analytic H=1 and H=i
quarter-cycle reconstruction signs; one nonzero actual harmonic FK compared with
independent centroid incident-pressure integration rho*g*exp(k*z+i*k.xy), including
all six signed normal/moment integrals about meanCG; multidirectional nonzero-phase
scalar sum against matrix reconstruction; global/local eta closure at three points
and times using copied R/t and phase exactly once. Verify every serialized index,
bank amplitude/phase/frequency/direction, transform/DOF/sourcehash and finite output.
Native applied-force histories are a later consumer runtime acceptance item, not
part of this coefficient batch. Source math/frame review PASS at3c15c37e,
contractb992fe7f (researcher); joint WEC approval remains required before READY.

Stopping point: reviewed ready proposal/source/schema/resources; await coordinator
execution allocation. No dynamics, sea changes, extra heading screen or renderer.
