# ONE initialized Jan90 native60s case: source ready, unrun

Later status: ONE allocated actual60s case completed with all runtime state/
force/history assertions passing; globalroll target provisionally met. Coordinator
accepted it for geometry-audit preparation ONLY. Independent actual physical
review/geometry/accuracy acceptance remain separate. Actual receipt below.

2026-10-24. Independent fixture EXECUTED PASS and four new case SOURCE PASS.
No actual body allocation is implied. Rejected held-release/calm-control/bank/
operator/fixture outputs stay immutable. No automatic refined mesh or encoder.

Reviewed sources: run_seeded_jan90_case.m
`8d8b8115417075206c4d7545c4b29d2810d5ad5a717a670e0f1df6be13b96ed1`,
prepare_seeded_jan90_inputs.m
`f5d0782315b2d349f5909c4d5a8179f328f979688b9c12706dcc019bd5180642`,
apply_seeded_case_initialization.m
`d33cdfa5068c5ea31d93204bc743a38ef038220e3d993f4c34c8b801d5c4d3c5`,
verify_seeded_cic_output.m
`b783e5f4da2edc1bcee54b92611f0ef2fcbe3a24c1921237326ea1d06bf37fd5`.
All seven independently fixture-validated helpers, old input config/builder and
proper case force hook are reused unchanged. House MATLAB checks pass.

Same explicit Jan05seed42/384whole-band MEM realization and headingoffset90,
same original17panel H5/M/K/Ainf/Ar/40sCIC,ode4dt.02, same unscaled exact native
component H forcing/no ramp. No coefficient repair or damping/restraint/filter.

Initial q0 and v0 come deterministically from the existing finite40sK particular
Q amplitudes, not a phase/trajectory fit. Worldomega/worldalpha use exact finite
Euler kinematics, with native Bushing initial targets and CG-pivot correction
verified by the fixture. Pre0 history extends -40.06..0, and the correct native
three force-history rows enter before RateTransition/PCHIP. Post0 velocities
are measured native dynamics, never a prescribed trajectory.

This is an approximate MEAN-LINEAR harmonic state/history, not an exact free
stationary native rotating-Solid/gyro/discrete-solver solution or nonlinear fluid
state. Native acceleration delay/default0/extrap remains UNCHANGED. Actual initial
acceleration minus harmonic seed and initial native wrench residual are reported
without compensation, cancellation or invented acceptance gate. Native actual
initial acceleration is solved by the original equations. Nonpassivity/mesh/
finitebody geometry limitations remain; state selection does not cure them.

All raw60bodycolumns and three CIC history/output logs are saved before assertions.
Known actual q0/vCG/worldomega, full H5/mass parameters and3001native output grid
must agree. All six signed applied excitation histories must equal exact H table.
Independent postprocessing convolution of specified pre0 world velocities plus
ACTUAL post0 world velocities verifies all three native force-history rows,
native PCHIP and the radiation force bus on every sample. Current endpoint is
counted once. No initial-sample rebasing, transient discard, drift/mode removal.
Failures preserve complete raw data and stop; no automatic second attempt.

Exact archive preparation (no numerical launch):

```powershell
.venv/Scripts/python.exe reference/wecsim/prepare_seeded_jan90_manifest.py
```

After root numerical grant, authoritative engine readiness, exclusive claim and
SEPARATE successful claim readback:

```powershell
.venv/Scripts/python.exe reference/capytaine/geometry_execution_guard.py run --proposal artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1-source/proposal.json --owner /root/wecsim_reference
```

Child: run_matlab.py run_seeded_jan90_case.m, fresh case-jan90-seeded-attempt1.log.
Output same-named directory under generated-bank-consumer-20261024. One engine/
oneCPU,60s/3001samples/dt.02/CIC40, expected30–180s, cap1800s. Exact process/child
terminal proof, idle-under-lease, separate idle readback, final release, then
independent actual state/force/history/dynamics review and geometry assessment.
No accepted-case/accuracy claim from successful execution alone.

Generated renderer is separately source-reviewed/owned and binds only a future
accepted trajectory after terminal. The seeded case must not wait on rendering,
nor does this proposal authorize an encoder or a second body run.

## Actual initialized case receipt

Output `artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1`.
Same independently reviewed source8d8b8115/f5d07823/d33cdfa5/b783e5f4 and validated
helpers/operator/bank. Total38.318s: preprocess26.316/compile4.381/dynamics3.306.
All3001samples0..60/.02 retained. Correct known q0 and nonzero CGvelocity/worldomega
readbacks passed; all3001x6 applied excitation errors EXACTZERO. Full independent
past+actualpost0 convolution, all3force-history rows and native PCHIP agree within
1.10e-11N/Nm; native radiation bus error EXACTZERO. Missing past fades as specified;
no measured future velocity is replaced by a harmonic trajectory. Ar default
startup/delay remained unchanged; no source/coefficients/input fitting.

| DOF |Maximum absolute native motion |Peak-to-peak native motion |
|---|---:|---:|
| Surge m |.36788|.67589|
| Sway m |1.54928|2.42297|
| Heave m |1.52026|2.54052|
| Roll rad |.08948494|.15683768|
| Pitch rad |.03007128|.05936681|
| Yaw rad |.03839487|.05735205|

Angular absolute peaks are approximately5.127deg/1.723deg/2.200deg. Final neutral
frame displacement[.123767,-1.549284,-1.039922]m and Euler
[-.08948494,-.01234669,.01458180]rad are retained, without recentering/detrending.
This is a provisional finite60s selected17panel model trajectory, not mesh/
passivity/geometry validation or proof of longer-time boundedness.

Actual q0=[.23520883,.13698402,-.25815063,-.00429300,-.00748473,-.00012470].
Initial CGvelocity[-.04924571,-.19377643,-.05531517]m/s and worldomega
[-.00233215,.00372168,.01119880]rad/s match the independent pivot-target mapping.
Initial radiation[-33.509,-48.205,-1045.390,46.235,-237.904,-33.215]N/Nm is seeded
and nonzero. Initial actual-minus-harmonic acceleration
[.254321,.124449,.174959,.046705,-.003252,.011289]m/s2/rad/s2 and native wrench
residual[-6404.65,-3134.03,-4406.05,-109.32,13.55,-124.45]N/Nm are explicitly
reported, not cancelled. Thus exact native stationary initialization is NOT
claimed; original numerical Ar startup and finite-pose/discrete differences remain.

Raw hashes:
- motion.mat d69652e3ca838afd2df6ddc249762d58321d833324038c51798399a36d055a95
- native-output.mat fe30d8eae1e082a44538cb8a40a5300098becd1cd63dbb4c63a07dbc998758b6
- native-cic-output.mat 2d889ef72e1a896fe6d10de9fbc0d826819b8fcf96ff4547bf0655c76c6d7d1e

Exact client7044/guard30896 terminal0(session17512/chunke6a8de), authoritative
both absent; engine9168child only nonnumerical mwdocsearch33628. All frozen
sources/inputs/dependencies matched through terminal. Idle recorded under lease,
separate idle readbackb436e0 passed, FINAL release. Independent actual stored
state/force/history/Newton-Euler review requested. Coordinator owns geometry
audit allocation and later artifact binding; no encoder/another body job follows.
