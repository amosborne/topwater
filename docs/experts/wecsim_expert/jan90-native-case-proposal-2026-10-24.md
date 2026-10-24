# ONE provisional native Jan05offset90 case: ready, not launched

Later status: ONE coordinator-allocated actual60s run completed. Runtime wiring
checks passed, but large motion REJECTS the intended linear candidate. See the
execution section below; the readiness wording preserves the prior proposal.

Updated2026-10-24. Final independent CASE SOURCE PASS at run script
`b58b5bb60ba1338c9f699361d033040f47564978869581591e3559dc36306ec4`.
Other reviewed source identities: prepare inputs `c3c40b9e...`, input config
`f3e2e32e...`, case force adapter `65650af0...`, model builder `b25e972d...`.
These are new case files; REF-001 and all stored harmonic data remain untouched.

Same immutable Jan05 whole384component realization/seed42/MEM spread and phases,
explicit heading offset90deg from incoming-aligned default. The initial mean
global bow yaw is the heading-file value4.41947980851rad, not a rendering trick.
Physics stays in unchanged mean-boat axes with the corresponding local wave
headings/phases. Global rendering maps position/orientation/field by the heading
file's R/t exactly once. Native roll is about the initial longitudinal direction.

## Physical initial-value assumptions

At generated t=0, etaCG=-.250572507061335m, sx=.032204329109661,
sy=.015621417478215. Deterministic pose is q0=[0,0,etaCG,roll0,pitch0,0],
roll0=atan(sy/sqrt(1+sx^2)), pitch0=-atan(sx); all velocities/worldomega zero.
Quaternion qy(pitch0)*qx(roll0) is converted to native axis-angle, checked against
Ry(pitch0)Rx(roll0). Initial pose is checked against this independent known pose,
never defined by subtracting the first simulated sample.

Body was idealized as held stationary at this pose before0, then holding forces
removed at0; radiation memory0 is consistent with that zero-velocity past. The
full stationary incident sea is already present. No ramp/prerun/transient discard.
Prior holding load is reported as the LINEAR estimate K*q0-Fexc0, not claimed as
an exact finite-wave load. The steady native mean-body diffraction field is an
explicit linear assumption, not an initialized nonlinear field at the actual pose.

Initial sampled footprint departure from the CG tangent plane is.183m and is not
small against.402m draft; plane-following pose does not guarantee exact curved
surface wet set/equilibrium. The actual trajectory is explicitly provisional
selected17panel model output, not trusted geometry/mesh/radiation validation.
Future nonlinear initialization would need a consistent incident/scattered field,
body/fluid initial state and holding/release history on actual geometry; no such
new model/adapter/trajectory comparison is requested or implemented here.

## Native implementation and assertions

`build_generated_model.m` repeats reviewed native World->official initialize->
native free6DOF/body->explicit reload order. Body.initial carries the pose; the
selected library uses separate meanCG translation and initial-displacement
axis-angle transform. Native ExternalForce/Torque and MotionSensor are World
resolved. Geometry remains mean-oriented STL minus CG once; no second heading
or equilibrium rotation is applied in physics.

Native official H5/M/K/Ainf/full inertia products and40sCIC are reused. Source
roundtrip Solidlocal+residualAinf is labelled PARAMETER identity. Preflight also
records R0*Isolid*R0.T, the world-resolved residual matrix and their actual
initial combination separately; it does not call fixedmean M+Ainf the effective
world tensor at nonzero R0. Preserve initial acceleration/net forces for later
native Newton-Euler review, including stock residual acceleration delay.

`prepare_jan90_case_inputs.m` prepares exact native H force for0:.02:60 and checks
independent scalar component sums. Scalar total eta history-20:.02:80 is only
official initialization compatibility, not the applied force. After native
preprocessing, `apply_jan90_force_table.m` disables the copied model-instance Boat
link via proper API, verifies one active positive imported branch and points its
FromWorkspace to jan90ForceTable. Maintained library remains unchanged; unused
compatibility excitation is archived separately. All signed3001x6 ACTUALLY APPLIED
force columns must equal the independently prepared full-bank table at native times.

Run assertions include exact H5 metadata, tensor parameter mapping, actual known
q0/v0/radiation0,3001native samples/0..60/.02, no PTO/mooring/damping/constraints,
and resolved movement relative to KNOWN q0. Raw native60columns and full trajectory
are saved before result assertions, with failure artifacts retained. Output q
stays relative to neutral meanCG; no recentering, mode removal or detrending.

## Exact command, resources and checkpoint

Prepare/archive source/input/dependency hashes (no engine/numerics):

```powershell
.venv/Scripts/python.exe reference/wecsim/prepare_jan90_case_manifest.py
```

After root grant, authoritative engine readiness and serial claim/separate verify:

```powershell
.venv/Scripts/python.exe reference/capytaine/geometry_execution_guard.py run --proposal artifacts/reference/generated-bank-consumer-20261024/case-jan90-attempt1-source/proposal.json --owner /root/wecsim_reference
```

Its exact child is run_matlab.py run_jan90_case.m with fresh case-jan90-attempt1.log.
Output is case-jan90-attempt1 under the consumer root. One MATLAB engine/oneCPU,
ode4dt.02/CICdt.02/memory40s,60s record; expected30–180s, declared1800s cap. Runtime
observation timeout is not terminal proof. On failure preserve whole run/report,
no automatic retry. On success exact terminal proof/idle-under-lease/separate
readback/final release, then coordinator review of actual output.

Renderer is a SEPARATE later artifact task, not a prerequisite to this run. It
will consume the entire native60s (all six histories), exact same generated sea,
native pose/world heading once, equal metres and no motion amplification. Start
with0/30/60 previews, fixed centered isometric view and native telemetry; choose
a truthful viewport from actual drift extent rather than silently removing drift.
No expensive encoder before native acceptance/preview review. No automatic second
mesh, wave/input change, extra BEM or accuracy claim.

## Actual execution: computational pass, physical candidate rejected

Exact reviewed/frozen source b58b5bb6 ran once in original engine9168. Output:
`artifacts/reference/generated-bank-consumer-20261024/case-jan90-attempt1`.
Full3001samples/60s/.02 trajectory and raw60columns were preserved. No retry,
prerun, initial-sample subtraction, detrending or encoder. Total43.762s including
preprocess29.624s/compile2.978s/dynamics4.352s. Source/inputs/dependencies match
their case-jan90-attempt1-source/proposal.json archive through terminal.

Runtime known q0=[0,0,-.2505725071,.01561205461,-.03219320279,0] matched;
initial six velocities and radiation force exactlyzero. H5/full parameter mass
checks passed. Prepared-versus-ACTUALLY-APPLIED excitation errors are EXACTZERO
on every axis/all3001samples. This rules out excitation replication or a wrong
force table at native output times as an explanation for the result.

Maximum absolute q: surge2.53339m/sway11.47412m/heave1.52325m,
roll1.876815rad(about107.53deg), pitch.155172rad(about8.89deg),
yaw.486592rad(about27.88deg). Final neutral-frame displacement is
[2.53339,9.76574,-1.04197]m; final Euler[-.0698141,-.0272406,-.2665042]rad.
All surge/sway/yaw drift is retained. The actual initial-value trajectory is
clearly different from the stationary particular response and is outside the
intended5–10deg linear-operating case. No trusted physical acceptance follows
from successful numerical/assertion completion. Instability/source attribution
remains for bounded independent stored-output review; no coefficient/input
repair or another simulation is authorized automatically.

Initial holding-load estimate [705.683,538.820,6899.327,544.361,-8782.670,-42.741]
N/Nm is retained with its approximate linear-model meaning. Native initial
total force is its negative to roundoff. Initial net force/acceleration, body
world inertia estimate, native residual added mass and complete force components
are available for independent Newton-Euler/feedback review at nonzero pose.

motion.mat SHA256190deabd054b255218f27ad602271ccaf2ce39274ea067f6a4a4009824457cd3;
native-output.mat SHA2560193c73974e68e47ea14a7cd05cd5610babe7267510c8aed879bd50312e51127.
Exact client16928/guard31416 waitexit0 (session9437/chunk458aa0), authoritative
both absent; only engine child was nonnumerical mwdocsearch33628. Idle written
under lease, separate idle readback0f04ba passed, lease released FINAL. No other
numerical job launched. REF-001 trusted baseline remains unchanged.
