# ONE seeded-CIC/joint interface fixture: source ready, unrun

Later status: allocated attempt3 completed and independently reviewed PASS.
The source-only/readiness wording below records the earlier boundary; actual
receipt follows at the end. No initialized60s body run is authorized by this note.

2026-10-24. Parent authorized authoring/source review ONLY. No numerical launch
or actual initialized60s case is authorized by this note. Same failedcase/bank/
operator remain immutable. No maintained fork/library edit.

Implementation files: cic_seed_table.m, cic_seed_lookup.m, harmonic_seed_state.m,
make_harmonic_seed.m, set_initial_joint_state.m, apply_cic_history_adapter.m and
check_seeded_initialization.m. All house MATLAB source checks pass. Independent
source/math/API review is requested before source freeze/execution allocation.
Current fixture SHA256
`5385c4c7685c0cf65ba444a95a9bbe17fdd8c5fef63ef0f1f1476881d3c94f4d`.

The seed builder deterministically reads existing Q_kernel40/omega, computes
q/dq/ddq and exact worldomega/worldalpha from finite Euler kinematics, and freezes
q0/v0/history. It does not fit data or prescribe post0 motion. It labels this a
meanlinear approximation, not exact native stationary dynamics. Primitive initial
targets are explicit m/m/s/rad/rad/s; joint rotation0 and CG pivot map are used.
Force/motion actuation remains NoTorque/ComputedMotion, with no damping/restraint.

The CIC adapter adds the exact missing negative-time trapezoid contribution to
all three native force-history rows BEFORE the native force RateTransition and
PCHIP extrapolator. Actual function sample times are the lookup input. Important
source detail: paired RateTransitions initially expose time[-dt,-2dt,-3dt], so
force RateTransition InitialCondition must be fullpast forces at those same times.
The analytic history extends through -T-3dt (one sample more than the original
proposal). No private persistent-state mutation or postPCHIP addition.
Native actualpost0 velocities remain in the nativeCIC buffer. Original Ar delay,
its defaultzero output and extrapolation are UNCHANGED.

The ONE fixture job contains:

1. Actual copied native linear CIC subsystem, three lag points0:.02:.04,
   asymmetric full6x6 values, nonzero/distinct negative/current/post0 velocities.
   Only its kernel bus source is replaced by an explicit constant for the fixture.
   Log actual postRateTransition time/force windows, compare with independent
   scalar indexed fullhistory integration. Actual native PCHIP output AND both
   input-history logs must cover17fast timestamps0:.005:.08 while nativeCIC
   remains.02; compare every actual radiation output to the direct fullhistory
   PCHIP reference. Also query native extrapFrad at four fractional times per window.
   Check halfweight/currentv0 is counted once and correction is zero at N+2.
2. Kinematic derivative test: finite rotation derivative versus worldomega,
   worldalpha versus finite derivative, plus seed-file q0/v0 roundtrip.
3. Actual native Bushing Joint with copied native CG/initial transforms/Solid/
   world MotionSensor, nonzero three-angle pose/CG offset and all six velocity
   targets. StopTime0 INITIALIZATION ONLY; actual CG pose, vCG and worldomega
   must agree with independent targets. This is not a hydro body trajectory.

No new ODE solver/pole search/BEM/nonlinear adapter/physics correction is present.
Fixture evidence covers the interface, not accuracy/passivity/physical suitability
of the Jan90 sea or guarantee that initialized forced motion will remain bounded.

Prospective preparation command (copies/hashes only):

```powershell
.venv/Scripts/python.exe reference/wecsim/prepare_seed_fixture_manifest.py
```

After independent source PASS and coordinator numerical allocation, authoritative
engine readiness, atomic claim and SEPARATE successful claim verification:

```powershell
.venv/Scripts/python.exe reference/capytaine/geometry_execution_guard.py run --proposal artifacts/reference/generated-bank-consumer-20261024/seed-initialization-fixture-attempt3-source/proposal.json --owner /root/wecsim_reference
```

Exact child: run_matlab.py check_seeded_initialization.m, fresh
seed-initialization-fixture-attempt3.log. Output same-named directory under
generated-bank-consumer-20261024. One engine/oneCPU, expected10–60s, cap600s.
Preserve raw fixture models/logs before assertions. Failure stops before next
fixture phase, no automatic retry. Exact terminal proof, idle-under-lease,
separate readback and final release. Stop after independent fixture evidence;
no initialized actual60s case, renderer or extra control follows automatically.

## Attempt1 failure and bounded source recovery

Attempt1 exact client36832/guard4592 terminal1 before ANY simulation or interface
assertion: deleting the kernel bus block left a dangling line attached to native
Function1 input2, so add_line refused. Source98747cfa, log/source archive and
stopped partial cic_history_fixture.slx remain under attempt1 paths. Inspection
confirmed that line SrcPortHandle=-1. Auxiliary nonnumerical preservation commands
initially had shell quoting/eval-output/invalid-hyphen-model-name errors; the
model was ultimately saved under its native valid name, with no simulation retry.
All frozen identities matched, exact client/guard absent, idle-under-lease
separate readbackdb3644 passed and lease released FINAL.

Coordinator authorized ONLY six-line proper-API harness repair: identify the
known selectorOut1->Function1In2 line, assert both endpoints, delete that line
BEFORE deleting its blocks. New output/archive/log use attempt2. No helper,
kernel math, initial-state map, operator, bank or failed output changed. New
source passes m_check; independent tiny-diff review requested. No retry allocation
or simulation is implied by this source recovery. Before any allocated attempt2,
close ONLY the owned stopped partial cic_history_fixture after confirming its
saved attempt1 identity; do not close other models or delete failure files.

Attempt2 (sourcee257036f) fixed that exact line issue but compilation rejected
the fixture Velocity FromWorkspace: Interpolateoff cannot use its default
OutputAfterFinalValue=Extrapolation. Exact parameter enum readback confirmed
Holding final value is supported. NativeCICmodel/source/log were saved before
simulation; no interface/physics assertion ran. Client34992/guard38852 terminal1,
all frozen identities matched, separate idle readback67c9df passed and lease
released FINAL. Coordinator authorized source-only one-property harness repair
and freshattempt3. Source5385c4c7 changes ONLY outputpath2->3 and explicit hold
final value; input samples/interpolation/nativeCIC/math/helpers remain unchanged.
No numerical retry follows without independent tiny review and a new allocation.

## Executed attempt3 PASS

Final reviewed source5385c4c7/helpers unchanged. Exact output
`artifacts/reference/generated-bank-consumer-20261024/seed-initialization-fixture-attempt3`.
One engine9168/oneCPU/22.951s. All17fast timestamps0:.005:.08 aligned, CIC sample
time.02 preserved. Force history row error2.55e-15; actual compiled native PCHIP
error1.89e-14; additional native query error4.59e-14; correction afterN+2 exactly0.
Native mechanical StopTime0 actual pose/CGvelocity/worldomega agree to about
2e-16, with nonzero three-angle pose and pivot arm[.5,-.1,.05]m. E/Edot finite
rotation/derivative errors1.45e-11/2.27e-12. No Ar or operator/input/fork change.

Independent executed review verified all frozen hashes, all actual17fast grids,
direct scalar fullhistory rows (4.44e-16), SciPy PCHIP versus actual compiled
native output (7.44e-15), native joint/CG map and zero late correction. No
simulation rerun by reviewer. Interface evidence only: no exactnative stationary
solution, passivity, mesh/geometry adequacy or forced-case boundedness claim.

Client39780/guard1420 exact waitexit0(session13858/chunk4be76f), authoritative
both absent, synchronous fixture returned/models closed; all frozen identities
unchanged. Idle-under-lease recorded, separate readback18ec53 passed, FINAL
release. Attempts1/2 and their original source/log/partial models remain unchanged.
Coordinator conditionally authorized NEXT CASE SOURCE preparation only after
this checkpoint; actual initialized60s run still requires independent case-source
review and a new coordinator grant. No renderer or auto-body execution.
