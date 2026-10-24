# Jan05 heading90 internal screen: source-ready, 2026-10-24

Later status: the coordinator allocated ONE stored screen after independent
source/math PASS. It executed successfully; actual evidence is appended below.
The source-only language records the earlier proposal boundary, not a live hold.

Scope is SOURCE/PROPOSAL ONLY. No consumer numerical launch, canonical heading
adoption, new BEM/ODE/MATLAB/body work or input alteration. The goal remains ONE
credible generated case with5–10deg global roll and appreciable other five DOFs,
with its a-priori applicability and uncertainty stated honestly.

Prepared `reference/wecsim/compare_jan90.py`, SHA256
`0eff175461f59761f610c80ae74d68f4e278a64a4fbdb1e4f072883b46e90aeb`.
AST check passed; execution/independent review pending. It reuses the previously
reviewed normalization, native versus40sK harmonic Z solves, finite-band PV and
time-quadrature checks. It consumes ONLY the Jan05offset90 pair table from
`artifacts/reference/radiation-heading-scout-20261024/attempt1`; all bank entries,
index, phase, amplitude, frequency, heading and global/heading-file hashes must
match. New scout radiation17/68 data are NOT substituted into the existing kernel.

Exact proposed numerical command, after allocation and source review:

```powershell
.venv/Scripts/python.exe reference/wecsim/compare_jan90.py --scout artifacts/reference/radiation-heading-scout-20261024/attempt1 --output artifacts/reference/generated-bank-consumer-20261024/jan90-attempt1
```

One CPU, expected2–10s, maximum120s with source deadline. Output is one summary
and stored all-six native/CIC particular histories, differences,
complex response components and raw diagnostics. These are not actual trajectories.
The prospective manifest command is
`.venv/Scripts/python.exe reference/wecsim/prepare_consumer_manifest.py jan90`.
It archives exact scout/native/BEMIO/incidentJan90/bank inputs and numerical
dependency files. Do not claim/run until coordinator allocates after review;
then record-claim, separate verify, guarded run, exact terminal proof,
idle-under-lease, separate readback and final release remain mandatory.

New geometry proxy uses exact `Rrelative=Rplane.T*Rbody`, where both native body
and the joint local-plane orientation use extrinsicxyz=RzRyRx. Plane pitch is
`-atan(eta_x)` and roll `atan(eta_y/sqrt(1+eta_x^2))`. Identity-relative fixtures
check a plane-following body. All three relative angles are reported; the old
independent Euler subtraction remains separately labelled as a proxy. Actual
relative rotation is exact with respect to this FIXED-MEAN-CG tangent plane,
not to a curved surface at an advected/rotated hull. Full57point footprint
departures, initial eta/slopes/Hessian/eta_t and graph Hessian components are
retained as sampled extrema, not continuous/principal-curvature bounds. The
original incident diagnostic's separate Taylor bound remains separate. No
wet-set/flooding/physicalIC/linearity qualification follows.

## Practical bounded next step, rather than an open-ended band study

Use this ONE Jan90 diagnostic to decide whether its unchanged sea/heading even
supports the requested global-roll/relative-pose hypothesis. If it fails, report
that rejection and ask the coordinator to make the concrete target-versus-
applicability tradeoff; do not invent a sea or run a heading suite.

If a candidate is retained, the next dynamic step should be ONE actual native
WEC-Sim run with declared physical initial state/duration, all transient/drift
kept and all-six applied-force/mass/frame/H5 checks reused. This establishes a
provisional selected-model trajectory, not trusted mesh accuracy. It should not
be delayed by an unbounded frequency-tail campaign. Geometry applicability must
be explicitly limited rather than asserted from harmonic screening.

The minimum meaningful mesh-motion variation then consists of ONE matched
second native trajectory on one refined, independently checked producer mesh,
with the SAME incident bank, initial state, duration, timestep, radiation band/
memory and physical metadata. Only the six complete histories establish motion
variation; sparse scout matrices cannot replace the full paired excitation and
radiation operator for that second run. No automatic three-mesh convergence
campaign or per-coefficient acceptance threshold is proposed.

Cap's executed scout reports strong incident-band pitch A sensitivity and
negative high-frequency diagonal damping on68panels, plus Ainf changes. Thus
68panel coefficients are not automatically a qualified refined reference. Before
relying on its trajectory, producer/reviewer must decide whether that full fixed
operator is numerically usable; preserve raw pathology and stop/report invalid
data. Do not repair/symmetrize/fit/filter coefficients or call a difference from
an invalid solve an accuracy bound. A paired native run is conditional on that
bounded usability decision, not an authorization for more solver investigation.

Stopping point: source, exact command/resources and this next-step proposal.
Coordinator owns numerical allocation, canonical selection and any user tradeoff.

## Executed ONE stored screen

Output `artifacts/reference/generated-bank-consumer-20261024/jan90-attempt1`;
source/input/dependency archive `jan90-attempt1-source/proposal.json` under the
same parent. Source0eff1754 unchanged. OneCPU/1.828s/exit0. All input/indices/
phase/hash/normalization/rotation-fixture/ZQ checks pass; maximum ZQ residual
7.72e-12. No scout radiation coefficient replaced the original17panel operator.

| DOF | Native harmonic peak |40sK harmonic peak |Maximum history difference |
|---|---:|---:|---:|
| Surge m |.33890|.35495|.02591|
| Sway m |1.55635|1.60865|.05230|
| Heave m |1.51218|1.52103|.03020|
| Roll deg |5.11879|5.22092|.11082|
| Pitch deg |1.72912|1.71021|.10644|
| Yaw deg |2.31102|2.06576|.29871|

Exact rotation relative to fixedmean-CG local plane has native/CIC peak
roll2.10540/2.13313deg, pitch2.48406/2.50798deg and yaw2.30954/2.06414deg;
relative-heave proxy peaks.01804/.02977m. These are stationary particular
responses; Euler angles come from Rplane.T*Rbody. Sampled fullfootprint tangent
plane departure peaks.30671m, initial.18304m; graph Hessian component peaks
[.04039,.02230,.03750]permetre. These are sampled geometric diagnostics, not
continuous wet-set/flooding bounds or actual-body initial conditions.

The global5deg target hypothesis is RETAINED as a provisional candidate for
coordinator selection. It is NOT a qualified linear canonical case: relative
roll still exceeds the exploratory+/-1.5deg flatwater shape label; curvature,
actual initial state, radiation history, free-mode drift and motion uncertainty
remain unresolved. Do not tune geometry/input/criteria to turn this into a pass.

Exact producer39260/guard25492 terminal0 (chunk967da7); authoritative CIM both
absent; all frozen identities match. Idle written under lease, separate
idle-under-lease readback77fb0c passed, then lease released as final operation.
Independent executed review requested. Stopping point reached; no actual WEC
body run, new BEM/ODE/MATLAB work or canonical adoption is authorized by this run.
