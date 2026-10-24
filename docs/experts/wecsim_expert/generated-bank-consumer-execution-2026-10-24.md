# REF-002 official import and harmonic ranking, 2026-10-24

Status: allocated consumer preparation completed; independent executed review
requested. No REF-002 body trajectory, canonical sea selection or dynamic
linearity qualification. No custom time-domain integrator was implemented.

## Inputs and official conversion

Native source `artifacts/reference/generated-bank-native-20261024/attempt1`
completed2838native solves on the unchanged17panel mesh. Three DEFAULT incoming
aligned bank tables preserve384components/record, seed42, full observed spectrum,
MEM directions and phases. Mesh warnings/uncertainty remain unqualified.
Maintained fork inspected clean at `f962103a348d63101ce669c2611f7d7c109755c3`.

`reference/wecsim/prepare_generated_bemio.m` executes official readCAPYTAINE,
radiationIRF40s/2001time points/2401integration-frequency points, compatibility
excitationIRF+/-20s, and writeBEMIOH5. It restores exact native M/CG/CB/volume/K
and direct Ainf, retaining rounded legacy metadata separately. The compatibility
heading0 excitation is not used as the generated-bank force.

Output `artifacts/reference/generated-bank-consumer-20261024/bemio-attempt1`:
canonical.h5, bemio.mat, interface-check.json and h5-readback.json. Single MATLAB
engine9168; conversion11.616s; phase real/imaginary residual1.455e-11. Subsequent
MATLAB h5read verifies exact CG/CB/volume, normalized full K/Ainf, and all36
radiation kernel entries against saved BEMIO data. No library/fork edits.

## Harmonic-only screen

`screen_generated_banks.py` builds dimensional total force from exact native
complex H using real(sum a H exp(i*phase-i*omega*t)). It checks three independent
analytic signs with distinct six-axis scales, actual exported harmonic quarter
cycles, and independently summed scalar versus vectorized totals. No extra
imaginary conjugation or CG phase. Synthetic error<=3.68e-16; actual quarter
error<=4.55e-13; whole-bank scalar/vector force error<=3.21e-10N/Nm.

Stationary harmonic response solves full coupled
`[K-w^2(M+A)-i*w*B]Q=a*H*exp(i*phase)`. A/B are linearly interpolated in omega
between the new native .05:.05:12rad/s grid, never treated as direct pair solves.
Dimensional scaling preserves the full generally nonsymmetric matrix; maximum
scaled condition2284.05. No pseudoinverse, damping addition or coefficient repair.

Finite-window harmonic absolute peaks, t0:.02:60 (not native trajectories):

| Record | Surge m | Sway m | Heave m | Roll deg | Pitch deg | Yaw deg |
|---|---:|---:|---:|---:|---:|---:|
| Apr24 2022 |1.31144|.72534|1.58150|7.55092|3.913|2.500|
| Mar07 2025 |1.03195|.89591|1.83010|10.20068|5.078|2.408|
| Jan05 2025 |1.52526|.35740|1.50863|3.80976|2.790|2.830|

Full peak-to-peak/RMS, complex component responses, prepared force tables and
stationary histories are stored in
`artifacts/reference/generated-bank-consumer-20261024/harmonic-attempt1`.
Screen elapsed2.359s. Its summary ranks Mar, Apr, Jan by harmonic roll; this
ranking alone does not select the case or establish admissible relative heel.

Finite40s kernel diagnostics reconstruct B=int K cos and
A=Ainf-int K sin/w at actual incident frequencies. All36 absolute A/B/impedance
errors and low/high-frequency B endpoints are reported in summary.json; these
are diagnostics, not a passed accuracy threshold. Last5s absolute kernel
integral fraction ranges.003663 to.012074 across entries. In particular, the
heave-heave impedance discrepancy reaches about34555 in dimensional impedance
units for one bank; no claim that the finite-kernel model equals the direct
frequency-domain model is made. Tail size alone does not bound response error.

Zero-restoring surge/sway/yaw homogeneous startup and drift remain unresolved.
The harmonic response is not a held-release initial-value solution. Initial
immersion, actual curved-wave geometry, relative attitude and geometric linear
operating envelope must be considered before choosing ONE native WEC-Sim run.
Global7.55deg harmonic roll is not proof of small relative heel or safe geometry.

## Execution and recovery evidence

Separate manifests/source snapshots:
`bemio-attempt1-source/proposal.json` and
`harmonic-attempt1-source/proposal.json` under the consumer output root.
Each used exclusive lease, record-claim, successful separate verification,
guarded launch, exact child wait/authoritative terminal proof, idle-under-lease,
successful separate idle readback, and final lease release. BEMIO session43550,
client32972/guard27360 terminal0. Harmonic guard39296/client32308, sequential
H5client27612 and screen35320 all terminal0 (chunkac9db3), all four authoritative
CIM identities absent. Final idle readback10496a passed; lease released. No
further numerical execution is allocated by this completed sequence.

Preparation deviations preserved:

- Restricted Get-Process hid original MATLAB9168. An unnecessary owned hidden
  duplicate bootstrap25648->40536->30448 was created, then exactly identified
  and terminated BEFORE numerical work. Authoritative engine query proved
  original9168/stopped libraries, which was reused. Do not infer replacement
  from restricted process visibility; authoritative query must come first.
  Records: engine-bootstrap.json and duplicate-bootstrap-reconciliation.json.
- Optional h5py was absent during harmonic manifest preparation. No lease or
  numerical spawn occurred. Incomplete archive was preserved as
  harmonic-attempt1-source-incomplete-manifest. Actual H5 readback moved to
  built-in MATLAB h5read; no dependency installed. The sequential wrapper waits
  each exact child and fails closed before the next stage on a nonzero exit.

Next: independent executed consumer review, then coordinator's ONE-case initial
state/relative-geometry decision. No body launch, wave changes, nonlinear adapter
or Topwater solver work follows automatically.

## Subsequent allocated stored-data comparison

`compare_generated_radiation.py` source
`e65ba863ba6b13d0e47a3e9287f645ea88306226b49c18acc0a2c49f9b8f2c95`
passed independent source/math review before freeze. Its ONE allocated run
completed1.859s/oneCPU/exit0, no BEM/ODE/body integration. Output:
`artifacts/reference/generated-bank-consumer-20261024/compare-attempt1`.
The source/inputs/dependency archive is `compare-attempt1-source/proposal.json`.

Executed checks: dimensional native A/B/Ainf/K versus official normalized BEMIO;
actual kernel0:.02:40 and integration-frequency endpoints; constant finite-band
PV fixture; exact removable finite-time cosine/sine weights; three joint RyRx
plane-normal closure fixtures; full coupled native and finite-kernel ZQ=F
residuals<=5.76e-12. No input/Ainf/damping fitting or modification.

The analysis reproduces Bstar(nu)=rho*nu*linear_interp(normalized B), keeping its
difference from linear interpolation of dimensional B explicit. It compares
the stored40s time-grid kernel transform with exact-in-time40s cosine/sine
weights on the same2401frequency trapezoid, and an infinite-time finite-band
PV Kramers-Kronig transform using pole subtraction plus the analytic constant
term. Raw transformed matrices and effective Ainf are retained. Native A12 is
not assumed asymptotic: A12(heave,heave)-directAinf=-8230.20071kg. Band cutoff,
coarse mesh/native consistency and frequency/time quadrature effects cannot
all be distinguished from this evidence; no coefficient acceptance gate adopted.
Optional actual virtual-mass-delay transfer was excluded pending its own source
derivation; neither harmonic model reproduces the native feedback time stepping.

All six particular histories use identical H, generated a/omega/theta/phi, dryM
and K. Units below are metres for surge/sway/heave, degrees for roll/pitch/yaw.
For April, native peaks[1.31144,.72534,1.58150,7.55092,3.91297,2.50019] become
finite40sK peaks[1.29706,.74118,1.59584,7.80214,3.99503,2.26153]. Maximum absolute
history differences are[.02104,.01820,.05393,.39567,.26906,.32952]; RMS difference
over native RMS is[.01903,.02906,.03823,.05273,.07691,.11831]. March roll is
10.20068 versus10.22845deg, but maximum roll-history difference.58060deg and
heave-history difference.07441m. Complete three-record/all-six peak/peak-to-peak/
RMS/difference tables and component complex responses are retained in summary/MAT.

The fixed-mean-CG local-plane comparison uses joint plane pitch=-atan(eta_x),
roll=atan(eta_y/sqrt(1+eta_x^2)), verified against normal[-eta_x,-eta_y,1]. Euler
angle subtraction is explicitly a SCREENING PROXY, not exact finite relative
rotation, moving/rotated-hull geometry or wet-set assessment. Native/finiteK
relative-roll proxy peaks: April5.92894/6.06868deg, March5.89417/6.33391deg,
January3.51997/3.55535deg. April relative-heave proxy peaks.03487/.04843m;
its full-footprint plane departure reaches.37225m (initial.06884m). These values
do NOT support claiming global5–10deg roll with demonstrated small relative
heel. The static righting departure at5deg remains material, and no actual
initial pose, trajectory, startup momentum/drift or linear operating envelope
has been qualified. No criteria relaxed to reach the requested roll.

Exact producer32716/guard22168 terminal0 (chunkbd61ae), authoritative CIM both
absent, all frozen source/input/dependency identities unchanged. Idle written
under lease, separate idle readback99bbad passed, then lease released as final
operation. Independent executed review PASS received: all frozen/output hashes
match and reviewer independently recomputed three-MAT history differences.
Reviewer decomposed shared-frequency heave A discrepancies: sampledK versus
exact finiteT maximum2.95kg, exact finiteT versus truncated-band PV9.00kg,
truncated-band PV versus native3573kg. Thus40s/time quadrature alone does not
explain the dominant discrepancy; missing frequency tails versus native coarse
mesh/irregular-frequency inconsistency remain unresolved. Stopping point reached;
next selection/initial-state/tradeoff decision belongs to the coordinator.

Coordinator subsequently requested ONLY signed reductions of the completed
stored arrays, explicitly without a new solver/job. `report_signed_band.py`
preserves exact source-array hashes in `signed-band-report.json` under the
consumer root. All three banks share actual omega values: A33native-A33band
ranges from-3573.03645kg at3.10898474rad/s(.49481029Hz,index384) to+375.78539kg
at.190015218rad/s(.030241861Hz,index1),314negative/70positive entries. Raw native
B33 is positive on all240computed frequencies (min1.056109 at.05rad/s,
max28536.0999 at4.25rad/s); incident-interpolated B33 is positive at all384pairs.
For passive omitted upper-frequency damping, its KK contribution would be
positive throughout this incident band, so an upper-tail-only explanation cannot
explain the dominant negative difference. Positivity at sampled frequencies does
not prove uncomputed-tail passivity; lower-tail/native numerical causes remain
unresolved. No physics correction, new BEM/ODE/body run or ledger/lease mutation.
