# Current state: WEC-Sim / BEMIO expert

Updated2026-10-24. Authorization/phase: [handoff](../../checkpoint-handoff.md).
Coordinator owns selection/acceptance/allocation. No automatic next numerical job.

REF-001 trusted minimum reference is COMPLETE/frozen:10s/501native samples,
actualrest and full6force/mass/frame checks, final centered equal-scale video
independently audited. [Record](minimum-reference-2026-10-24.md).
Maintained user fork pin f962103a348d63101ce669c2611f7d7c109755c3 (officialv7.1.0
base): reviewed import history/component sum/Solid inertia order/selected rigid
force dispatch corrections. [Pin/scope](../wecsim-fork.json). No upstream/main push.

REF-002 goal remains ONE fresh buoy-generated case, roll5-10deg plus appreciable
other5DOF, a-priori geometry applicability stated. No spectrum/phase/amplitude
fit, suite, nonlinear adapter or custom integrator. Topwater solver remains held.
Three immutable seed42/full MEM384component banks are internal candidates only.

Official new BEMIO40sCIC/H5/fullM/K/Ainf readbacks and harmonic-only comparisons
passed integrity review. Native/CIC particular responses differ materially, but
DO NOT qualify startup/free-mode drift. [Stored evidence/limits](generated-bank-consumer-execution-2026-10-24.md).
Jan05offset90 harmonic roll5.119/5.221deg and fixedmean-plane relative roll2.105/
2.133deg retained one provisional candidate; sampled initial/peak plane departures
.183/.307m do not establish actual wet set. [Screen](jan90-screen-proposal-2026-10-24.md).

ONE actual60s Jan90held-release case completed43.762s: all3001signed6appliedforce
errors EXACTZERO; independent knownq0/v0/memory0/H5/parameterM and native physical
balance checks pass. BUT roll~107.53deg/sway11.47m/yaw27.88deg REJECT intended
linear candidate. Entire raw60s preserved, no crop/rerun/encoder. Neutral-world
free-mode drift/heading departure precedes late runaway; mechanism not proved.
[Actual case/source/history](jan90-native-case-proposal-2026-10-24.md).

Stored full radiation matrices have guaranteed nonpassive harmonic directions,
even in incident band (symmetric scaled native B min-184.99,40sK min-196.28).
Sampled diagonal B33 positivity was insufficient. Symmetric M+Ainf is positive;
actual residual gain/frozen feedback checks do not establish a simple mass-loop
instability. Cutoff/coarsemesh/native consistency remain unresolved; no coefficient
repair/damping fit accepted. Report: case-jan90-attempt1/stored-passivity-diagnosis.json.

ONE calm-water control, same unchanged17operator/CIC40/dt.02, initial0.5degroll,
zero allv/memory/excitation, completed30.555s. Independent executed review PASS:
roll decays to.00004258deg, all3001zero force/wrench/Ar checks pass, sway drift
.01576m retained. This ONE small perturbation does not prove passivity/all-state
stability or forced-case validity. [Executed receipt](free-decay-control-proposal-2026-10-24.md).

CURRENT assignment: SOURCE/PHYSICAL IC PROPOSAL ONLY for future same Jan90sea
with deterministic particular q/v, native worldomega and consistent CIC plus
acceleration-delay prehistory. No implementation/numerical grant. Cannot merely
set q/v while zeroing memory, prescribe a trajectory, suppress an unstable mode,
or use initialized state as proof of passive/accurate hydrodynamics.

All owned jobs terminal/leasefree; shared execution ledger owns current liveness.
Before engine replacement query authoritative outside-sandbox IPC/exactPID:
restricted process visibility is insufficient. Historical leads/setup live in
linked records; no untested marine-gyro workaround is accepted.
