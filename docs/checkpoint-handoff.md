# Checkpoint handoff

Updated 2026-10-24. Read before choosing work; latest user instructions control
scope. Stable routing is in [AGENTS.md](../AGENTS.md).

## Authorization and phase

User authorized completion of the minimum-mesh end-to-end reference. REF-001,
CAP-001 and WEC-001 in [assignments](experts/assignments.md) are active: ten
seconds of actual-rest six-DOF motion, upstream incoming waves and isometric
video. Mesh refinement is subsequent work. Topwater implementation, tests,
numerical runs and optimization remain HELD.

## Current evidence and next priority

Minimum Capytaine export completed: 17 panels, 1086 solves, 69.9800752 seconds,
terminal exit0, session84329. All14 producer finite/identity/metadata checks pass.
Local outputs: artifacts/reference/minimum-20261024/native-attempt2. The original
exporter is archived; failed attempt1 (pre-solve ImportError) remains preserved.
[Producer evidence](experts/capytaine_expert/minimum-producer-2026-10-24.md).

Native centroid quadrature's spurious hydrostatic cross-couplings are corrected
in artifacts/reference/minimum-20261024/native-gl2 using native GL2 and zero new
BEM solves. All8 independent analytic checks pass, fullK error4.37e-11; every
original frequency variable/coordinate is unchanged. Independent executed review
passes. Preserve the original centroid output as diagnosed historical evidence.

Stock WEC-Sim v7.1.0 imported-excitation cropping defect reproduced independently:
future/past/interior errors [2,16,0]. Multiple-import omission independently
reproduced with error32. Both generic corrections passed independent executed
review: scalar dt1/.5, boundary/interior, insufficient support, multiple headings,
full6DOF and asymmetric lag checks. Maintained fork project branch is published
at 2f43f547301b223cd9a07e84dc88c9375249af0e; forkpin records identity.
[Startup review](experts/numerical_methods_researcher/startup-review-2026-10-24.md).
Finite-band packet tails must meet the prospectively frozen numerical-calm and
arrival gates in the [shared contract](experts/shared-contract.md); do not clip
signals or crop body startup. Official BEMIO import/H5 write passed with exact
metadata and phase residual1.455e-11; producer jointly approved stored import
normalization/CG/CB/fullM/K/Ainf. Runtime H5/body inertia mapping remains pending.
No motion/video result exists yet.

Next priority: incident-only startup diagnostics (single prospective spectral
packet proposal; direct spectral and retained-IRF initial force checks), focused
Simscape added-inertia product ordering reproduction, then ten-second WEC-Sim
motion/video with runtime mapping signoff. No full library or coefficient campaign
substitutes for this milestone. Packet/model files remain draft and unrun.

## Execution and ownership

[Execution state](experts/execution-state.json) owns numerical liveness. Producer
session84329 and hydrostatics correction were verified terminal and leases
released. Consumer fixed regressions and official BEMIO jobs are terminal exit0;
consumer owns the next serial packet diagnostic and focused inertia reproduction.
Recovered shared MATLAB engine9168/launcher37656 is the latest observed identity,
not a durable liveness claim. Reconcile actual process
identity and follow [claim/recovery](experts/coordination.md) before every launch.
Never infer terminal status from a timeout or an idle record alone.

Active branch baseline/fluid-body; preserve other worktrees and artifacts.
Coordinator owns handoff/assignments/integration, Capytaine owns producer,
WEC-Sim owns BEMIO/case/fork, researcher independently reviews startup and physical
corrections. Native roles are available; no activation campaign is needed.
Maintained fork/pin: [wecsim-fork.json](experts/wecsim-fork.json).
Two exercised-path corrections have passed numerical/independent review; four
historical marine-gyro leads remain separately recorded, without blind patches.

Earlier quadrature session36762 remains user-interrupted terminal exit1; partial
data is historical, not a task to resume or a numerical failure to relabel.
