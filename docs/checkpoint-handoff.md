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

Physical limitation: native centroid quadrature gives spurious hydrostatic
cross-couplings on the asymmetric triangle partition of the symmetric surface.
Producer owns a bounded exact degree2 hydrostatics-only correction with independent
researcher review, retaining the same hydrodynamic coefficients and previous
outputs. Do not treat the centroid hydrostatics as accepted motion physics.

Stock WEC-Sim v7.1.0 imported-excitation cropping defect reproduced independently:
future/past/interior errors [2,16,0]. A generic padding correction is authored;
fixed regressions and final independent review are pending. Multiple-import
omission is another focused reproduction, not yet a confirmed correction.
[Startup review](experts/numerical_methods_researcher/startup-review-2026-10-24.md).
Finite-band packet tails must meet the prospectively frozen numerical-calm and
arrival gates in the [shared contract](experts/shared-contract.md); do not clip
signals or crop body startup. No motion/video result exists yet.

Next priority: finish consumer reproductions/corrections, reviewed hydrostatics,
then official BEMIO import, incident-only startup diagnostics, and ten-second
WEC-Sim motion/video with executed producer/consumer signoff. No full library
or coefficient campaign substitutes for this milestone.

## Execution and ownership

[Execution state](experts/execution-state.json) owns numerical liveness. Producer
session84329 was verified terminal and its lease released. Consumer owns the next
serial stock-multiple-import/fixed-regression sequence; engine recovery may be
needed because the prior shared MATLAB engine exited. Reconcile actual process
identity and follow [claim/recovery](experts/coordination.md) before every launch.
Never infer terminal status from a timeout or an idle record alone.

Active branch baseline/fluid-body; preserve other worktrees and artifacts.
Coordinator owns handoff/assignments/integration, Capytaine owns producer,
WEC-Sim owns BEMIO/case/fork, researcher independently reviews startup and physical
corrections. Native roles are available; no activation campaign is needed.
Maintained fork/pin: [wecsim-fork.json](experts/wecsim-fork.json). No engine
correction has passed numerical review at this checkpoint.

Earlier quadrature session36762 remains user-interrupted terminal exit1; partial
data is historical, not a task to resume or a numerical failure to relabel.
