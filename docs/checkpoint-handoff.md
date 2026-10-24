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
full6DOF and asymmetric lag checks. These corrections were first published
at 2f43f547301b223cd9a07e84dc88c9375249af0e; forkpin records current identity.
[Startup review](experts/numerical_methods_researcher/startup-review-2026-10-24.md).
Finite-band packet tails must meet the prospectively frozen numerical-calm and
arrival gates in the [shared contract](experts/shared-contract.md); do not clip
signals or crop body startup. Official BEMIO import/H5 write passed with exact
metadata and phase residual1.455e-11; producer jointly approved stored import
normalization/CG/CB/fullM/K/Ainf. Runtime H5/body inertia mapping remains pending.
No motion/video result exists yet.

Packet attempt2 is rejected prospectively: eta12.346um and fluid speed42.824um/s
pass, arrival92.2% passes, but initial heave .37135N and pitch .13014Nm exceed
the frozen .1 limits (direct spectral checks agree). Results preserved in
artifacts/reference/minimum-20261024/packet-attempt2. No body run is authorized
from that packet. Attempt1 failed report serialization; raw artifacts retained.

Fixed-order scalar screen [11,11.5,12,12.5]s selected first passing .25Hz packet
focus11.5s; arrival59.4% canonical. Fullphysical-hull-footprint supplemental
check and independent executed review PASS: eta3.3236um, fluidspeed25.4955um/s
over587412samples plus valid Taylor bounds, unchanged35/43.9823um gates.
Exactfieldidentity with packet-attempt3 is verified; accepted incident inputs
are frozen. This proves input startup only, not body rest or applied model force.

Actual H5/body/Simscape inertia reproduction confirmed xy/yz ordering error
(native tensor error1071.6434kgm2). Three proper-API library boundary reorderings
give roundoff errors (fullmatrix1.82e-12). Independent executed/semantic review
PASS for selected rigid path, including unchanged connectivity/scripts and
documented R2025b migration/defaults. Published fork now
fcd7c365a1201bbb56ce7f571a5baf9b2a94d89b; forkpin owns exact dependency status.
Preserve original coefficients and nonsymmetric residual, never zero couplings.

Case/frame source review passes (CG/STL once, RzRyRx radians, worldomega,
all6DOFfree/noPTOmoor/empiricaldamping, rest/memory reset). Body remains HELD:
actual parent excitation mask clones a branch per importedwave, each consuming
the same precomputed total. Without correction it would apply3*total. Consumer
owns bounded actualmask reproduction, properAPI homogeneous-import singlebranch
correction and fresh/repeated/reconfiguration regressions; researcher independently
reviews. Prior preprocessing tests did not cover that compiledmask integration.

Next: resolve mask duplication with independent executed review, then first
ten-second WEC-Sim motion/runtime mapping checks and isometricvideo. No full library
or coefficient campaign substitutes for this milestone. Model remains unrun.

## Execution and ownership

[Execution state](experts/execution-state.json) owns numerical liveness. Producer
session84329 and hydrostatics correction were verified terminal and leases
released. Consumer fixed regressions and official BEMIO jobs are terminal exit0;
consumer owns the next serial actualmask reproduction/correction/regressions,
then conditionally authorized firstbodyrun after all independent gates pass.
Recovered shared MATLAB engine9168/launcher37656 is the latest observed identity,
not a durable liveness claim. Reconcile actual process
identity and follow [claim/recovery](experts/coordination.md) before every launch.
Never infer terminal status from a timeout or an idle record alone.

Active branch baseline/fluid-body; preserve other worktrees and artifacts.
Coordinator owns handoff/assignments/integration, Capytaine owns producer,
WEC-Sim owns BEMIO/case/fork, researcher independently reviews startup and physical
corrections. Native roles are available; no activation campaign is needed.
Maintained fork/pin: [wecsim-fork.json](experts/wecsim-fork.json).
Three exercised-path corrections have passed numerical/independent review; four
historical marine-gyro leads remain separately recorded, without blind patches.

Earlier quadrature session36762 remains user-interrupted terminal exit1; partial
data is historical, not a task to resume or a numerical failure to relabel.
