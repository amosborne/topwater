# CANCELLED / UNRUN: incremental official CIC example diagnostic

LATEST steering cancels this diagnostic path: no PTO-zero execution. Source and
readiness artifacts below are preserved UNRUN. Current goal is the milder buoy
physicalincoming train/full60s/video.

Historical requested priority: progressively modify the actual passing WEC-Sim provided
example toward our configuration. No new mesh study, coefficient repair, custom
sea search or root-cause assumption. Coordinator owns each next allocation.
This is a bounded PLAN and ONE first-step source/readiness proposal, not a suite.
No new simulation or model modification has been performed in this assignment.

## Starting point and first cheapest meaningful change

Actual passed official v7.1.0 `tests/RegressionTests/RegularWaves/regularCIC`,
maintained pinf962103a/basea23d299, supplied RM3 H5/STL/model and two heave-only
bodies, regularH2.5m/T8s,200s/dt.1/ramp100/CIC20/stateSpace0/b2b0. Existing
full-history supplied-golden heave errors <=4.16e-13m versus official1e-5m;
native forces/grid/initialv+Frad0 pass. [Official receipt](official-examples-readiness-2026-10-24.md).
This exercises two restricted heave channels, NOT arbitrary6DOF radiation coupling.

STEP1 ONLY: change `pto(1).damping` from1,200,000 to0Ns/m using the existing native
MCR interface. Keep PTO stiffness0, actual translational PTO joint, outer
Translational constraint, original waves/ramp/bodymass/inertia/H5/kernel/dt/duration,
zero initial velocity and native zero CIC history. Do not remove the PTO block:
zero damping removes mechanical dissipation while retaining the exact kinematics.
This directly tests whether the substantial provided PTO dissipation masks a
problem under the SAME excitation and radiation consumption, without exchanging
our boat coefficients or changing startup simultaneously.

There is no supplied golden trajectory for changed physics. Do not apply the
baseline1e-5m equality test to this control. Compare all2001 native state/force
histories, peaks/drift/late response and measured work scales with the preserved
baseline. A larger response is not automatically an instability or energy bug.

## Measurements, signs and open mechanical budget

Save native bodies/PTO output plus original APPLIED addedmass/total force from
body.hydroForce.hf1.storage before assertions. WEC-Sim postprocessing replaces
reported addedmass/total force; do not use that reconstructed force as the applied
residual acceleration-feedback term. Actual S=adjustedMass in heave remains in
hf.adjustedMass; hf.mass is dry/equilibrium mass. Both examples use identical
geometry/M/H5, so those same S/m/K properties also apply to the stored baseline.

Report in J/W, for BOTH passed baseline and PTO0 control:

- Incident force work Wexc=int sum(v*Fexc)dt and radiation resisting work
  Wrad=int sum(v*Frad)dt, cumulative minima/total and instantaneous minima.
- Applied residual addedmass work War=int sum(v*FAr_applied)dt, separate from
  the transferred-addedmass kinetic storage already included in the Solid.
- PTO absorption Wpto=int c*(vfloat_z-vspar_z)^2dt; also retain the actual signed
  PTO internal-mechanics power, plus/minus sign residuals and magnitude closure.
  PTO0 native internal-mechanics power must be exactlyzero. Ideal joint constraint
  work is zero only on its constrained motions: verify actual off-axis velocities
  remain negligible rather than assuming this is a free6DOF energy system.
- E_S+U = sum[.5*S*vz^2+.5*K33*(z-CGz)^2], explicitly decomposed into dry kinetic,
  transferred-addedmass kinetic and hydrostatic storage. No additional fullAinf
  kinetic term is added; that would doublecount the transferred subset.
- Mechanical budget residual R=Delta(E_S+U)-[Wexc-Wrad-War-Wother-Wpto]. Report
  full R(t), max/end J, actual incident-work scale and peak mechanical energy;
  also pointwise power closure sum(S*vz*az)-sum(v*Fhydro_applied)+c*vrel^2.

This is an OPEN forced mechanical budget. Frad work includes fluid memory exchange,
not independently measured farfield outgoing energy. Nonzero R includes sampled
work quadrature/native time integration and any unaccounted terms; report it,
do not declare energy conservation/neutrality or hide it behind bounded motions.
No arbitrary energy-relative threshold or eigenvalue tolerance is invented here.
If pointwise closure fails, resolve sign/quantity/constraint accounting BEFORE
attributing any defect. If integrated residual is large against work/storage
scales, report an unclosed budget and stop for a focused reviewer decision.

## Decision branches and possible later single-feature steps

If STEP1 fails native finite/grid/input/force/PTO0 checks, preserve all outputs and
stop; no automatic source correction/rerun. If mechanically consistent but strongly
growing, the contrast with the same damped baseline establishes damping sensitivity,
not a library/BEM culprit. Do not add damping to make it pass. If bounded, it only
shows this provided heave-only system operates without that PTO dissipation.

Subsequent changes are proposals for separate decisions, not preallocated jobs:

| Axis to isolate | Keep fixed and what it would exercise |
|---|---|
| Solver dt.1 to.02 | Explicitly keep CICdt=.1 first; isolate solver/Ar-feedback time discretization. Changing defaultdt and CICdt together would confound them. |
| CICdt.1 to.02 | Keep solver.02 and CIC horizon20; isolate convolution quadrature/RT/PCHIP sampling. |
| CIC horizon20 to10 | Keep temporal grids/coefficient dataset fixed; isolate finite-memory truncation. |
| Outer constraint heave to provided Floating3DOF, then6DOF | Keep relative translational PTO joint. Pair gains collective motion, but bodies are still kinematically tied in five relative coordinates. NOT two independently free6DOF bodies. |
| Relative PTO kinematics to free relative6DOF | PTO already has zero mechanical damping/stiffness. This separate topology change finally exposes independent body rotations/coupled self-radiation; source/frame/constraint-work checks needed first. |
| True calm incoming-wave startup | Stationary regular wave/ramp100 cannot be made a physical calm-arrival field by simply setting ramp0. Replace incident-field/startup/ramp as ONE physically coherent input package, source-reviewed with full-hull eta/velocity/F0 gates; retain v0/CIC0. No stationary-history workaround. |
| RM3 bodies/H5 to our boat | Replace geometry/fullM/CG/K/Ainf/radiation/excitation as ONE coherent producer-consumer package. Mixing an RM3 mass/H5 with our hull is not a meaningful one-coefficient experiment. Preserve matched axes/units/phase, actual loaded metadata and full6force/NE checks. |

No need to perform all listed steps. Choose the next cheapest contrast from actual
STEP1 evidence. b2b0 suppresses interbody hydrodynamic coupling in the official
case; enabling it would be a separate change and is not automatically needed for
our eventual one-body case. No six-DOF conclusion from heave-only success.

## First-step source/freeze and prospective command

NEW `run_official_pto_zero.m` SHAe3e15b252cdf438af67a4535d28db2ef8dec9213048c307261b947c8f87057f6;
NEW `audit_official_pto_energy.m` SHA126abfc69e64b64454753bc0485875260e32b18da62f13d690696b668d83eff9.
MATLAB syntax lint passes. `prepare_official_pto_zero.py` copies byte-identical
official assets into fresh `artifacts/reference/official-pto-zero-20261024/cic`;
physical input/model/H5 files unchanged, sole MCR physical override is PTO0.
Already declared display-only explorerOFF/invisible figures retained.

Prospective manifest `artifacts/reference/official-pto-zero-20261024/cic/proposal.json`,
SHAb84490356dd986c3a0b3cfca7a09525825a9690045b74b7db57704806c7df523.
Exact proposed producer:

```
python reference/wecsim/run_matlab.py reference/wecsim/run_official_pto_zero.m --log artifacts/reference/official-pto-zero-20261024/cic/engine.log
```

One existingMATLABengine/oneCPU/max600s, expected30-240s. Native wecSim entry,
original200s/2001grid, no new BEM/import model or custom integrator. Independent
source/math/freeze review required BEFORE allocation. After any future grant:
authoritative engine/model/process readiness, atomicclaim -> SEPARATEverify -> run;
exactterminal/frozeninputs -> idleUNDERlease -> SEPARATEreadback -> FINALrelease.
STOP after this ONE control. No automatic suite/nextDOF/input change or encoder.
