# Owned bounded assignments

Task statuses are `queued` (owned but not authorized), `authorized` (within the
current user request), `running`, or `completed`. Current project authorization
is recorded only in [the handoff](../checkpoint-handoff.md); an execution lease
schedules authorized work and never grants scope. Update status in place.

## CAP-001: first canonical producer/interface

Status: queued. Owner capytaine_expert; consumer partner wecsim_expert. Readiness inspection
completed in capytaine_expert/native-setup-2026-10-24.md. Subsequent execution is
queued, with no current numerical lease.

Scope: make the minimum-mesh native Capytaine export needed by the first WEC-Sim
case, preserving exact physical surface/full inertia/CG and agreed DOF/heading/
phase/units. Coordinate frequency and transient coverage; retain native files
and exact identity. Exclude Topwater code, mesh campaign, intermediate coefficient
acceptance thresholds, and unauthorized fixes to the consumer.

Deliver compatible native files, small jointly approved producer/consumer
contract and reproducible command/manifest. Acceptance: consumer can read the
physical and numerical data without missing/defaulted metadata, phase or unit
ambiguity; no coefficient threshold is required. Freeze controls before launch.
Stop at reviewed minimum-case export/interface; no automatic refinement.

## WEC-001: maintained fork, independent leads and exercised paths

Status: queued. Owner wecsim_expert. Fork/source-map/setup inspection completed
in wecsim_expert/native-setup-2026-10-24.md. Remaining reproductions await
user-authorized reference work and a bounded coordinator assignment, then a
sole numerical execution lease.

Scope: pinned v7.1.0 maintained fork; independently reproduce the four specific
marine-gyro leads (yaw-order/frame and velocity/acceleration mapping, unequal
postprocessing grids, angular-bin energy normalization, direction row/column
shape) using independent expectations. Inspect the canonical producer/BEMIO/
startup/radiation/hydrostatic/frame/output/video paths actually exercised.
Review source map first; choose the smallest decisive reproductions. No full
library audit, automatic carry-over patches, string-replaced generated models,
Topwater changes, broad refactor or upstream posting/PR without authorization.

Deliver per lead: selected/affected versions, exact paths, minimal frozen
reproduction, independent expected versus actual, physical consequence, cause
confidence, confirmed/unaffected/unconfirmed status, focused regression and
minimal generic fork correction ONLY if confirmed. Case-specific configuration
stays here. Missing evidence is not a confirmed defect. Record exact pins and
remaining limitations in the KB; concise report to coordinator and Capytaine
partner. Independent reviewer distinct from author before physical reliance.

Stop at bounded lead/exercised-path findings and reviewed minimal corrections
or explicit limitations. Reproduction commands/resources/expectations must be
declared before their execution lease. Result-invalidating confirmed issues
must be corrected or the affected result withheld. Speculative unrelated issues
must not delay the first working minimum-mesh10s simulation indefinitely.

## TOP-001: held-state ownership

Status: completed (setup only; implementation held). Owner topwater_expert.
Activation/documented equation/code/evidence map completed
in topwater_expert/held-source-map.md. No solver implementation, simulations,
test runs, optimization or method experiments. Stop pending future user authority.

## NUM-001: focused research support

Status: completed (setup only). Owner numerical_methods_researcher. Proposal template and targeted initial
wave/IRF-padding question completed in its KB. No broad literature/optimization
survey or numerical lease. Next only a focused question requested by reference
owners; proposed experiments go through coordinator selection and implementation
owner. Stop at supported answer and cheapest rejection plan, not a new solver.

## Downstream milestone (not this setup's completion criterion)

The first freely floating WEC-Sim ten-second incoming-wave motion and isometric
video, followed by motion variation with mesh. Keep shared forcing/initial state
and exact producer identity; no PTO/mooring/invented constraints. Coefficients
remain diagnostics. Setup, fork creation and role confirmations do not claim
that simulation/video or the four defect reproductions already exist.
