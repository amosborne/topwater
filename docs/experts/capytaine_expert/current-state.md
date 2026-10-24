# Current state: Capytaine expert

Setup and read-only producer readiness inspection are complete. Reference work
is queued, not authorized; follow the [handoff](../../checkpoint-handoff.md).
CAP-001 scope and acceptance are in [assignments](../assignments.md). BEMIO
belongs to wecsim_expert. Numerical scheduling follows coordination.md.

Recorded setup: Capytaine 3.0.0 in the project .venv. Neutral input is
reference/capytaine/canonical-boat.json; geometry, wet17, dry mass and CG were
checked against the frozen input before e972c21. Recheck relevant versions and
hashes before production. No accepted hydro output, motion result or video exists.

Preserve the unfinished, uncommitted reference/capytaine/export_case.py. It was
never launched; incident_packet/export_controls are missing and unapproved.
Producer frequency/heading coverage, hydrostatic layout, full inertia/frame
mapping and transient arrival/timeline need joint consumer agreement before use.
No automatic mesh or coefficient campaign.

[Setup/interface evidence and source map](native-setup-2026-10-24.md) records
exact paths, versions, hashes and limits. Native exports use Mesh/FloatingBody,
RadiationProblem/DiffractionProblem, io.xarray.export_dataset and
io.legacy.export_hydrostatics. Prior studies remain historical, not evidence of
six-motion accuracy. Next: await user-authorized reference work and a bounded
coordinator assignment; do not repeat completed role confirmation.
