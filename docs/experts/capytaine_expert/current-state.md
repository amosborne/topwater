# Current state: Capytaine expert

Native custom-role activation confirmed 2026-10-24; delivered developer
instructions match the producer role. Setup/read-only readiness inspection is
complete; see [native setup and interface note](native-setup-2026-10-24.md) for
exact paths, versions, hashes and blockers. No numerical lease or jobs. Exporter
and canonical input preserved. Missing packet/export controls and consumer
hydrostatic-file layout/startup timeline agreement remain genuine prerequisites;
no completed hydrodynamic or motion result is claimed. Stop pending coordinator
assignment; BEMIO remains with wecsim_expert.

Confirmed: Capytaine3.0.0 in project analysis .venv. Neutral geometry and physical
inputs are in ../../../reference/capytaine/canonical-boat.json. Shared
contract and current10s motion/video goal supersede older coefficient-first
campaigns. All old failures remain historical; do not automatically resume them.

Source map: installed .venv/Lib/site-packages/capytaine; native Mesh/FloatingBody,
RadiationProblem/DiffractionProblem, io.xarray export_dataset and io.legacy
export_hydrostatics. Pinned WEC-Sim upstream example
artifacts/dependencies/wec-sim-v7.1.0/examples/BEMIO/CAPYTAINE demonstrates netCDF
plus hydrostatic-file export. WEC-Sim/BEMIO details belong to wecsim_expert.

Reproduce physical input without Topwater imports: parse the neutral JSON.
Its full hull/wet17, dry M and CG were verified exactly against the existing
frozen input before e972c21. Existing reference/capytaine/export_case.py is
unfinished and never launched: retain it, review before use, do not treat its
missing incident_packet/export_controls fields as frozen or accepted.

Known evidence: docs/capytaine-wecsim-workflow.md and setup JSON; old studies
under validation/evidence and resultdocs. Coefficient reproduction is not proof
of six-motion accuracy. No active boat result or video exists.

Current bounded assignment: confirm this role and ownership, review the proposed
producer interface READ-ONLY, and return the minimum export/frequency/heading/
hydrostatic requirements for the first minimum-mesh WEC-Sim case. Coordinate
with the WEC-Sim expert. Exclude hydrodynamic jobs, new mesh studies, source edits
and Topwater code. Deliver a concise readiness/interface note and update this
knowledge base; stop pending execution-owner assignment. No numerical lease.
