# Native activation and minimum-case export readiness, 2026-10-24

Native `capytaine_expert` developer instructions were delivered on activation,
including producer ownership, WEC-Sim ownership of BEMIO, no numerical execution
without coordinator lease, immutable running inputs, independent physical-result
review, and held Topwater implementation. Selected parent settings were inherited.
Durable context is these files, not permanent model learning.

Scope: read-only producer/interface inspection and this KB update. Acceptance is
an evidence-backed readiness note, not a completed export. Stop here pending an
explicit assignment and execution lease. No MATLAB, Python solver, tests,
hydrodynamic runs or mesh studies were launched; exporter and canonical input
were not modified.

## Source identity and inspected paths

Project HEAD: `e972c217577f0865b6315b5aa521c95cd326532e`.
Capytaine installed metadata explicitly reports **3.0.0**; source inspection did
not import or run it. WEC-Sim example checkout HEAD independently read as
`a23d2993807410fe78cc63546c1ba9d81912024f` (v7.1.0 per recorded environment).
`docs/capytaine-wecsim-environment.json` records MATLAB R2025b and prior library/
license checks; those checks were not rerun and are not boat-motion evidence.

Exact paths read (relative to project root):

- `.codex/agents/capytaine_expert.toml`
- `docs/experts/capytaine_expert/README.md`
- `docs/experts/capytaine_expert/current-state.md`
- `docs/experts/shared-contract.md`
- `docs/experts/coordination.md`
- `reference/capytaine/canonical-boat.json`
- `reference/capytaine/export_case.py`
- `docs/capytaine-wecsim-environment.json`
- `.venv/Lib/site-packages/capytaine-3.0.0.dist-info/METADATA` (header)
- `.venv/Lib/site-packages/capytaine/io/legacy.py` (export source sections)
- `.venv/Lib/site-packages/capytaine/io/xarray.py` (export source sections)
- `artifacts/dependencies/wec-sim-v7.1.0/examples/BEMIO/CAPYTAINE/sphere/main.py`
- `artifacts/dependencies/wec-sim-v7.1.0/examples/BEMIO/CAPYTAINE/sphere/bemio.m`

An initial guessed example filename `.../CAPYTAINE/capytaine_test.py` did not
exist; it was not read. File enumeration located the actual sphere example above.

SHA256 at inspection: exporter
`a1feb6e8581695c1a7562bf1e129bbd609f96b54ab7dcac25da9b5ecdaa92942`;
canonical input
`ea18844904e03bbfd72abc3396240c97b1994b3584e36c8149f8ae5f60a99e4f`.

## Producer interface for subsequent joint approval

The physical input separates full 24 ft hull/material geometry from 17 wet world
triangles. Keep equilibrium, world CG `[0,0,0.15660636967711072]` m, mass 2500 kg
and full dry 6x6 matrix, including roll/yaw entry -70.03519269408712 kg m².
Keep rho=1025 kg/m³, g=9.81 m/s², infinite depth and zero speed.
Modes are Surge/Sway/Heave/Roll/Pitch/Yaw, rotations about world CG; world axes
x forward, y port, z up. Propagation headings are CCW from +x: source radians,
canonical directions 0/45/90 degrees. Complex source convention is
Re(amplitude exp(-i omega t)). The consumer owns phase/normalization conversion;
do not add a second conjugation to the official importer conversion recorded in
the shared contract.

Native deliverables must include netCDF with complete six-mode A/B and complex
diffraction plus Froude-Krylov excitation at agreed frequencies/headings,
hydrostatics (`Hydrostatics.dat`, `KH.dat`), physical dry mass/inertia and CG/CB/
volume identities, exact mesh and input/producer/dependency hashes, versions and
reproducible command. Preserve A-infinity explicitly with agreed consumer mapping.
The installed export source splits complex netCDF fields and emits these native
hydrostatic filenames. The official sphere example exports both netCDF and
hydrostatics into the **same** outputs directory, then runs readCAPYTAINE,
radiationIRF, excitationIRF and writeBEMIOH5. That example's density/depth/
sampling and optional state-space fit are examples, not canonical requirements.
No BEMIO implementation ownership is assumed here.

## Genuine readiness blockers and limits

1. The unfinished exporter reads absent `incident_packet` and `export_controls`
   keys. Its packet call occurs before any solve, so it is not a launchable
   canonical producer as written. Do not invent accepted controls or run it.
2. Prototype netCDF is at output root, hydrostatic files under `hydrostatics/`;
   this differs from the official example layout. WEC-Sim expert must specify
   the actual readCAPYTAINE lookup/location contract before producer acceptance.
3. Canonical prospective sampling is omega=0.1:0.1:12 rad/s and three headings;
   sea carriers are 0.1/0.25/0.5 Hz. This is not evidence of adequate excitation
   IRF bandwidth, window/padding or ten-second startup accuracy. Agree actual
   packet spectrum, per-direction histories, convolution reference/timeline and
   A-infinity/full-inertia mapping with the consumer before a launch.
4. Compact upstream initial wave support does not itself prove initially zero
   excitation at a resting hull when an excitation IRF has negative-time support.
   Arrival physics and rendered incident surface must use the same timeline;
   preserve t=0 rest, no imposed ramp or global switch. No numerical or physical
   correctness claim is made from this read-only source inspection.

Next action: consumer confirms file layout and mapping, both experts freeze the
bounded first-case interface/startup decision, then coordinator assigns execution
owner and a separate authorized task. No coefficient gate or mesh study here.
