# Current state: Capytaine expert

CAP-001 minimum producer is authorized by the 2026-10-24 reference goal; follow
the [handoff](../../checkpoint-handoff.md). BEMIO and incident-packet/startup
generation belong to wecsim_expert. Numerical scheduling follows coordination.md.

Recorded setup: Capytaine 3.0.0 in project .venv. Physical canonical input
reference/capytaine/canonical-boat.json remains unchanged (SHA256
ea18844904e03bbfd72abc3396240c97b1994b3584e36c8149f8ae5f60a99e4f).
Native export completed; no accepted motion result or video exists yet.

Original unfinished exporter preserved at
artifacts/preserved-producers/export_case-a1feb6e8581695c1.py with original SHA256
a1feb6e8581695c1a7562bf1e129bbd609f96b54ab7dcac25da9b5ecdaa92942.
Active exporter is hydrodynamics-only; consumer owns forcing generation.
Consumer approved reference/capytaine/interface-contract.json: omega0.1:0.1:12
rad/s, propagation headings0/45/90deg, 17 panels, 1086 solves, native root
netCDF/legacy files plus full dimensional M/K/Ainf and exact CG/CB/volume in
native.mat/.npz. Official legacy text rounds metadata; retain exact auxiliary
values downstream. This is producer signoff only; startup/motion remains pending.

First allocated attempt (artifacts/reference/minimum-20261024/native) terminated
exit1 before solves: prototype imported Delhommeau from an unexported package
namespace. Corrected to capytaine.green_functions.delhommeau after installed-source
inspection. Failed manifest/log preserved. Attempt2 completed terminal exit0,
1086solves in69.9800752s with one thread. All14 finite/identity/metadata checks passed
in artifacts/reference/minimum-20261024/native-attempt2/producer-checks.json.
Ledger was set idle under lease before final release. Source/input/dependency
identity and command are preserved in producer-run-proposal-attempt2.json.

Physical limitation requiring resolution before motion reliance: native one-centroid
triangle quadrature integrates degree2 CB and hydrostatic integrands inexactly.
It produces CBxy=[-0.0130379264,0.00111213018] and spurious roll/pitch hydrostatic
coupling19607.59 despite the symmetric physical surface. Source evidence and
proposed exact hydrostatic-only triangle integration are in
[minimum producer evidence](minimum-producer-2026-10-24.md). No output was altered
or symmetry imposed. Await coordinator scope and independent correction review.
High-frequency native wavelength/resolution/irregular-frequency warnings also
remain diagnostic limits; no automatic refinement/coefficient campaign.

Historical [setup/interface evidence](native-setup-2026-10-24.md) is background,
not a motion validation. Next: resolve bounded hydrostatic quadrature issue if
assigned and independently review consumer matrix mapping when ready. Do not
repeat role confirmation or launch additional numerical work without allocation.
