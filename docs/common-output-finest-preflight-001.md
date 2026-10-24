# Finest-source resource preflight accepts a short cost screen only

The read-only preflight at source `0b8dd9a` passes its prospectively declared
resource gates. It does NOT authorize the full field run or qualify any source/
reference result. No native evaluation or BEM solve ran; complete=true,
passed=false, resource_preflight_accepted=true, full_run_authorized=false.
Exact command chunk6da3b5 is terminal exit0, elapsed4.682s below60s preflight cap.

The [protocol](common-output-preflight-protocol.md) and
[controls](../validation/common-output-preflight-controls.json) declare a future
OFFLINE full-job ceiling of3600s,2GiB estimated peak memory and2GiB additional
disk. These prospective resource limits do not change any existing physical
gate, prior result, or Topwater's200wall seconds per simulated second requirement.

| Complete future work | Count |
|---|---:|
|Common targets at orders12/16/20|124416 /221184 /345600|
|All field source-target pairs|2,388,787,200|
|Native calibration pairs|11,943,936|
|All pairs|2,400,731,136|
|100-target native blocks, including calibration|6948|
|New six-mode complex fields, including calibration|66,686,976 bytes|

Every measured54/216 per-order rate and the independently measured864/order20
rate was retained. Complete-work estimates range2173–2578s before safety and
preparation margins. Upper estimate=max(rate projections)×1.25+60s=3282.320s,
below3600s. This linear-work estimate is not measured finest-source throughput.
Order20 alone would omit roughly half the field work.

| Resource | Model | Prospective ceiling |
|---|---:|---:|
|Peak memory|1,325,831,756 bytes (1.23GiB)|2GiB|
|Additional disk|1,847,910,862 bytes (1.72GiB)|2GiB|

The memory model includes273,656,016 resident input-array bytes, two177,265,854
byte NPZ verification buffers, new fields/native points, largest load buffers,
both complex S AND K allocations for100×3456 pairs with4×scratch allowance,
and512MiB runtime overhead. The disk model includes242,663,160 input-bundle bytes,
immutable/checkpoint/audit/packaging copies and actual producer bytes. Both are
explicit engineering models, not measured RSS or guaranteed upper bounds.

Recorded host availability was15,151,890,432 physical-memory bytes and
202,161,901,568 disk bytes. Required1.25×memory/2×disk availability margins passed.
Those observations expire; the native screen must recheck host resources and
the inspected Capytaine allocation-path hashes before evaluation.

All source/common-target/physical/input identities match their SHA-pinned
archived graphs. Installed Capytaine3.0 allocation-path Python hashes are
recorded without importing or executing Capytaine. Independent review reproduced
every work/rate/resource/host gate and provenance check. All278 tests pass;
exact test session9837 is terminal exit0, JUnit included. Fresh evidence
restoration verified the complete graph and all four logical input bundles.
No native numerical job is live. Source/full-reference qualification remain false.

[Preflight producer/data](../validation/evidence/common-output-finest-preflight-001.zip),
[summary and reconstructed model inputs](../validation/evidence/common-output-finest-preflight-001-summary.json).
SHA256 `11ea95f6e3c3a7a0249521b82f734bde21cb85cb679898f9e97c6f5d80dbcb3f`.
Required external tracked inputs are common-output-sources-001.zip,
frozen-source-gauss-continuation-001.zip and frozen-source-gauss-001.zip; their
SHAs are in the summary. After extraction run `python restore_inputs.py
--common-bundle PATH --quadrature-bundle PATH --historical-bundle PATH`.
The paired-study input is restored exactly from the historical bundle.

Reproduce with a fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.common_output_preflight --common artifacts/validation/common-output-sources-001.json --common-bundle validation/evidence/common-output-sources-001.zip --study artifacts/validation/paired-body-study-002.json --study-bundle validation/evidence/paired-body-001.zip --quadrature artifacts/validation/frozen-source-gauss-continuation-001.json --quadrature-bundle validation/evidence/frozen-source-gauss-continuation-001.zip --historical-bundle validation/evidence/frozen-source-gauss-001.zip --output artifacts/validation/common-output-finest-preflight-next.json
```

Next bounded goal: separately declare a representative3456-source native cost
screen,60s total including preparation. Repeat full native centroid potential/
all36-load calibration first at unchanged1e-10 scaled tolerance. Evaluate actual
full-source fields on deterministic100-target blocks at first/middle/last
positions of EACH reflected half for EACH common order; verify physical-facet
coverage and add deterministic blocks if needed, at most24blocks/order. Keep all
sampled points/owners/weights/fields and exact timing data. No mirrored fields.
Record actual process peak working set and committed memory, and refresh host/dependency
identity gates. No source solve, remeshing, optimization or full field job.

Before launch declare the complete probe/sampling/projection policy. Require
native/finite/ownership controls,60s budget, actual peak-memory/commit within2GiB,
refreshed resource margins and a conservatively projected complete job within
3600s before accepting feasibility. A pass only supports a separately reviewed
full-job decision; never auto-launch. Deliver immutable evidence, independent
review and checkpoint/tracker/handoff, then stop. Full4level source convergence,
finite A/B/excitation/phase and all physical gates remain unresolved.
