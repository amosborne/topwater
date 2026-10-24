# Full-output consumer rehearsal passes; stop the rejected family's long extension

At source `d0786ad`, the actual minimal output consumer completed its resource
rehearsal in24.503s below60s. Exact session15358 is terminal exit0. All691200×6
complex synthetic field values were allocated, touched and retained together:
66,355,200 bytes plus331,776 native-field bytes. Native3456 potentials/all36
loads again reproduced at unchanged1e-10 scaled tolerance. No new BEM solve or
long field evaluation ran. Synthetic fields are RESOURCE-ONLY, never reference
potentials; complete=true, passed=false, rehearsal_passed=true, resource_only=true,
synthetic_fields=true, full_run_authorized=false. Physical qualification staysfalse.

The [prospective protocol](finest-consumer-rehearsal-protocol.md) and
[controls](../validation/finest-consumer-rehearsal-controls.json) were reviewed
before launch. One review blocker was corrected before launch: the generic
checkpoint helper read whole NPZ bytes for hashing. The consumer now uses the
same streaming checkpoint writer intended for actual fields. Historical payload
SHA/producer checks precede selected-array loading; unneeded old fields remain
in pinned archives rather than being deleted or relabeled.

All eight required stages passed: inputs, native calibration, checkpoint12,
checkpoint16, checkpoint20, final checkpoint, array verification and package
verification. The exact consumer keeps caller-owned complex buffers, uses the
full-size pressure-load calculation, saves actual NPZ checkpoints, streams raw
hashes, compares restored arrays one at a time, and writes/verifies actual ZIP
CRC, metadata snapshot and streamed embedded NPZ SHA. High-entropy deterministic
real/imaginary data exercise retained memory and pessimistic packaging size.
The immutable final NPZ SHA does not change during metadata-only stage updates.

| Measured resource | Result | Declared ceiling |
|---|---:|---:|
|Peak working set|533,397,504 bytes|2GiB|
|Peak committed memory|1,926,377,472 bytes|2GiB|
|Accounted additional disk, including64KiB metadata allowance|434,347,893 bytes|2GiB|

Peak committed memory occurred during input loading and did not grow during
the retained-field, load, save, verification or packaging stages. Current private
memory grew from1,565,564,928 input bytes to1,677,119,488 final bytes. All stage
counters are retained. This closes the specific full-output lifetime concern
for this consumer and its tested counts; it is not a numerical accuracy result
or a universal memory bound for another workload.

Independent postrun review reproduced every seeded field value, target identity/
ownership/normal/weight, all three synthetic load matrices, native calibration,
producer/input/dependency/native-library hashes, stage order and resource/file
checks. All284 tests pass; exact test session84363 is terminal exit0, JUnit
included. Fresh restoration verified the complete graph and exact inputs. No
native job remains live. Offline costs never establish Topwater200x throughput.

[Complete producer/resource data](../validation/evidence/finest-consumer-rehearsal-001.zip),
[summary and all36 necessary-magnitude diagnostics](../validation/evidence/finest-consumer-rehearsal-001-summary.json).
SHA256 `74dc090f8b4dfc1a5aaedd1c5340f4342e1e3ae2b2e76862c9769493564bb657`.
The immutable payload ZIP is included once; it contains the exact NPZ and actual
producer graph. Required tracked external inputs are common-output-sources-001.zip,
frozen-source-gauss-continuation-001.zip and frozen-source-gauss-001.zip. After
extraction run `python restore_inputs.py --common-bundle PATH --quadrature-bundle
PATH --historical-bundle PATH`. This restores source/NPZ paths with SHA checks;
all synthetic/resource-only flags remain intact.

Reproduce with a fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.finest_consumer --probe artifacts/validation/finest-cost-probe-001.json --probe-bundle validation/evidence/finest-cost-probe-001.zip --output artifacts/validation/finest-consumer-rehearsal-next.json
```

Do NOT spend the projected52-minute3456-source full field run merely to seek
qualification of this paired uniform family. The already failed216→864 increment
would still be one of the two required finest increments. Under the unchanged
charged0.001+1e-6 policy, that increment alone requires hypothetical finest
|heave added mass|>=457,647kg and |pitch added inertia|>=922,284kgm²; current864
values are8246.6kg and16275.4kgm². All36 necessary (not sufficient) magnitudes are
derived in the archived packaging program/summary from the exact common-output
increments/charges/scales. A normal finer result cannot plausibly rescue this
qualification strategy. This is an engineering decision to stop expensive work,
not an extrapolation, continuum bound, new physical pass or claim that slivers
have been proven causal. Existing native four-level and controlled three-source
failures remain failed. Generic mathematical coverage and all data are preserved.

Next bounded goal: READ-ONLY asset/physical/mesh-quality preflight of the existing
fixed-seed study003 family on the SAME canonical boat. Inspect actual archives
reference-study-003.zip, frozen-fixed-001.zip and periodic-sources-003.zip with
their manifests/raw producer graphs. Verify original physical identity, all four
meshes/facets/moments/conformity, radius/vertex-order controls and quality before
claiming this is a better candidate. Inventory reusable finite-frequency sources
separately from absent A∞ sources, which would need fresh prospectively declared
export. Assess a common physical output grid and complete resource/uncertainty
requirements without any Green evaluation or source solve. Declare60s read-only
budget/2GiB memory and additional-disk ceilings before execution. Deliver the
inventory, quality/cost decision, independent review, checkpoint/tracker/handoff,
then stop. No silent geometry/tolerance changes, new formulation, runtime edit,
remeshing or optimization campaign. Old study003 gates stayfailed; all36 source
and finite A/B/excitation/phase gates still require qualification.
