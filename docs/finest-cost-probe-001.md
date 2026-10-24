# Finest-source timing probe passes; full-output memory still needs rehearsal

At source `ebaf281`, the bounded native probe completed in21.056s below its
declared60s cap. Exact session4378 is terminal exit0. Native potential and all36
loads reproduced at unchanged1e-10 scaled tolerance; tolerance ratios are
1.95418e-5 and1.69874e-6 (ratios OF that tolerance, not physical relative errors).
All18 sampled100-target blocks cover all seven physical facets, both reflected
halves and each common order12/16/20.5256 points were newly evaluated including
3456 native calibration points; zero BEM solves ran.

The [prospective protocol](finest-cost-probe-protocol.md) and
[controls](../validation/finest-cost-probe-controls.json) were reviewed and
committed before launch. Physical facets, source/target/CG/inertia/Green identity,
owned interior points, input/producer/installed dependency and native-binary
hashes passed. Both halves used independent full-source evaluation. Raw indices,
points/owners/weights/facet IDs, six complex fields and times are retained.

Native calibration took10.401s; preparation4.754s. Sample blocks took0.286–0.349s.
The declared maximum sampled pair rate projects2414.747s of COMPLETE future
field work. Including native calibration,1.25 safety and60s reserve gives
3091.435s, below the prospective3600s offline ceiling. This is a sample-based
projection, not measured full-job timing or Topwater simulation throughput.

| Measured process resource | Bytes | GiB |
|---|---:|---:|
|Peak working set|532,029,440|0.50|
|Peak committed memory|1,925,001,216|1.79|

Both actual probe caps passed at2GiB. However, measured committed memory exceeds
the earlier engineering byte estimate1,325,831,756 by599,169,460 bytes. Its peak
occurred during input loading. Current private memory was1,698,713,600 bytes at
the input snapshot and1,764,352,000 bytes at completion. The probe does not retain
all future output fields or exercise full-size save/verification stages.

Adding future fields66,686,976 bytes, largest load workspace49,766,400 bytes and
64MiB reserve to measured peak gives2,108,563,456 bytes: only38,920,192 bytes
(37.1MiB) below2GiB. This candidate additive estimate is not sufficient full-job
memory evidence; newly introduced verification/serialization copies could exceed
it. Independent review requires a bounded rehearsal of the ACTUAL minimal
consumer and its output stages before the long field run. Do not load an extra
95MB quadrature-array graph already represented by SHA-bound common-grid inputs.
No memory ceiling or physical gate was relaxed to hide this finding.

The declared probe itself is feasible: complete=true, passed=false,
probe_feasible=true, full_run_authorized=false. Source/full-reference qualification
remain false and all old failures remain failed. All280 tests pass; exact test
session77091 is terminal exit0, JUnit included. Independent postrun review
reconstructed native loads, sample geometry/coverage, identities, counters and
projection. Fresh archive restoration verified the complete graph and inputs.
No native numerical job remains live.

[Complete producer/data](../validation/evidence/finest-cost-probe-001.zip),
[summary](../validation/evidence/finest-cost-probe-001-summary.json).
SHA256 `cd5756817dcd38a99ae4931fa9f3befe5c8ba733a4cb1ace74cbe9852bb084cc`.
Exact preflight archive is embedded. Required tracked external inputs are
common-output-sources-001.zip, frozen-source-gauss-continuation-001.zip and
frozen-source-gauss-001.zip; SHAs are in the summary. After extraction run
`python restore_inputs.py --common-bundle PATH --quadrature-bundle PATH
--historical-bundle PATH` to restore all four exact logical input bundles.

Reproduce with a fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.finest_cost_probe --preflight artifacts/validation/common-output-finest-preflight-001.json --preflight-bundle validation/evidence/common-output-finest-preflight-001.zip --output artifacts/validation/finest-cost-probe-next.json
```

Next bounded goal: implement the actual minimal common-grid finest-source
consumer with a separately declared60s output-memory rehearsal. Use the SAME
input/import path and array lifetimes intended for the real run; allocate AND
touch all full six-mode output and load buffers, then exercise actual checkpoint
save/hash, verification and packaging while monitoring peak working set/commit.
Refresh host/dependency identities and retain exact immutable producers/inputs.
Keep2GiB process and2GiB additional-disk caps; failed/missing stages cannot pass.
Declare any small native warmup/calibration workload before launch. Rehearsal
fields must be marked resource-only and can never pass a physical gate. No full
field run, source solve, remeshing, new physics or optimization campaign.
Deliver actual lifetime/stage memory evidence, independent review and checkpoint,
then make the separately reviewed full-job decision. All four source levels,
finite A/B/excitation/phase and full-reference physical gates remain required.
