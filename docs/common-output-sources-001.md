# Common output grid confirms unresolved source discretization

At source `1005feb`, the existing 54/216/864-panel Capytaine source solutions
were evaluated on IDENTICAL physical pressure-integration targets, normals and
weights: retained 864-source Gauss12/16/20 grids. All three source solutions
pass their pressure-integration controls. Source convergence still fails all18
symmetry-allowed entries for both increments; charged maximum tolerance ratios
are7541.47/2656.03. The quadrature-robust source trend passes. No numerical
threshold was relaxed. This is a three-source diagnostic, not qualification or
replacement of the original failed four-level or finite-frequency study.

The [prospective protocol](common-output-source-protocol.md) declares the source
comparison and its explicitly charged integration uncertainty. For each source,
Uq=max(abs(A20-A16),abs(A16-A12)). Each source increment includes the sum of its
two adjacent Uq envelopes; the robust trend compares the upper latest increment
against the lower previous increment plus the original absolute floor. These
empirical envelopes are not certified continuum error bounds. All integration
controls must pass before the source-policy diagnostic can pass.

| Source panels | Quadrature previous max ratio | Latest max ratio | Integration |
|---:|---:|---:|---|
|54|0.013451|0.003595|pass|
|216|0.065159|0.021549|pass|
|864|0.246235|0.082461|pass|

| Diagonal | 54-source output | 216-source output | 864-source output |
|---|---:|---:|---:|
|Surge (kg)|90.795191|91.371704|91.570110|
|Sway (kg)|535.943152|566.871780|571.206087|
|Heave (kg)|9463.055173|8704.232501|8246.591058|
|Roll (kg m²)|710.361150|933.144655|982.408698|
|Pitch (kg m²)|5980.206276|15353.161464|16275.381398|
|Yaw (kg m²)|520.372604|1137.944973|1282.109775|

For example, the latest heave source change is -457.641442 kg with adjacent
integration charge0.007870 kg; pitch changes+922.219934 kg with charge0.073400.
Using common output nodes makes the unresolved source effect visible separately
from pressure-integration changes. It does not localize continuum error or prove
which source-mesh correction will resolve it. All36 signed/absolute source
increments, Us envelopes, separate Uq charges, ratios and failure masks are saved.
Full complex fields, all nine36-entry matrices, native controls and actual
source densities are retained. Reused864 fields remain bound to their actual
historical producer rather than being relabeled as newly computed.

Native calibration of54/216 reproduced all six potentials/all36 loads at the
unchanged1e-10 scaled tolerance. Complete per-facet physical surface, moments,
conformity, radius, common-target identity and independent reflected-half field
evaluation checks passed. Source, physical and Green identities match. No BEM
solve, remeshing, matrix symmetrization, fitted field or new formulation was used.

Exact session56519 is terminal exit1, complete=true/passed=false, elapsed185.941s
below the declared360s offline cap. All276 local tests pass; exact test session
28461 is terminal exit0, JUnit included. Independent prelaunch/postrun review
reconstructed all nine matrices, geometry/fields/provenance, uncertainty charges
and every gate/mask. Fresh archive restoration passed. No native job remains
live. Source/four-level/full-reference qualification remain false; offline cost
does not establish Topwater simulation speed or coupled-motion accuracy.

Evidence consists of ALL THREE tracked archives:

- [Common-output data/producers](../validation/evidence/common-output-sources-001.zip),
  SHA256 `c208dffe223dffb1db87af8b8bf0eabd96d8ad6a40b358e924c26208dd7fbef1`.
- [Quadrature input](../validation/evidence/frozen-source-gauss-continuation-001.zip),
  SHA256 `7e999761722f68b52c9675abfb32bc02aff101c55f03ad76babf749154ab7a63`.
- [Historical input](../validation/evidence/frozen-source-gauss-001.zip),
  SHA256 `28a35794d8c20d35e59f955e0b4ced8fec3cdebd4c6f8e01beffd1cb205fdbb0`.

After extracting the common-output archive, run `python restore_inputs.py
--quadrature-bundle PATH --historical-bundle PATH`. This SHA-verifies and
restores all three exact logical input bundles, including the paired-study
payload contained in the historical archive. The complete owned graph was
verified after fresh restoration. [Summary](../validation/evidence/common-output-sources-001-summary.json).

Reproduce with a fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.common_output_sources --study artifacts/validation/paired-body-study-002.json --study-bundle validation/evidence/paired-body-001.zip --quadrature artifacts/validation/frozen-source-gauss-continuation-001.json --quadrature-bundle validation/evidence/frozen-source-gauss-continuation-001.zip --historical-bundle validation/evidence/frozen-source-gauss-001.zip --output artifacts/validation/common-output-sources-next.json
```

Next bounded goal: resource preflight ONLY for the archived3456-source solution
on this same physical output grid. Include ALL124416+221184+345600 targets:
2,388,787,200 source-target pairs plus native calibration, preparation and
storage. Order20 alone is only half that work. Archive the exact graph/input
identities and calculate explicit memory/disk/time models against prospective
ceilings. A separately declared short representative native cost screen may
follow if needed, but no full field job launches automatically. Deliver the
preflight/review and accept-or-reject decision, then checkpoint and stop. The
measured864-source rate suggests roughly2500s for all three orders before
additional overhead; this is offline reference work, not a runtime projection.
Use that evidence to justify the next source-control choice. The fourth source
level, finite A/B, complex wave responses and all full-reference gates remain
required. The broader time-domain simulation goal remains unachieved.
