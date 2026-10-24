# Fixed-source pressure integration stabilizes; reference remains unqualified

The prospective order-20 continuation at source `81d42a4` passes the unchanged
all-36-entry pressure-integration stabilization gate. Both 12→16 and 16→20
increments pass, maximum tolerance ratios 0.246235 and 0.082461, and the trend
passes. The policy remains 0.001 relative +1e-6 normalized absolute. This is
conditional on the SAME frozen 864-panel source solution, not a continuum bound
or source-mesh/reference qualification. All earlier failed studies stay failed.

Exact native session 81910 is terminal exit 1: complete=true, passed=false,
quadrature_stabilization.passed=true. Overall passed remains false deliberately.
Elapsed 318.337 s is below the prospectively declared 400 s offline cap. No BEM
solve ran. Native centroid and 345600 order-20 targets were newly evaluated;
orders 12/16 were reused only after pinned historical bundle/producer/input and
exact regenerated target/owner/weight checks. Reconstructed reused-load ratios
are zero. All 692064 retained points have positive per-face area weights and
strict ownership. Both reflected halves were independently evaluated, with
geometry-only target reflection. Source, normal, sign and physical inputs stayed
unchanged. Native potential/all-36-load controls again pass 1e-10 scaled checks.

| Order-20 diagonal | Value | Units |
|---|---:|---|
| Surge | 91.570110 | kg |
| Sway | 571.206087 | kg |
| Heave | 8246.591058 | kg |
| Roll | 982.408698 | kg m² |
| Pitch | 16275.381398 | kg m² |
| Yaw | 1282.109775 | kg m² |

Measured pitch output change is -672.247631 kg m²; unresolved remainder is
+242.705250 kg m², telescoping exactly to native finer-minus-coarse change
-429.542381 kg m². The remainder still includes changed source/collocation
solution AND finer-native output quadrature, not pure source error. All complex,
signed and absolute entries, field arrays, masks and measured uncertainty
envelopes remain available. No assumed spectral rate or extrapolation was used.

The [prospective protocol](frozen-source-gauss-continuation-protocol.md) and
[controls](../validation/frozen-source-gauss-continuation-controls.json) were
reviewed and committed before launch. All 273 local tests passed, JUnit included;
exact test session 83436 is terminal exit 0. Independent prelaunch/postrun review
and fresh archive restoration cover the numerical identities and limited gate.
No live native jobs remain. This is offline reference cost, not simulation speed.

Evidence consists of BOTH tracked archives:

- [Continuation data and actual producers](../validation/evidence/frozen-source-gauss-continuation-001.zip),
  SHA256 `7e999761722f68b52c9675abfb32bc02aff101c55f03ad76babf749154ab7a63`.
- [Pinned historical input bundle](../validation/evidence/frozen-source-gauss-001.zip),
  SHA256 `28a35794d8c20d35e59f955e0b4ced8fec3cdebd4c6f8e01beffd1cb205fdbb0`.

After extracting the continuation, run `python restore_shared_inputs.py
--reuse-bundle PATH_TO_HISTORICAL_ZIP`. It verifies the historical archive SHA,
restores its exact paired-study parent and shared logical input paths. Both
complete archived producer graphs were verified after fresh restoration.
[Machine-readable summary](../validation/evidence/frozen-source-gauss-continuation-001-summary.json).

Reproduce with a fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.frozen_source_quadrature --controls validation/frozen-source-gauss-continuation-controls.json --study artifacts/validation/paired-body-study-002.json --bundle validation/evidence/paired-body-001.zip --reuse-study artifacts/validation/frozen-source-gauss-001.json --reuse-bundle validation/evidence/frozen-source-gauss-001.zip --output artifacts/validation/frozen-source-gauss-continuation-next.json
```

Next bounded validation work: evaluate the existing 54/216 source solutions on
ONE COMMON physical output grid: the retained 864-source Gauss12/16/20 targets,
normals and weights. Reuse the 864-source results with exact archive/source/
physical identity checks. Repeat each source solution's own native calibration
first. Comparing source solutions at identical targets separates changes in
source resolution from changes in output nodes. The physical hull is identical;
verify its complete facet/surface identity, not just a common mesh label.

Keep output-integration and source-discretization uncertainty separate and charge
both before any source gate can pass; declare the complete policy before launch.
The observed order-20 rate projects about195 s for all new coarse-source field
work (an estimate, not guaranteed time); declare360 s total serial offline cap.
No new BEM solve or remeshing is needed. Deliver all36 source increments,
conditional quadrature gates and attribution, review and a fresh checkpoint,
then stop. This three-source diagnostic cannot qualify or replace the original
failed four-level study. Preflight the3456-source evaluation on this COMMON grid
separately: its order20 alone projects about1247 s offline, plus any required
quadrature controls. Do not multiply both source and target counts without a
need. Complete declared source levels, finite-frequency A/B, complex excitation/
phase and all required physical gates remain necessary. No Topwater runtime or
coupled-motion qualification follows.
