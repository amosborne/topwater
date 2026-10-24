# Frozen-source Duffy–Gauss pressure integration: not yet stabilized

At source `d7af5e0`, native centroid plus Gauss orders 8/12/16 completed
401,760 target evaluations in 324.115 s, below the declared 600 s offline cap.
Exact session 64310 is terminal exit 1: complete=true, passed=false. No new BEM
solve ran. This is external reference work, not Topwater simulation speed.
The [prospective protocol](frozen-source-gauss-protocol.md) and
[controls](../validation/frozen-source-gauss-controls.json) were frozen before
launch and independently reviewed.

All 864 normalized Capytaine sources, physical facets, world CG/full dry inertia,
density and infinite-frequency Green settings stayed fixed. Native evaluation
reproduced all six potentials and all 36 loads at the unchanged 1e-10 scaled
tolerance. Higher-order targets stayed strictly inside their owned facets;
positive weights preserve each actual face area. Only target geometry was
reflected; the complete source field was evaluated independently on both halves.
No field reflection, matrix averaging, forced zero coupling or new source solve
was used.

The unchanged 0.001 relative +1e-6 normalized integration diagnostic still fails.
The 8→12 increment fails only entry [1,5], sway force from yaw acceleration,
at 1.562861 times its tolerance. All 12→16 entries pass, maximum ratio 0.246215;
the trend passes. The previous centroid rule and native source-mesh studies
remain failed. One passing increment does not establish stabilization.

| Integration | Sway–yaw entry (kg m) | Pitch diagonal (kg m²) |
|---|---:|---:|
| Native centroid | -13.64577667 | 16947.629029 |
| Gauss 8 | -0.79123661 | 16275.715292 |
| Gauss 12 | -0.78195545 | 16275.435746 |
| Gauss 16 | -0.78049329 | 16275.390212 |

At order 16, measured pitch output change is -672.238818 kg m², against the
native finer-minus-coarse change -429.542381 kg m². Their unresolved remainder
is +242.696436 kg m². These terms exactly telescope; the remainder includes
changed source/collocation solution AND finer-native output quadrature. It is
not pure source error, and neither term is a continuum uncertainty bound.
All 36 complex matrices, targets/owners/weights/six potentials, signed and
absolute terms, opposition masks and gate ratios are retained.

[Complete evidence](../validation/evidence/frozen-source-gauss-001.zip),
[summary](../validation/evidence/frozen-source-gauss-001-summary.json).
SHA256: `28a35794d8c20d35e59f955e0b4ced8fec3cdebd4c6f8e01beffd1cb205fdbb0`.
After extraction run `python restore_shared_inputs.py`; fresh restoration of
both complete producer/input graphs was verified. Independent postrun review
reconstructed every target, weight, pressure matrix, telescope and failure mask.
All 272 local tests passed; JUnit is included. Test session 85636 is terminal
exit 0. No native numerical job remains live. Source/reference qualification,
coupled nonlinear motion and practical runtime remain unqualified.

Reproduce with a fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.frozen_source_quadrature --controls validation/frozen-source-gauss-controls.json --study artifacts/validation/paired-body-study-002.json --bundle validation/evidence/paired-body-001.zip --output artifacts/validation/frozen-source-gauss-next.json
```

Next bounded goal: prospectively append one order-20 evaluation on the SAME
frozen sources (345,600 targets), after native calibration; reuse orders 12/16
only with exact producer/input/array identity checks. The measured order-16 rate
projects about 279 s for order 20, offline only. Declare a 400 s total ceiling
including preparation before launch. Retain the original failed study; judge
12→16 and 16→20 with the same all-36 criteria and trend against order 20.
No inferred spectral rate or extrapolation. Deliver immutable evidence,
independent review and a checkpoint, then stop. Even a passing conditional
integration gate leaves source-mesh and finite-wave qualification outstanding;
no remeshing, major optimization or full motion run belongs to this step.
