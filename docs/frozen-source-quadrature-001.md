# Frozen-source pressure integration: composite centroid rule rejected

The diagnostic at source `a2a8100` completed all four integration levels in
60.216 s, with zero new BEM solves. Exact session 82940 is terminal exit 1 by
design: complete=true, passed=false. This is external reference work, not
Topwater simulation time. The [prospective protocol](frozen-source-quadrature-protocol.md)
and [controls](../validation/frozen-source-quadrature-controls.json) were frozen
before launch.

The same 864-panel Capytaine sources, physical facets, world CG/full dry inertia,
density 1025 and infinite-frequency Green settings were used throughout. Native
centroid evaluation reproduced all six potential traces and all 36 original loads
within 1e-10 scaled tolerance. Refinements used 4/16/64 equal-area children per
source triangle: 3456/13824/55296 strict owned-interior trace points, plus 864 native
points. Every target/owner/weight/potential is retained. No interior-volume
continuation, new source solve or fitted pressure was introduced.

Pressure integration is a material contributor to the measured mesh increment.
The final terms below obey exact telescoping; units are kg for translations and
kg m^2 for rotations. Neither term is a certified continuum error.

| Diagonal | Measured finer-native minus coarse-native | Measured output-quadrature change | Unresolved remainder |
|---|---:|---:|---:|
| Surge |-2.180|-2.274|+0.094|
| Sway |-17.092|-22.725|+5.633|
| Heave |-311.678|-182.989|-128.689|
| Roll |-31.146|-57.642|+26.496|
| Pitch |-429.542|-636.547|+207.005|
| Yaw |-29.505|-78.373|+48.868|

Several terms oppose and partly cancel. The remainder includes changed source/
collocation solution AND finer-native output quadrature; it cannot be labeled
pure source error. All 36 signed matrices, absolute entry magnitudes and opposing-
component masks are preserved in the primary and derived attribution artifacts.

The composite-centroid rule remains under-resolved: the original 0.001 relative
+1e-6 normalized stabilization diagnostic fails 18 previous-increment entries and
17 latest-increment entries. Maximum tolerance ratios are 596.21/282.18; the trend
passes. No threshold was relaxed, no native failed gate was relabeled passed,
and source/reference qualification remains false.

[Complete evidence](../validation/evidence/frozen-source-quadrature-001.zip),
[summary](../validation/evidence/frozen-source-quadrature-001-summary.json).
SHA256: `c3e93235371cf72ef296dab1245afacb9a7e5e9ec47ca9621643717e0a053ada`.
Run `python restore_shared_inputs.py` after extraction to restore both complete
graphs and shared parent/data payloads with SHA verification. Fresh extraction
was verified. All 268 local tests passed, including ownership, weights, moments
and telescoping; JUnit is included. Independent review verified all 73440 points,
four pressure matrices, stabilization masks and actual provenance. No numerical
job remains live; tests session 45615 is terminal exit 0.

Reproduce with a fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.frozen_source_quadrature --study artifacts/validation/paired-body-study-002.json --bundle validation/evidence/paired-body-001.zip --output artifacts/validation/frozen-source-quadrature-next.json
```

Next bounded correction: a prospectively declared Duffy–Gauss pressure rule on
the same frozen 864 sources, using mapped [0,1] Gauss nodes, barycentric coordinates
((1-u)(1-v),u,(1-u)v) and relative weights 2(1-u)wu*wv. Proposed orders 8/12/16 plus
native reproduction need 55296/124416/221184 targets and a separately declared
600 s offline ceiling. Test polynomial moments and strict ownership. Construct
mirror-consistent TARGET geometry on the verified paired mesh, then evaluate
the full source field independently on BOTH halves; never mirror potentials or
average/zero matrices. Retain unchanged all 36 stabilization gates. No spectral
rate, Richardson extrapolation, pure-source remainder or full-reference claim.
