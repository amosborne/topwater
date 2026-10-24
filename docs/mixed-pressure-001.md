# Closed-boundary pressure component and discrepancy attribution

This checkpoint adds conforming fluid geometry and an owned acceleration-potential
pressure/body calculation on the existing P0 Rankine solver. It does **not** add
a complete nonlinear RHS or qualify a fluid/body simulation. Source checkpoints:
pressure probe `6213d96`, discrepancy attribution `902c2be`.

The [pressure protocol](mixed-pressure-protocol.md) and
[attribution protocol](pressure-discrepancy-protocol.md) were committed before
their separate serial runs. Both artifacts completed with passed=false and exit 1
intentionally. The calm pressure probe uses empty waves, zero velocity/Hessian
and the exact canonical equilibrium. The saved mild state is geometry-only;
missing nonlinear derivatives are never replaced with zero for a pressure call.

Fluid geometry retains every contact vertex and actual wet facet, including
concave contours, then closes the free surface with a wall and bottom. Meshes
are connected, closed, consistently oriented manifolds. Artificial hull closure
faces carry no pressure. Shallow-bottom and opening contact are rejected.
PressureSystem owns geometry, body/pose, boundary assignment and factors together.
Seven mixed chi solves supply known pressure and six acceleration columns; dry
world inertia, gyro and weight enter once. Hydrostatic moment is integrated exactly.

The calm mesh has 762 panels: 588 surface/48 body/84 wall/42 bottom. The actual saved
mild surface closes into 2,616: 2,030/151/290/145. Calm acceleration is <=3.12e-15 SI;
combined-chi body pressure reconstructs the saved load. These checks validate
algebra and geometry, not fluid response accuracy. The six auxiliary chi traces
were not saved: this bundle cannot independently reconstruct all 36 raw A entries
from them. Preserve all seven traces in the next numerical producer.

The measured calm component costs were warm-up 0.279 s once, geometry 0.100 s,
assembly/factor 0.073 s, pressure 0.0062 s. Its provisional midpoint dt=.01 projection
is 36.17 wall s/s, excluding mandatory derivative recovery and all remaining RHS
work. This is a component estimate, **not** a practical-runtime pass. No mild
pressure cost or one-second simulation was measured.

The all-six inertia discrepancies are large. A follow-up isolates the same 48
body panels in Topwater's exact flat halfspace and Capytaine at infinite frequency.
Farthest-first triangle ordering, all seven physical facets, shared edges and area
moments were verified (maximum moment difference 7.11e-15). Six Capytaine calls
completed in the bounded diagnostic. All 144 entry comparisons reproduce the four
failed diagnostic outcomes under the existing 1%+1e-5 normalized candidate policy.
These are diagnostics; the finest reference itself remains unqualified.

| Diagonal | Closed P0 | P0 halfspace, same body | Capytaine, same body | Capytaine finest |
|---|---:|---:|---:|---:|
| Surge, kg |112.92|120.03|114.18|87.61|
| Sway, kg |749.14|841.33|667.56|570.06|
| Heave, kg |9190.36|11370.27|10094.21|8007.12|
| Roll, kg m^2 |1077.16|1256.94|1006.36|994.68|
| Pitch, kg m^2 |7802.21|8242.25|7128.24|16296.33|
| Yaw, kg m^2 |838.14|930.56|811.25|1312.53|

There are multiple unresolved contributions: body resolution sensitivity,
coarse direct/indirect discretization differences and explicit surface/domain
representation. Do not assign their entire difference to one owner or treat
same-body correlation as convergence. Preserve all 36 entries and the fine
reference's separate uncertainty envelope, retained in the bundle.

[Evidence bundle](../validation/evidence/mixed-pressure-001.zip),
[summary and shared-input bindings](../validation/evidence/mixed-pressure-001-summary.json).
SHA256: `ec77f90b85b6a101efa81b65cf714100ab95b9c4f7ca44ce4b9f20ce66914236`.
The two complete graphs share one identical 37 MB parent payload. After extracting,
run `python restore_shared_inputs.py`; both original paths are restored with SHA
verification. Fresh extraction, both owned artifacts and restored payload hashes
were checked. All 255 local tests passed, with test source and JUnit evidence saved.

Independent review verified signs, known normal derivative, physical/provenance
identity, pressure-load reconstruction and all 144 comparison entries. No blocker
remains for preserving this explicitly unqualified component investigation.
Both numerical commands were terminal exit 1; tests session 32046 terminal exit 0.
No numerical jobs remain live. No reference, wave, motion or speed gate is promoted.

Reproduce using fresh output stems and the locked environment:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.mixed_pressure --fields artifacts/validation/periodic-fields-001.json --sources artifacts/validation/periodic-sources-003.json --bundle validation/evidence/periodic-fields-001.zip --output artifacts/validation/mixed-pressure-next.json
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.pressure_discrepancy --probe artifacts/validation/mixed-pressure-next.json --output artifacts/validation/pressure-discrepancy-next.json
```

Next bounded validation: paired body-only refinement in Topwater halfspace and
Capytaine on the unchanged physical facets. Preflight controlled nested meshes
and actual counts/resources, then declare the complete study and all 36-entry
criteria before launch. Keep surface/domain fixed out of this comparison; no
Richardson assumption, same-mesh qualification or tolerance relaxation. Only
after attributing/correcting these errors return to explicit-surface/domain
controls and boundary-derivative recovery. Known-incident field representation,
initial compatibility, nonlinear/time controls and full runtime remain open.
