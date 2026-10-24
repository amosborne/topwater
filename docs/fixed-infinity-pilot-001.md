# Fresh A∞ sources and common-grid integration pass the coarsest pilot

At source `cde8f2c`, all six fresh infinite-frequency radiation cases on the
1232-panel fixed-seed mesh pass native reproduction, and all36 pressure-output
integration entries pass the unchanged12→16/16→20 stabilization and trend gates.
Maximum tolerance ratios are0.554736/0.139239. Exact session59843 is terminal
exit0, elapsed222.557s below600s. complete=true, passed=false, pilot_passed=true,
source_export_verified=true. One source mesh does NOT establish source-mesh or
finite-wave convergence; all original study003 failures remain failed.

The [prospective protocol](fixed-infinity-pilot-protocol.md) and
[controls](../validation/fixed-infinity-pilot-controls.json) were reviewed and
committed before native work. Actual sources/potentials/BCs were normalized by
(-i*encounter_omega) before symbolic infinity cancellation. All six real/finite
bases, body rigid normals, pristine Sσ=χ and Kσ=q, native pressure integration
and all36 original raw A∞ entries reproduce at1e-10 scaled tolerance. Native
solver A equals the original matrix EXACTLY. Largest source/trace tolerance ratio
is0.000112691 OF that1e-10 tolerance, not a physical relative error.

Pristine S/K were copied before LU ownership changes; hashes/shapes, normalized
bases, operator-produced traces and residual vectors are retained. Dense matrices
and solve caches were released before field evaluation. Export took0.618s.
The same source surface/worldCG/full dry inertia, depth/axes, density1025/gravity
9.81, Green settings and pinned dependency/native-library identities were used.
No fitting, mirrored field values, zeroed coupling or matrix averaging occurred.

Actual fields on COMMON308-seed grids44352/78848/123200 points took39.565/71.757/
104.234s. Both reflected halves were independently evaluated through the checked
caller-owned buffer/streamed checkpoint/verification/package pipeline. Every
point/owner/normal/positive weight and all complex fields are retained. The Uq
envelope and all36 signed/absolute output changes, increment ratios/failure masks
and trend diagnostics remain available; no extrapolation or continuum bound.

| Diagonal | Original native pressure output | Gauss20 output |
|---|---:|---:|
|Surge (kg)|91.115183|88.688485|
|Sway (kg)|577.612392|568.894688|
|Heave (kg)|8153.427509|8111.035902|
|Roll (kg m²)|1004.069656|985.037507|
|Pitch (kg m²)|16574.584853|16469.405633|
|Yaw (kg m²)|1329.854125|1306.984294|

These measured integration changes are not source-discretization error estimates.
Source uncertainty still needs the other three meshes and separate charged Uq.
The conditional gate remains0.001relative+1e-6 normalized, both finest increments
and nonincreasing trend. No old native result was relabeled as converged.

Peak working set289,247,232B, peak committed1,732,276,224B and accounted disk
112,069,679B pass the declared2GiB caps. All288 local tests pass; exact test
session73676 is terminal exit0, JUnit included. Independent postrun review
reconstructed native/source/target identities, loads, traces, Uq/masks, file
roundtrips and resource arithmetic. Fresh complete graph/input restoration
passed. No native job is live. Offline costs are not Topwater simulation speed.

[Complete actual-source/field producer data](../validation/evidence/fixed-infinity-pilot-001.zip),
[summary](../validation/evidence/fixed-infinity-pilot-001-summary.json).
SHA256 `e4d43719636db7ac67731b5f98d418b247371b585f2bd89069b3405eaaee2e79`.
Immutable payload contains the NPZ and actual code/physical graph; the admitted
asset archive is embedded. Required external inputs: frozen-fixed-001.zip,
reference-study-003.zip, periodic-sources-003.zip. Restore after extraction with
`python restore_inputs.py --frozen-bundle PATH --raw-bundle PATH --sources-bundle
PATH`. The payload-manifest snapshot precedes final package metrics; both it and
final metadata are SHA-bound, and the final NPZ stays immutable.

Reproduce with a fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.fixed_infinity_pilot --preflight artifacts/validation/fixed-asset-preflight-002.json --preflight-bundle validation/evidence/fixed-asset-preflight-001.zip --output artifacts/validation/fixed-infinity-pilot-next.json
```

Next bounded goal: SOURCE-EXPORT ONLY for the remaining2772/4928/7700 meshes,
six modes each (18 solves). Declare600s total, prospective8GiB NATIVE process and
2GiB disk caps before launch; this separate dense-export workload is different
from the2GiB cold field consumer. FIRST inspect actual installed allocation/LU/
cache ownership, model every dense copy/workspace and host margin, and reject
before launch if the explicit model cannot fit. At7700, EACH complex matrix is
948,640,000 bytes. Match warmup/solve cache keys; row-stream canonical operator
hashes instead of allocating giant tobytes copies. Preserve pristine checks and
all36 native/raw-A reproduction at unchanged1e-10 per mesh. Retain normalized
σ/χ/BC, operator-produced traces/residuals/hashes, then release each dense state.

Archive actual producer/input/model/host identities, run serially with exact
terminal-process discipline, monitor actual memory, and deliver reviewed complete
sources or concrete rejection with checkpoint/tracker/handoff. No Green field
grids, finite-wave solves, remesh/formulation/runtime/optimization work. Do not
promote full/source convergence. Later common-grid field evaluations must use
separately reviewed COLD processes with2GiB caps and fresh cost/quadrature gates;
native peak counters must not be confused with field-consumer peaks. No automatic
field job follows this export. The24missing finite diffraction cases and every
full36-entry/phase/reference gate remain open.
