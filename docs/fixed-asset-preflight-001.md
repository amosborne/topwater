# Fixed-seed assets verified; fresh A∞ source pilot required

Corrected preflight002 at source `5c55945` completed read-only in3.364s below
60s. Exact command chunkfb7674 is terminal exit0. Current complete hull, calm
InitialState including zero velocities, CG, orientation, full dry inertia,
restoring and wetted surface match the archived canonical problem BITWISE.
Serialized physical input bytes also match. All four mesh copies agree across
freeze, raw study and periodic source recovery; physical conventions and named
Green profiles match their actual recorded settings.

The [prospective protocol](fixed-asset-preflight-protocol.md) and
[controls](../validation/fixed-asset-preflight-controls.json) were frozen before
execution. Attempt001 at `1035672` stopped on a reader error: the named finite
Green profile was incorrectly treated as a dictionary. It is preserved with its
actual producer/input graph, complete=false/passed=false. The corrected reader
compares recovery settings to the original recorded finite dictionary and checks
the infinite0/0 profile. Neither attempt imported/evaluated Capytaine or solved
sources. Independent review also required an explicit full calm-state guard,
which was added and tested before execution.

| Mesh panels | Seed subdivision | Max aspect indicator | Min angle (degrees) |
|---:|---:|---:|---:|
|1232|2|315.12536|0.181820|
|2772|3|315.12536|0.181820|
|4928|4|315.12536|0.181820|
|7700|5|315.12536|0.181820|

All seven physical facets, containment/moments/conformity,154 half-seed affine
patches, reflection and first-vertex maximum-radius convention pass unchanged
geometry controls. Max lattice errors stay below2.09e-9 versus1e-7 control.
Archived actual Capytaine radius proofs are zero; the new audit independently
checks the geometric convention without claiming a new native radius call.
Full per-facet/global shape/edge/area statistics are retained. The aspect range
is far smaller than the rejected paired-family value50600.99, but geometry quality
is not hydrodynamic convergence, and slivers remain unproven as the error cause.

Source inventory prevents two important false reuse claims:

- Existing periodic archive has six radiation modes at each of three frequencies
  and one paired diffraction heading per frequency at each mesh:84 source solves.
  Six of nine original frequency-heading diffraction cases per mesh are missing,
  requiring24 additional finite-frequency exports for complete original coverage.
- raw_Ainf contains4×6×6 coefficients, NOT A∞ densities. Fresh six-mode A∞ export
  on each of four meshes requires24 source solves. Finite densities use their
  recorded finite Green function/unit-displacement normalization; they cannot be
  relabeled as infinite-frequency acceleration-potential densities.

Common physical output grids from the same308 seed faces pass strict owned-
interior/area/normal checks:44352/78848/123200 targets for Gauss12/16/20. All three
orders across all four source meshes would cost4,098,124,800 source-target pairs;
the finest alone1,897,280,000. This is a work count, not authorization or timing
proof. Original dense source-export workspace estimates are0.399/1.633/4.891/
11.703GB under the original12GiB cap, requiring a fresh host check before solves.

Peak working set222,322,688B, peak committed1,625,776,128B and accounted disk
56,442,082B are below the declared2GiB caps. All287 local tests pass; exact test
session81145 is terminal exit0, JUnit included. Independent postrun review
reconstructed physical identities, geometry/quality, inventories, common grids,
work/resource counts and every archive/producer binding. Fresh restoration of
both attempt graphs and exact inputs passed. No native numerical job is live.

complete=true/passed=false/assets_usable=true means ASSET ADMISSION ONLY. Full
reference/source qualification remain false; all original study003 coefficient
failures and incomplete nonlinear-restart status remain unchanged.
[Complete evidence](../validation/evidence/fixed-asset-preflight-001.zip),
[summary](../validation/evidence/fixed-asset-preflight-001-summary.json).
SHA256 `cd66598c22357706b8517065a886e79742e887500a63c388c559fc70e7fc924f`.
Restore with `python restore_inputs.py --frozen-bundle PATH --raw-bundle PATH
--sources-bundle PATH` using frozen-fixed-001.zip, reference-study-003.zip and
periodic-sources-003.zip. Both historical attempts and all raw inputs remain bound.

Reproduce with a fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.fixed_asset_preflight --frozen-bundle validation/evidence/frozen-fixed-001.zip --raw-bundle validation/evidence/reference-study-003.zip --sources-bundle validation/evidence/periodic-sources-003.zip --output artifacts/validation/fixed-asset-preflight-next.json
```

Next bounded goal: fresh SIX-MODE A∞ source export and pressure-output quadrature
PILOT on the COARSEST1232-panel fixed-seed mesh only. Prospectively declare600s
total/2GiB process and additional-disk caps, one native thread and100-target
blocks. Archive actual producer/physical input before launch. Use the original
recorded infinite Green settings, same surface/CG/modes and explicit symbolic
infinity normalization. Save all six normalized densities/potentials/BCs and
operator-produced traces/residuals before releasing pristine S/K and source-
solve caches for output evaluation. Native potentials/traces and all36 original
raw A∞ entries must reproduce at1e-10 scaled tolerance before comparison.

Evaluate actual frozen sources on the308-seed COMMON Gauss12/16/20 grids:
303,564,800 field pairs. Reuse the checked consumer buffer/checkpoint/verification
path with ACTUAL fields; no synthetic data. Declare unchanged all36 0.001relative
+1e-6 normalized quadrature increments/trend and separate Uq before launch. Native
reproduction or integration failure cannot pass. Deliver raw field/source/phase/
load/uncertainty evidence, independent review, checkpoint/tracker/handoff and
stop after the pilot. No all-level source or finite-wave promotion, no remeshing,
new formulation, runtime optimization or long paired-family continuation. The
other18 A∞ exports and24 missing finite diffraction exports remain future work;
all original36-entry/finite-wave/phase gates still require qualification.
