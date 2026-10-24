# Fixed-topology reference mesh preflight

The overarching reference-qualification goal remains active and incomplete.
Previous checkpoint ec7b472 made progress: two complete failed studies are
preserved, vertex-order sensitivity is measured and its reviewed numerical
control correction is implemented. No hydrodynamic job remains live.

## Bounded scope and stopping checkpoint

Design and geometry-preflight a fixed-topology graded reference mesh family on
the exact frozen physical surface. No hydrodynamic solves, Topwater runtime
changes, convergence tolerance changes, new solver formulations, optimization
or timing work are part of this checkpoint. Deliver the archived seed, patch
definitions, shared-curve distributions, four complete meshes, per-facet shape
statistics, surface/input equality proofs, resource preflight and independent
review. Stop at a reviewed geometry decision; commit/push, update issue #30 and
write the next bounded convergence-study handoff. Geometry success does not
qualify the hydrodynamic reference.

## Reviewed construction

Independent reviewer fixed_topology_review rejected adoption of a centroid-quad
fan: opposite-edge constraints would tie short chine-edge subdivisions to long
edge counts and preserve poor shapes over the strip. No numerical job ran for
that candidate. Use the reviewer's simpler alternative instead:

1. Generate ONE graded seed with the existing pinned Gmsh 4.15.2 mesher and
   verify the complete seven-facet physical surface. Initial declared controls
   are bulk size 1.2 m, edge size 0.3 m and distance growth 0.3. Preserve all
   narrow strips. Archive the actual seed before family generation.
2. Treat each half-surface seed triangle as a fixed planar transfinite patch.
   Reuse shared seed edges, each with uniform integer subdivisions N. Use Gmsh's
   structured triangular algorithm, with smoothing/optimization disabled. The
   grading map is the seed's frozen piecewise-affine map, sampled more finely;
   it does not change across levels. No independent unstructured remeshing.
3. Initial declared N are 2, 3, 4 and 5. Full panel counts must equal seed_count
   times N squared. The target is approximately the previous final 7,600 panels,
   but actual counts, not estimates or a silently dropped level, control preflight.
4. Reflect the verified half-surface and apply farthest-vertex-first ordering.

Gmsh supports explicit transfinite curve counts and the structured triangular
algorithm ([Gmsh manual](https://gmsh.info/doc/texinfo/gmsh.html)). Every actual
patch must match the declared affine barycentric lattice and fixed triangle
connectivity; a Gmsh fallback to unstructured meshing cannot pass unnoticed.

## Acceptance declared before geometry jobs

- Existing geometry controls: coordinate containment within 1e-11 m; per-facet
  area, oriented area and first/second moments within 1e-10 relative + 1e-12 SI
  absolute; outward winding, conforming opposing edges, waterline-only open
  boundaries and no added closure/lid. Reflection must preserve the actual
  surface. All original physical arrays and input JSON must match the reset.
- Each patch has exactly N squared positive-area children, affine barycentric
  coordinates within 1e-7 dimensionless of the integer lattice, and the exact
  declared lattice connectivity. Child areas agree with patch_area/N squared
  within 1e-8 relative + 1e-14 m squared absolute. These are geometry-roundoff
  checks, not hydrodynamic accuracy tolerances.
- Record the seed, child and per-physical-facet edge/angle/aspect distributions;
  review avoidable slivers separately from the genuine acute physical corners.
  Do not infer hydrodynamic suitability solely from a surface proof.
- Every declared level must fit the unchanged canonical workspace estimator
  and explicit 12 GiB estimate cap. Keep the existing 1800-second elapsed limit.
  This estimator is not an OS memory limit; no resource relaxation is implied.
- Preserve the unchanged full hydrodynamic case product and acceptance policy
  for the next study. This geometry-only report is explicitly unqualified and
  cannot satisfy a reference or physical-promotion gate.

If the initial seed exceeds resource limits or has an unresolved quality defect,
preserve the failed candidate, diagnose it and explicitly declare/review the next
geometry correction before another job. Do not tune meshes against hydrodynamic
results or weaken a convergence gate.

## Candidate 001 and declared correction 002

Candidate 001 completed seed generation, then correctly failed resource preflight:
382 seed panels imply 1528/3438/6112/9550 panels. No family or hydrodynamic solve
ran. The complete seed/producers and bound freeze are in
[candidate 001 evidence](../validation/evidence/transfinite-preflight-001.zip).
Broad facets had minimum angles above 38 degrees, but the chine strips had only
seven seed panels per side, aspect ratios 79–315 across the strip. This is a
seed-quality issue to correct before freezing a family, not a tolerance failure.

Reviewer fixed_topology_review approved another geometry-only candidate before
launch: bulk 1.5 m / edge 0.5 m / growth 0.3 coarsens broad facets, while each planar facet
receives a local maximum size 20 times its minimum support width. Width is
computed in the facet plane normal to each polygon edge, retaining the minimum;
it is not a coordinate bounding box. Gmsh Restrict fields include facet boundaries
so adjoining patches share the finer feature-edge distribution. Actual widths
and caps must be recorded. N remains 2/3/4/5; resource and accuracy gates are
unchanged. Candidate 001 remains preserved and unqualified.

The new optional seed control changes a validator source, so regenerate a fresh
canonical freeze before this candidate. Verify all eight original physical arrays
and input JSON against the reset. Candidate 002 must state canonical-contract
SHA256 separately from its geometry protocol SHA256 and embed complete bound
and reset freezes, resolving the review's provenance findings.

## Candidate 002 geometry decision

Candidate 002 passed at source checkpoint `46d35f6`. Its complete, self-contained
[evidence bundle](../validation/evidence/transfinite-preflight-002.zip) embeds
both bound/reset freeze JSON, NPZ and producer archives. The
[summary](../validation/evidence/transfinite-preflight-002-summary.json) records
the bundle SHA256, canonical/protocol identities, complete physical equality,
seed and per-level/per-facet statistics.

The 308-panel seed produces 1,232 / 2,772 / 4,928 / 7,700 panels. All actual Gmsh patch
connectivity, affine lattice, winding/Jacobian, child area, shared-curve and full
surface checks passed. Maximum lattice error is 2.09e-9 versus 1e-7 allowed;
maximum moment difference is 1.28e-13 SI. Finest workspace estimate is
11,703,102,464 bytes, below the unchanged 12 GiB cap of 12,884,901,888 bytes.
Original physical arrays/dtypes and input JSON match reset exactly.

Thin-facet resolution improved from 7 to 31 seed panels per side; median aspect
fell from 157.6 to 39.4 (area-weighted median 27.2). The genuine acute apex still
has about 0.18 degree angle and 315 maximum aspect. Main-facet maximum aspect is
4.06 and stern area-weighted 99th percentile 4.27; the finer shared curves did
not introduce a broad neighboring sliver population. Globally, the area
fraction with aspect above 50 is 3.42e-5. These are shape observations, not a
physical accuracy argument.

Independent reviewer fixed_topology_review recomputed every patch and full
surface proof and verified archived/current producer hashes, explicit canonical
contract hash, both embedded freeze identities and physical equality. No
blocking geometry/provenance finding remains. Approval is limited to geometry
readiness for a full unchanged convergence study, not reference qualification.

## Selection for the next bounded study

The canonical contract now explicitly selects fixed-seed-transfinite, with the
same seed controls and N=2/3/4/5. The canonical freeze stores the seed and every
actual mesh, including farthest-first ordering, geometry proofs and shape
statistics. New producer identity includes reference_transfinite.py. Missing
selected meshes cannot fall back to midpoint subdivision, and stale owners
cannot pass. Missing seeds and incorrect seed_count times N squared mappings
also fail. The focused validation suite passes 51 tests. The independent reviewer
found no remaining blocking integration defect before the fresh canonical freeze.

Create a fresh canonical freeze after this integration, verify its original
physical arrays/input JSON and each mesh against candidate 002, and archive it
before the geometry checkpoint closes. Do not run hydrodynamics in this
checkpoint. The next bounded study then runs all 132 calls and unchanged
entrywise convergence/trend/uncertainty gates using those exact frozen meshes,
after a fresh resource-headroom check. No Richardson estimate or new reference
qualification follows from this geometry decision.

## Completed canonical adoption proof

The fresh canonical freeze passed at source checkpoint
96925446db194ed4f46e45ae74a3a62cf2d27525. Its seed and all four actual meshes
are bit-for-bit equal to candidate 002. All recorded current producers match
their archives. Original physical arrays/dtypes and input JSON remain exactly
equal to reset, and every non-reference-numerical contract section is unchanged.
Capytaine's actual radii equal maximum vertex-to-centroid distances on every
panel (zero difference in this verification).

[Complete canonical freeze/proof bundle](../validation/evidence/frozen-fixed-001.zip)
includes the freeze, producers, input JSON, proof/verifier and the complete
geometry candidate bundle. [Proof summary and SHA256](../validation/evidence/frozen-fixed-001-summary.json).
The archive README explains extraction and proof reproduction. No hydrodynamic
solve ran. Current reference qualification remains unproven; this completes
only the declared geometry checkpoint and prepares study 003.

Final independent adoption audit by fixed_topology_review verified actual
dtype/shape/bytes for physical arrays and seed/mesh comparisons, current
producer identities, physical identity and proof hashes. No blocking adoption
finding remains. Study 003 still requires its fresh OS-memory check and every
unchanged hydrodynamic convergence gate.
