# Geometry tooling boundary

The workflow ingests one convex physical triangular hull with one uncapped
simple opening. `boat_geometry.physical_hull` infers its boundary and
creates a nonphysical internal cap for the existing `HullMesh` API. Saved boats
and results contain physical faces only. No implicit repair or welding occurs.
General `HullMesh` and its existing mathematical coverage are unchanged.

## Current implementation

NumPy arrays remain the exchange representation. SciPy supplies convex hulls,
spatial search, connected components and linear algebra. Gmsh 4.15.2 is already
a correlation dependency and is used by historical validation mesh generators.
Capytaine owns its external-reference mesh and hydrodynamic calculations.
No general mesh-processing library has been adopted as native runtime geometry.

## Executed checkpoint

The installed Gmsh constrained-facet probe passed on the ordinary and rigidly
rotated/scaled fixtures: the hole, collinear boundary vertices, segment chains,
orientation and area were preserved. The probe is part of the passing interface
suite (188passed/1optional skip). libigl is not installed and its gated probe
was explicitly skipped; no capability or adoption claim follows from that skip.

The current export helper uses bounded convex-halfspace triangle clipping,
preserves physical planes, and propagates cut vertices through neighboring
exterior faces before triangulation. Its flat/sloped/wavy contact, scale and
opening tests passed. This does not establish a general mesh Boolean kernel,
nonlinear fluid-state remeshing or arbitrary degeneracy support. Gmsh remains
an optional correlation dependency, not a newly required native runtime backend.

The new ingestion boundary checks edge and vertex manifoldness, one connected
disk, consistent outward physical normals, convex supporting planes, geometric
triangle intersections/T-junctions, and one simple rim. It preserves all
physical vertices/faces. Geometry predicates operate after translation and span
normalization; the contact tolerance is 1e-10 of maximum hull span. The internal
cap for a planar opening is a fan about a new rim-centroid vertex. A nonplanar
opening is supported when the missing patch comprises whole convex-envelope
facets; SciPy identifies the facets and the existing constrained ear routine
triangulates them while retaining collinear rim vertices. Physical faces and
their geometry remain unchanged. A valid nonplanar opening cutting a partially
occupied convex facet raises `UnsupportedClosureError`: that is an internal
closure-backend limitation, not a requirement that user inputs be planar.
No generic polygon-Boolean engine is added. Caps receive neither pressure nor
shell mass. The constrained single-opening contract is not a general self-
intersection repair or unrestricted CAD importer.

`workflow_data.Result` validates numeric storage, absolute poses, hull validity,
packed surface indices and nondegenerate triangles. Its surface topology may
vary and its sampling times may differ from body times. It does not certify a
physical free-surface solution, watertight contact, flooding support or accuracy.
Those checks belong to the surface producer and solve acceptance boundary.

## Replace operations, preserve physics

Candidates for off-the-shelf operations are robust mesh validity and ordered
boundary-loop extraction; hull/surface intersection and trimming with shared
contact topology; and constrained polygon triangulation preserving vertices,
holes and physical-face lineage. Keep NumPy inputs/outputs and explicit tags.
Do not export package object types or silently replace an invalid hull by its
convex envelope.

Keep pressure/mass/inertia integration, frame ownership, physical/artificial
panel selection, wave construction, fluid boundary conditions and state transfer
in Topwater. A geometry package does not choose the fluid domain or flooding
policy. Existing exact flat-water integration is a useful comparison oracle.

## Bounded executable feasibility probes

Source is in `tests/test_mesh_tooling.py`. These are decision gates, not a
package-adoption claim. Run only after coordinator
allocation of the sole numerical lease. Check the installed package identity and
API first; do not install a package merely to satisfy an unallocated probe.

**Gmsh:** the executable probe uses its Python curve-loop surface interface on one
planar facet containing an outer boundary and inner hole. Include collinear
boundary vertices and fixed contact segments. Generate a triangular surface,
then extract arrays and verify every supplied constrained edge is preserved,
there are no triangle centroids in the hole, total area matches the outer area
minus the hole, all face orientations agree, and there are no unused vertices
or duplicate faces. It repeats after rigid transform/scaling. Original point
and curve tags recover supplied boundary identity. Near-coincident constraints,
general polygon differences and dynamic hull/water trimming remain further
capability questions. Stop on the first unsupported operation or failed
geometric invariant; no boat simulation or broad mesh sweep.

**libigl:** the executable check reports whether the `libigl` Python distribution
is installed and skips without installing it. Execution requires separate
distribution-review approval (`TOPWATER_LIBIGL_PROBE=1`). It then checks the actual
boundary-loop and edge/vertex-manifold bindings and extracts one known boundary.
It reports a candidate self-intersection binding separately; mere presence does
not qualify robustness. After API/distribution review, extend the accepted API
probe to an open box, a second opening, a bow-tie vertex and a non-manifold edge,
plus crossing/coplanar triangles and permitted shared edges. If a binding is
absent, report it as unavailable rather than substitute an unreviewed remesher
or claim C++ capabilities are available in Python. No automatic repair.

The Gmsh probe reports complete call wall time and Python allocation peak;
native-package peak memory is explicitly unmeasured. Record preserved constraints/
identity, failure behavior and API availability.
They are feasibility checks, not dynamic-meshing or physical-runtime benchmarks.

## Acceptance before adoption

- Closed/capped input, extra openings, invalid winding, pinched vertices,
  non-manifold edges, geometric crossings and T-junctions fail explicitly.
- Tangent, edge/vertex and coplanar contact have defined behavior; ambiguous
  contacts are rejected rather than invented.
- A tilted convex hull under a nonflat graph has matching wet-hull/free-surface
  edges, no cracks, no interior water triangles and consistent orientations.
- Concave wet facets retain all collinear boundary constraints. Rigid transforms,
  face reordering and supported scale changes preserve classifications.
- Existing box/canonical calm cases retain displacement, wet area and integrated
  loads; artificial caps remain excluded. Variable surface connectivity remains
  supported by the result archive.

Package adoption requires a separate coordinator decision using these results
and the independent tooling review. No replacement or speed claim follows from
this document alone.

## Export allocation evidence (2026-10-25)

The exporter now materializes the completed cut-vertex registry once and reuses
it for triangulation. Shared ear clipping, predicates, tolerances and point
indices are unchanged; local projection changes were deferred. Independent
source review approved this bounded allocation change.

The coordinator's serial check completed successfully: 16 tests passed in
0.40 s, with frozen input identities verified. A separate whole-export comparison
loaded the actual old module from `a29600b` and the revised module on one 76x76
flat-grid/box fixture (5,776 nodes; 11,250 faces). Vertices and faces were bitwise
identical in both calls per version. Old whole-export times were 1.861/1.627 s;
registry-reuse times were 1.336/1.422 s. See the local
[report](../artifacts/workflow/export-allocation-attempt2/report.json),
[terminal evidence](../artifacts/workflow/export-allocation-attempt2/terminal.json)
and [test log](../artifacts/workflow/export-allocation-attempt2/output.log).

This is one geometry fixture, with two timings per version. It establishes no
physical accuracy, wet-area cost, or complete canonical-workflow speed claim.

## Constrained-contact triangulation correction (2026-10-25)

The second canonical attempt stopped at saved frame 277 (22.16 s). One
four-vertex export polygon contained an intermediate cut vertex on a constrained
contact edge. Projected ear clipping accepted a roundoff-scale final sliver
whose world-coordinate triangle failed the unchanged degeneracy validator. The
saved diagnostic found one bad triangle with area 3.04944854e-20 m2; no body,
wave, or hydrodynamic failure was inferred. See the
[one-frame diagnosis](../artifacts/workflow/surface-frame-diagnostic-attempt1/frame/report.json).

The export-only correction flips an unconstrained internal diagonal when both
replacement triangles pass the existing world-coordinate predicate, preserve
orientation and signed/unsigned area, and retain every vertex and directed
boundary edge. Otherwise it fails explicitly. A reconstructed four-vertex
exception path additionally verifies oriented parent coverage. Shared generic
ear clipping is unchanged. No threshold was relaxed, no positive-area triangle
was blindly dropped, and no hull/contact point was moved. The Capytaine owner
independently approved the source correction.

The coordinator's current interface suite passed 228 tests with one optional
skip in 14.98 s. A geometry-only verification of actual saved frame 277 and
cached neighbors 276/278 passed in 8.423 s. Reconstructed frame-277 inputs were
bitwise identical to the saved failure geometry; corrected output vertices and
directed boundary chains were unchanged, with exactly two triangles over four
vertices changing their internal diagonal. Area/contact constraints were
preserved and all three frames passed topology checks. See the
[three-frame report](../artifacts/workflow/surface-frame-verification-attempt1/frames/report.json).

The all-frame cached geometry sweep helper is frozen and awaiting execution.
Three accepted frames do not establish complete-case export, video completion,
physical accuracy, or native-solver readiness.
