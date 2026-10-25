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
