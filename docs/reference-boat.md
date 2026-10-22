# Illustrative 24 ft boat

`tests/reference_boat.py` defines a simple open hard-chine V hull for the new
boat application. It is an authored example, not a measured center-console
design or a reconstruction of a manufacturer's boat. The canoe remains a
separate numerical and stress-test fixture.

The hull is 24 ft long, 9 ft wide and 1 m deep, with 20 degree bottom deadrise.
The chine spans 80% of the maximum beam and the topsides flare outward. The
aft three quarters have a constant section; the forward quarter tapers to a
raised stem. Its physical transom and sides close below an open sheer rim.
A marked artificial cap closes the mesh for geometry validation; it carries
neither mass nor pressure. Every physical facet is planar, so triangle
subdivision preserves the same hull exactly.

Round loading assumptions total **1250 kg**:

- **600 kg shell**, uniform mass per area over the physical hull and transom.
- **400 kg amidships point load**, representing console, adults and gear,
  at body x=0, y=0, z=0.8 m.
- **250 kg motor point load**, at the transom centerline, z=0.7 m.

Body x points forward, y to port and z upward from the keel; x=0 is amidships.
The shell area determines its uniform areal mass from the stated 600 kg budget.
No thickness, laminate density or hidden extra payload is inferred. Triangle
integration gives exact mass, CG and the full inertia tensor, including the
longitudinal/vertical product. Point masses have offset inertia but no
additional intrinsic inertia. This deliberately omits detailed loading,
fuel/slosh and equipment geometry.

`reference_body(refinement=0)` returns the complete `RigidBody`.
`equilibrium(refinement=0)` independently solves flat-water heave and pitch
using Topwater's actual clipped triangles, with water density 1025 kg/m³ and
gravity 9.81 m/s². The returned pose has zero velocities. Tests check force and
moment balance, freeboard, local hydrostatic stability, independent quadratic
mass quadrature and invariance under geometry-preserving refinement. Computed
CG, draft, trim and inertia are derived outputs, not additional precise inputs
or claims of real-boat validation.
