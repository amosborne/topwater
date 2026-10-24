# Illustrative 24 ft boat

`tests/reference_boat.py` defines a simple open hard-chine V hull for the new
boat application. It is an authored example, not a measured center-console
design or a reconstruction of a manufacturer's boat. The 24 ft boat is the sole application reference; the former canoe case is retired.

The hull is 24 ft long, 9 ft wide and 1 m deep, with 20 degree bottom deadrise.
The chine spans 80% of the maximum beam and the topsides flare outward. The
aft three quarters have a constant section; the forward quarter tapers to a
raised stem. Its physical transom and sides close below an open sheer rim.
A marked artificial cap closes the mesh for geometry validation; it carries
neither mass nor pressure. Every physical facet is planar, so triangle
subdivision preserves the same hull exactly.

Round loading assumptions total **2500 kg**:

- **1600 kg hull/structure**, uniform mass per area over the physical hull and transom.
- **600 kg amidships point load**, representing console, adults, gear and fuel,
  at body x=0, y=0, z=0.8 m.
- **300 kg motor point load**, at the transom centerline, z=0.7 m.

Body x points forward, y to port and z upward from the keel; x=0 is amidships.
The shell area determines its uniform areal mass from the stated 1600 kg budget.
No thickness, laminate density or hidden extra payload is inferred. Triangle
integration gives exact mass, CG and the full inertia tensor, including the
longitudinal/vertical product. Point masses have offset inertia but no
additional intrinsic inertia. This deliberately omits detailed loading,
fuel slosh and equipment geometry. The distributed mass represents the boat's
structure with one thin-shell distribution; it is not a literal outer-laminate
thickness or a separate deck model.

The total is informed by marine-gyro's collected specification catalog and its
standard loading rules (full fuel, installed engines, two adults and gear).
The 33 center-console entries from 23–25 ft span about 1700–3900 kg with a median
near 2700 kg. They are not a sales-weighted sample and some engine/liquid bases
are interpreted. The round 2500 kg reference and its three-way allocation are
authored simplifications, not fitted measurements of a particular boat.
The source snapshots used were `data/boat_catalog.json` (SHA256
`11349cde737217c64bc8e9390fd64f1f56cba7bf7129637fe93f1fcbf6e862fc`)
and `data/loading.json` (SHA256
`0f7308902879938107b4ffbc0ed63b91254de8cd039debe1a729418a5232f7cd`).

`reference_body(refinement=0)` returns the complete `RigidBody`.
`equilibrium(refinement=0)` independently solves flat-water heave and pitch
using Topwater's actual clipped triangles, with water density 1025 kg/m³ and
gravity 9.81 m/s². The returned pose has zero velocities. Tests check force and
moment balance, freeboard, local hydrostatic stability, independent quadratic
mass quadrature and invariance under geometry-preserving refinement. Computed
CG, draft, trim and inertia are derived outputs, not additional precise inputs
or claims of real-boat validation.
