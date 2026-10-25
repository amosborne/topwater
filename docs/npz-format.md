# Canonical NPZ data contract

NPZ files contain only named numeric arrays. Load with `allow_pickle=False`.
Required members are exact: missing, extra, duplicate, object-valued or malformed
members are fatal. Units and conventions are fixed here; no metadata object,
solver marker, duration, version field or provenance report is embedded.

World axes are east/north/up; body axes are bow/port/up, both right handed.
All quantities use metres, kilograms, seconds and radians. Float arrays must be
finite. Connectivity and offsets are integer arrays.

## boat.npz

| Member | Shape | Meaning |
|---|---|---|
| vertices_m | N x3 | Complete physical hull, body coordinates |
| faces | F x3 | Outward physical triangles, zero-based indices |
| mass_kg | scalar | Total dry mass including the defined loading |
| cg_m | 3 | CG in the same body coordinates |
| inertia_cg_kg_m2 | 3 x3 | Full symmetric inertia tensor about CG, body axes |

Supply an uncapped convex shell with one opening. Its rim is inferred from
single-incidence edges; an internal closure never carries water pressure.
See [mesh validation and tooling](mesh-tooling.md) for current support limits.

## sea.npz

Four equal-length vectors: `amplitudes_m`, `angular_frequencies_rad_s`,
`directions_rad`, `phases_rad`. Amplitudes are nonnegative; angular frequencies
are positive. Empty vectors mean calm water. The stationary component phase is
`k*(x*cos(direction)+y*sin(direction)) - omega*t + phase`; solve supplies the
documented physical upstream finite train. Direction is propagation CCW from
world+x. Arrays are authoritative; no generator replay is required to solve.

## result.npz

| Member | Shape | Meaning |
|---|---|---|
| time_s | T | Motion times, strictly increasing from0 |
| position_m | T x3 | Absolute world CG position |
| orientation_wxyz | T x4 | Unit active body-to-world quaternion |
| hull_vertices_m | N x3 | Complete exported physical hull, body coordinates |
| hull_faces | F x3 | Exported physical connectivity |
| cg_m | 3 | CG in exported hull coordinates |
| surface_time_s | S | Surface times, same endpoints as motion |
| surface_vertices_m | sum(Ns) x3 | Packed world-space total water-air interface |
| surface_faces | sum(Fs) x3 | Packed triangles, indices LOCAL to each frame |
| surface_vertex_offsets | S+1 | Frame boundaries in packed vertices |
| surface_face_offsets | S+1 | Frame boundaries in packed faces |

Both offset vectors start at0, increase strictly, and end at their packed array
length. Counts and connectivity may change by frame. The format imposes no
height-graph or linear-superposition restriction on future nonlinear output.

The hull can be copied or remeshed but must preserve the simulated physical
shape. Pose mapping is `position + R @ (vertex - cg)`. The result's hull is the
authoritative geometry for rendering that result; editing boat.npz afterward
does not reinterpret it. Surface geometry is mandatory and describes the total
interface represented by the solver, not an incident-only cache. Export preserves
the producer's surface mesh; it does not cut or remesh triangles at the hull.
Surface portions inside the posed hull are invalid fluid and are implicitly
excluded by that hull. Finite stored values there are not physical predictions.
Consumers compute the physical-hull/surface intersection separately when they
need a waterline; no waterline array or validity metadata is stored.

This permits an uncut WEC-Sim sampling grid as well as a native solver's evolving
surface mesh. A native solver need not produce an incident/diffraction/radiation
decomposition, a height graph, or values inside the hull. Format
validation is not a physical-accuracy certificate or a restart capability.
